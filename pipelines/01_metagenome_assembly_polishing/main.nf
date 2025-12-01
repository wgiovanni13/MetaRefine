#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ===========================================
//               PARAMETERS
// ===========================================

params.outDir = 'results'
params.threads = 1
params.memory = '4.GB'
params.longReadID = 'SRR22273710'
params.shortReadIDs = 'SRR22273278,SRR22273277,SRR22273276'
params.assemblyDir = 'assembly'


// ===========================================
//             INPUT VALIDATION
// ===========================================

short_read_ids = params.shortReadIDs.split(/[\s,]+/)*.trim().grep { it }

short_read_ids.each { id ->
if (!(id ==~ /^SRR\d+$/)) {
error("Invalid SRA ID format: '$id'")
}
}

// ===========================================
//                  WORKFLOW
// ===========================================

workflow {
    def long_read = Channel.value(params.longReadID)
    def short_reads = Channel.fromList(short_read_ids)

    def long_downloaded = download_long_reads(long_read)
    quality_control_long_reads(long_downloaded)
    def long_trimmed = trim_long_reads(long_downloaded)
    quality_control_trimmed_long_reads(long_trimmed.trimmed_reads)
    def assembly_results = metaflye_assembly(long_trimmed.trimmed_reads)
    def longest_contig = extract_longest_contig(assembly_results)

    def short_downloaded = download_short_reads(short_reads)
    quality_control_short_reads(short_downloaded)
    def short_trimmed = trim_short_reads(short_downloaded)
    quality_control_trimmed_reads(short_trimmed.trimmed_reads)

    def short_reads_for_polish = short_downloaded
        .map { id, files -> files }
        .collect()

    def polished = polish_with_nextpolish(longest_contig, short_reads_for_polish)
    def formatted = format_fasta(polished)

    def final_assembly = polish_with_pilon_iterative(formatted, short_reads_for_polish)
}

// ===========================================
//                PROCESSES
// ===========================================

// ===========================================
//            DOWNLOAD PROCESSES 
// ===========================================

process download_long_reads {
tag "$long_read_id"
publishDir "${params.outDir}/long_reads", mode: 'copy'

cpus params.threads
memory params.memory

input:
val long_read_id

output:
tuple val(long_read_id), path("${long_read_id}.fastq")

script:
"""
ml build-env/f2021
ml sra-toolkit/3.0.0-centos_linux64

prefetch --progress --verify yes -O . ${long_read_id}
fasterq-dump --threads ${task.cpus} --progress -O . ${long_read_id}

# Check generated file and move it to the correct name

if [ -f "${long_read_id}_1.fastq" ]; then
mv ${long_read_id}_1.fastq ${long_read_id}.fastq
elif [ -f "${long_read_id}.fastq" ]; then
echo "File now have correct name"
else
echo "Error: Fastq file not found"
exit 1
fi
"""
}

process download_short_reads {
tag "$short_read_id"
publishDir "${params.outDir}/short_reads", mode: 'copy'

cpus 1
memory '32.GB'

input:
val short_read_id

output:
tuple val(short_read_id), path("${short_read_id}_{1,2}.fastq")

script:
"""
ml sra-toolkit/3.0.0-centos_linux64

prefetch --progress --verify yes -O . ${short_read_id}
fasterq-dump --split-files --threads 1 --progress -O . ${short_read_id}
"""
}

// ===========================================
//            QC PROCESS LONG READS
// ===========================================

process quality_control_long_reads {
tag "$long_read_id"
publishDir "${params.outDir}/QC/long_reads", mode: 'copy'

cpus params.threads
memory params.memory

input:
tuple val(long_read_id), path(fastq_file)

output:
tuple path("*html"), path("${long_read_id}_nanostat_summary.txt")

script:
"""
source activate nanopack_env

NanoStat \
--fastq ${fastq_file} \
--outdir . \
--prefix ${long_read_id}_ \
--name ${long_read_id}_nanostat_summary.txt

NanoPlot \
--fastq ${fastq_file} \
--outdir . \
--prefix ${long_read_id}_ \
--plots dot
"""
}

// ===========================================
//         TRIMMING PROCESS LONG READS
// ===========================================

process trim_long_reads {
tag "$long_read_id"
publishDir "${params.outDir}/trimmed_long_reads", mode: 'copy'

cpus params.threads
memory params.memory

input:
tuple val(long_read_id), path(fastq_file)

output:
tuple val(long_read_id), path("${long_read_id}_trimmed.fastq"), emit: trimmed_reads
path "*_trimming_summary.txt", emit: trimming_summary

script:
"""
source activate nanopack_env

echo "Pre-trimming stats:" > ${long_read_id}_trimming_summary.txt
NanoStat --fastq ${fastq_file} >> ${long_read_id}_trimming_summary.txt

cat ${fastq_file} | \
NanoFilt \
--quality 10 \
--length 1000 \
--maxlength 150000 \
--headcrop 50 \
--tailcrop 50 \
> ${long_read_id}_trimmed.fastq

echo -e "\nPost-trimming stats:" >> ${long_read_id}_trimming_summary.txt
NanoStat --fastq ${long_read_id}_trimmed.fastq >> ${long_read_id}_trimming_summary.txt
"""
}

// ===========================================
//    QC AFTER TRIMMING PROCESS LONG READS
// ===========================================

process quality_control_trimmed_long_reads {
tag "$long_read_id"
publishDir "${params.outDir}/QC/trimmed_long_reads", mode: 'copy'

cpus params.threads
memory params.memory

input:
tuple val(long_read_id), path(fastq_file)

output:
path("*{html,txt}")

script:
"""
source activate nanopack_env

NanoStat \
--fastq ${fastq_file} \
--outdir . \
--prefix ${long_read_id}_trimmed_ \
--name ${long_read_id}_trimmed_nanostat_summary.txt

NanoPlot \
--fastq ${fastq_file} \
--outdir . \
--prefix ${long_read_id}_trimmed_ \
--plots dot
"""
}

// ===========================================
//             ASSEMBLY PROCESS
// ===========================================

process metaflye_assembly {
tag "$long_read_id"
publishDir "${params.outDir}/${params.assemblyDir}/metaflye", mode: 'copy'

clusterOptions = '--qos=medium'

cpus 16
memory '128.GB'
time '47h'

input:
tuple val(long_read_id), path(fastq_file)

output:
path("flye_output")

script:
"""
mkdir -p flye_output
source activate metaflye_env

flye \
--nano-raw ${fastq_file} \
--out-dir flye_output \
--meta \
--threads ${task.cpus}
"""
}

// ===========================================
//              QC SHORT READS 
// ===========================================

process quality_control_short_reads {
tag "$short_read_id"
publishDir "${params.outDir}/QC/short_reads", mode: 'copy'

cpus params.threads
memory params.memory

input:
tuple val(short_read_id), path(fastq_files)

output:
path("${short_read_id}_QC_report.html")

script:
"""
ml build-env/f2021
ml fastqc/0.11.9-java-11
ml multiqc/1.11-foss-2020b-python-3.8.6

fastqc --threads ${task.cpus} -o . ${fastq_files}
multiqc . -o .
cp multiqc_report.html ${short_read_id}_QC_report.html
"""
}

// ===========================================
//            TRIMMING SHORT READS
// ===========================================

process trim_short_reads {
tag "$short_read_id"
publishDir "${params.outDir}/trimmed_reads", mode: 'copy'

cpus params.threads
memory params.memory

input:
tuple val(short_read_id), path(fastq_files)

output:
tuple val(short_read_id), path("*trim.fastq"), emit: trimmed_reads
path "*_fastp.{json,html}", emit: fastp_reports

script:
def read1 = fastq_files[0]
def read2 = fastq_files[1]
"""
ml build-env/2020
ml fastp/0.20.1-gcc-8.2.0-2.31.1

fastp \
-i ${read1} \
-I ${read2} \
-o ${short_read_id}_1.trim.fastq \
-O ${short_read_id}_2.trim.fastq \
--trim_front1 3 \
--trim_front2 3 \
--cut_tail \
--cut_window_size 3 \
--cut_mean_quality 20 \
--length_required 35 \
--correction \
--low_complexity_filter \
--thread ${task.cpus} \
--html ${short_read_id}_fastp.html \
--json ${short_read_id}_fastp.json
"""
}

// ===========================================
//        QC AFTER TRIMMING SHORT READS
// ===========================================

process quality_control_trimmed_reads {
tag "$short_read_id"
publishDir "${params.outDir}/QC/trimmed_reads", mode: 'copy'

cpus params.threads
memory params.memory

input:
tuple val(short_read_id), path(fastq_files)

output:
path("${short_read_id}_trimmed_QC_report.html")

script:
"""
ml build-env/f2021
ml fastqc/0.11.9-java-11
ml multiqc/1.11-foss-2020b-python-3.8.6

fastqc --threads ${task.cpus} -o . ${fastq_files}
multiqc . -o .
cp multiqc_report.html ${short_read_id}_trimmed_QC_report.html
"""
}

// ===========================================
//       EXTRACT LONGEST CONTIG PROCESS
// ===========================================

process extract_longest_contig {
tag "$long_read_id"
publishDir "${params.outDir}/${params.assemblyDir}/longest_contig", mode: 'copy'

cpus 1
memory '4.GB'

input:
path(flye_output)

output:
path("longest_contig.fasta")

script:
"""
ml build-env/f2021
ml samtools/1.14-gcc-10.2.0

# Find largest contig in assembly_info.txt

longest_contig=\$(awk 'NR>1 {print \$1,\$2}' ${flye_output}/assembly_info.txt | sort -k2,2nr | head -n1 | cut -d' ' -f1)

# Extract longest contig using samtools

samtools faidx ${flye_output}/assembly.fasta \${longest_contig} > longest_contig.fasta
"""
}

// ===========================================
//	  LR CURATION PROCESS WITH SR
// ===========================================

process polish_with_nextpolish {
tag "nextpolish"
publishDir "${params.outDir}/polished_assembly", mode: 'copy'

clusterOptions = '--qos=long'

cpus 4
memory '64.GB'
time '6h'

input:
path(longest_contig)
path(short_reads_fastq)

output:
path("output.genome.nextpolish.fasta")

script:
"""
source activate nextpolish_env

mkdir -p polished

# Create sgs.fofn file with absolut routes of fastq files

ls -d \$PWD/${short_reads_fastq} > sgs.fofn

# Create config file

cat <<EOF > run.cfg
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = ${task.cpus}
multithread_jobs = ${task.cpus}
genome = \$PWD/${longest_contig}
genome_size = auto
workdir = ./polished
polish_options = -p ${task.cpus}

[sgs_option]
sgs_fofn = \$PWD/sgs.fofn
sgs_options = -max_depth 400 -bwa
EOF

# Execute NextPolish

nextPolish run.cfg

# Copy file

cp polished/genome.nextpolish.fasta output.genome.nextpolish.fasta
"""
}

// ===========================================
//       CONVERT IN FASTA FORMAT PROCESS
// ===========================================

process format_fasta {
tag "format_fasta"
publishDir "${params.outDir}/final_assembly", mode: 'copy'

input:
path(input_fasta)

output:
path("formatted_assembly.fasta")

script:
"""
awk '/^>/ {

# Process header

header=\$0
gsub(" ", "_", header)
print header
next
}
{
sequence=toupper(\$0)
for(i=1; i<=length(sequence); i+=60) {
print substr(sequence, i, 60)
}
}' ${input_fasta} > formatted_assembly.fasta
"""
}

// ===========================================
//        ITERATIVE POLISHING PROCESS
// ===========================================

process polish_with_pilon_iterative {
    tag "pilon_iterative"
    publishDir "${params.outDir}/pilon_polishing", mode: 'copy'

    clusterOptions = '--qos=long'

    cpus 4
    memory '32.GB'
    time '13d23h'

    input:
    path(initial_assembly)
    path(short_reads_fastq)

    output:
    path("final_polished_assembly.fasta")
    path("pilon_iterations/*")

    script:
    """
    mkdir -p pilon_iterations
    cp ${initial_assembly} pilon_iterations/current_assembly.fasta

    iteration=1
    changes=1

    while [ \$changes -gt 0 ]; do

        # Load modules

        ml build-env/2020
        ml bwa/0.7.17-foss-2018b
        ml samtools/1.10-gcc-8.3.0
        ml pilon/1.23-java-1.8

        bwa index pilon_iterations/current_assembly.fasta

        # Map reads and process BAM files

        for reads in \$(ls -1 *_1.fastq | sed 's/_1.fastq//'); do
            bwa mem -t ${task.cpus} pilon_iterations/current_assembly.fasta \${reads}_1.fastq \${reads}_2.fastq | \
            samtools sort -n -o \${reads}_querysorted.bam
            samtools fixmate -m \${reads}_querysorted.bam \${reads}_fixmate.bam
            samtools sort -o \${reads}_coordsorted.bam \${reads}_fixmate.bam
            samtools markdup \${reads}_coordsorted.bam \${reads}_markdup.bam
            samtools index \${reads}_markdup.bam
        done

        # Combine BAM files

        samtools merge combined.bam *_markdup.bam
        samtools index combined.bam

        # Execute Pilon

        java -Xmx${task.memory.toGiga()}G -jar \$EBROOTPILON/pilon.jar \
            --genome pilon_iterations/current_assembly.fasta \
            --bam combined.bam \
            --output pilon_iterations/polished_assembly_\${iteration} \
            --changes \
            --fix all \
            --threads ${task.cpus}

        # Count changes
        
        changes=\$(wc -l < pilon_iterations/polished_assembly_\${iteration}.changes || echo "0")
        echo "Iteration \${iteration} completed with \${changes} changes"
        
        if [ \$changes -gt 0 ]; then
            cp pilon_iterations/polished_assembly_\${iteration}.fasta pilon_iterations/current_assembly.fasta
            iteration=\$((iteration + 1))
            
            # Cleaning temporal files for next iteration

            rm -f *_querysorted.bam *_fixmate.bam *_coordsorted.bam *_markdup.bam* combined.bam*
        fi
    done

    # Copy final result

    cp pilon_iterations/current_assembly.fasta final_polished_assembly.fasta
    """
}
