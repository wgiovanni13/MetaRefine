#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ===========================================
//               PARAMETERS
// ===========================================

params.outDir = '/scratch/Giovanni.Guzman/nextflow_pipeline_metatranscriptome/'
params.threads = 6
params.memory = '8.GB'
params.genomeDir = '/users/wagner.guzman/scripts/loki_b36/frameshift'
params.refGenome = '/users/wagner.guzman/scripts/loki_b36/frameshift/polished_assembly_8.fasta'

if (!params.sraIDs) {
    error "You must provide --sraIDs as a comma- or space-separated list (e.g., --sraIDs 'SRR123,SRR456')"
}

sra_ids = params.sraIDs.split(/[\s,]+/)*.trim().grep { it }

sra_ids.each { id ->
    if (!(id ==~ /^SRR\d+$/)) {
        error("Invalid SRA ID format: '$id'")
    }
}

// ===========================================
//                WORKFLOW
// ===========================================

workflow {
    ch_sra_ids = Channel.fromList(sra_ids)

    ch_downloaded = ch_sra_ids | download_sra
    ch_qc = ch_downloaded | quality_control
    ch_trim = ch_downloaded | trim_reads
    ch_trimmed = ch_trim[0]
    ch_reports = ch_trim[1]

    ch_genome = Channel.of(params.genomeDir)
    ch_index = ch_genome | build_index

    ch_aln_input = ch_trimmed.combine(ch_index)
    ch_sam = ch_aln_input | align_reads
    ch_bam = ch_sam | sam_to_bam
    ch_merged = ch_bam.map { it[1] }.collect() | merge_bams
    ch_filtered = ch_merged | filter_properly_paired

    ch_fastq = ch_filtered.name_sorted_bam | bam_to_fastq
    ch_assembly = ch_fastq | run_rnaspades

    ch_ref_genome = Channel.fromPath(params.refGenome)
    map_assembly_to_reference(ch_assembly.transcripts, ch_ref_genome)
}

// ===========================================
//                PROCESSES 
// ===========================================

// ===========================================
//            DOWNLOAD PROCESSES 
// ===========================================

process download_sra {
tag "$sra_id"
publishDir "${params.outDir}/data", mode: 'copy'

cpus params.threads
memory params.memory

input:
val sra_id

output:
tuple val(sra_id), path("${sra_id}_*.fastq.gz")

script:
"""
ml sra-toolkit/3.0.0-centos_linux64

prefetch --progress --verify yes ${sra_id}
fasterq-dump --split-files --threads ${task.cpus} --progress -O . ${sra_id}
gzip ${sra_id}_*.fastq
"""
}

// ===========================================
//                QC PROCESS
// ===========================================

process quality_control {
tag "$sra_id"
publishDir "${params.outDir}/QC/${sra_id}", mode: 'copy', pattern: '*'

cpus params.threads
memory params.memory

input:
tuple val(sra_id), path(fastq_files)

output:
path("${sra_id}_QC_report.html")
path("multiqc_data")

script:
"""
ml fastqc/0.11.9-java-11
ml multiqc/1.11-foss-2020b-python-3.8.6

fastqc --threads ${task.cpus} -o . ${fastq_files}
multiqc . -o .
mv multiqc_report.html ${sra_id}_QC_report.html
"""
}

// ===========================================
//              TRIMMING PROCESS  
// ===========================================

process trim_reads {
tag "$sra_id"
publishDir "${params.outDir}/trimmed", mode: 'copy'

cpus params.threads
memory params.memory

input:
tuple val(sra_id), path(fastq_files)

output:
tuple val(sra_id), path("${sra_id}_1.trim.fastq.gz"), path("${sra_id}_2.trim.fastq.gz")
path("${sra_id}_fastp.*")

script:
def fq1 = fastq_files.find { it.name.endsWith('_1.fastq.gz') }
def fq2 = fastq_files.find { it.name.endsWith('_2.fastq.gz') }

"""
ml build-env/f2022
ml fastp/0.23.4-gcc-12.2.0

fastp -i $fq1 -I $fq2 \
--trim_front1 5 \
--trim_front2 5 \
--cut_tail \
--cut_window_size 4 \
--cut_mean_quality 20 \
--length_required 50 \
--thread ${task.cpus} \
--html ${sra_id}_fastp.html \
--json ${sra_id}_fastp.json \
-o ${sra_id}_1.trim.fastq \
-O ${sra_id}_2.trim.fastq

gzip ${sra_id}_1.trim.fastq
gzip ${sra_id}_2.trim.fastq
"""
}

// ===========================================
//             BUILD INDEX PROCESS
// ===========================================

process build_index {
publishDir "${params.outDir}/aligned/loki_index", mode: 'copy'

cpus params.threads
memory params.memory

input:
val genomeDir

output:
tuple val("index"), path("loki.*")

script:
"""
ml build-env/f2022
ml bowtie2/2.5.1-gcc-12.2.0

mkdir -p "${params.outDir}/aligned/loki_index"

bowtie2-build "${genomeDir}/polished_assembly.fasta" "loki"
"""
}

// ===========================================
//            ALIGN READS PROCESS
// ===========================================

process align_reads {
tag "$sra_id"
publishDir "${params.outDir}/aligned", mode: 'copy'

cpus params.threads
memory params.memory

input:
tuple val(sra_id), path(read1), path(read2), val("index"), path(index_files)

output:
path("${sra_id}.sam")

script:
"""
ml build-env/f2022
ml bowtie2/2.5.1-gcc-12.2.0

bowtie2 -x loki \
-1 $read1 \
-2 $read2 \
-p ${task.cpus} \
--very-sensitive \
--no-unal \
-S ${sra_id}.sam
"""
}

// ===========================================
//             SAM2BAM PROCESS
// ===========================================

process sam_to_bam {
tag "$sam_file.baseName"
publishDir "${params.outDir}/bam", mode: 'copy'

cpus params.threads
memory params.memory

input:
path sam_file

output:
tuple val(sam_file.baseName), path("${sam_file.baseName}.sorted.bam"), path("${sam_file.baseName}.sorted.bam.bai")

script:
"""
ml build-env/f2022
ml samtools/1.18-gcc-12.3.0

samtools view -@ ${task.cpus} -bS ${sam_file} | samtools sort -@ ${task.cpus} -o ${sam_file.baseName}.sorted.bam
samtools index ${sam_file.baseName}.sorted.bam
"""
}

// ===========================================
//             MERGE BAM PROCESS
// ===========================================

process merge_bams {
publishDir "${params.outDir}/merged", mode: 'copy'

cpus params.threads
memory '16.GB'

input:
path bams

output:
path("merged_6samples.sorted.bam")

script:
"""
ml build-env/f2022
ml samtools/1.18-gcc-12.3.0

samtools merge -@ ${task.cpus} merged_6samples.sorted.bam ${bams.join(' ')}
"""
}

// ===========================================
//         FILTER PAIRED PROCESS
// ===========================================

process filter_properly_paired {
    publishDir "${params.outDir}/merged", mode: 'copy'

    cpus params.threads
    memory '16.GB'

    input:
    path merged_bam

    output:
    path "merged_properly_paired.sorted.bam"
    path "merged_properly_paired.sorted.bam.bai"
    path "merged_properly_paired.name_sorted.bam", emit: name_sorted_bam

    script:
    """
    ml build-env/f2022
    ml samtools/1.18-gcc-12.3.0

    echo "Filtering properly paired reads..."
    samtools view -@ ${task.cpus} -h -b -f 2 -F 4 -F 256 ${merged_bam} > merged_properly_paired.bam

    echo "Sorting by coordinate..."
    samtools sort -@ ${task.cpus} -o merged_properly_paired.sorted.bam merged_properly_paired.bam
    samtools index merged_properly_paired.sorted.bam

    echo "Sorting by name for FASTQ extraction..."
    samtools sort -@ ${task.cpus} -n -o merged_properly_paired.name_sorted.bam merged_properly_paired.bam

    echo "Filtering and sorting completed."
    """
}

// ===========================================
//             BAM2FASTQ PROCESS
// ===========================================

process bam_to_fastq {
    publishDir "${params.outDir}/fastq_assembly", mode: 'copy'

    cpus params.threads
    memory '16.GB'

    input:
    path name_sorted_bam

    output:
    tuple path("reads_R1.fastq"), path("reads_R2.fastq")

    script:
    """
    ml build-env/f2022
    ml samtools/1.18-gcc-12.3.0

    samtools fastq -@ ${task.cpus} \
        -1 reads_R1.fastq \
        -2 reads_R2.fastq \
        -0 /dev/null -s /dev/null -n \
        ${name_sorted_bam}
    """
}

// ===========================================
//          ASSEMBLE TRANSCRIPT PROCESS
// ===========================================

process run_rnaspades {
    publishDir "${params.outDir}/rnaspades_assembly", mode: 'copy'
    errorStrategy 'ignore'
    
    cpus 8
    memory '96.GB'

    input:
    tuple path(reads_1), path(reads_2)

    output:
    path "spades_output/transcripts.fasta", emit: transcripts 
    path "spades_output"

    script:
    """
    ml build-env/f2021
    ml spades/3.15.3-gcc-10.3.0

    rnaspades.py \
        -1 ${reads_1} \
        -2 ${reads_2} \
        -t ${task.cpus} \
        -m 96 \
        -o spades_output
    """
}

// ===========================================
//       MAP CONTIG TO REFERENCE PROCESS
// ===========================================

process map_assembly_to_reference {
    publishDir "${params.outDir}/mapping_assembly", mode: 'copy'

    cpus 4
    memory '8.GB'

    input:
    path transcripts
    path ref_genome

    output:
    path "assembly_mapping.sorted.bam"
    path "assembly_mapping.sorted.bam.bai"

    script:
    """
    ml build-env/f2022
    ml minimap2/2.24-gcccore-11.2.0
    ml samtools/1.18-gcc-12.3.0

    echo "Mapping transcript assembly to reference genome using minimap2..."
    
    minimap2 -ax asm5 ${ref_genome} ${transcripts} > assembly_mapping.sam

    echo "Mapping completed. Converting SAM to sorted BAM..."

    samtools view -@ ${task.cpus} -bS assembly_mapping.sam | \
    samtools sort -@ ${task.cpus} -o assembly_mapping.sorted.bam
    
    samtools index assembly_mapping.sorted.bam

    echo "Assembly mapping and indexing completed."
    """
}
