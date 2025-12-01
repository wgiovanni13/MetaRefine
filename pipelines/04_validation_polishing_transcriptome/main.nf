#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ========================================
// PARAMETERS
// ========================================

// Working directory

params.workdir = "/scratch/Giovanni.Guzman/nextflow_output"
params.outdir = "${params.workdir}/comparison_results"

// Input genomes

params.genome_original = "${params.workdir}/CP104013.1.fasta"
params.genome_polished = "${params.workdir}/polished_assembly_8.fasta"

// DIAMOND database

params.diamond_db_dir = '/groups/berger/lab/Giovanni.Guzman/DIAMOND_db'
params.diamond_db_name = 'uniref90'

// Scripts directory

params.scripts_dir = "${params.workdir}/scripts"

// Analysis parameters

params.transtable = 11
params.max_frameshift_gap = 200
params.min_frameshift_query_coverage = 80.0
params.max_frameshift_subject_ratio = 1.2

// ========================================
// WORKFLOW
// ========================================

workflow {
    
    // Create channels for both genomes
    genome_original_ch = Channel.fromPath(params.genome_original)
        .map { file -> tuple("original", file) }
    
    genome_polished_ch = Channel.fromPath(params.genome_polished)
        .map { file -> tuple("polished", file) }
    
    // Combine both genomes
    genomes_ch = genome_original_ch.mix(genome_polished_ch)
    
    // Stage 1: Annotate with Prokka
    PROKKA_ANNOTATION(genomes_ch)
    
    // Stage 2: Extract CDS sequences (nucleotides)
    EXTRACT_CDS_SEQUENCES(PROKKA_ANNOTATION.out.gff_faa)
    
    // Stage 3: DIAMOND vs UniRef90
    diamond_db_ch = Channel.fromPath("${params.diamond_db_dir}/uniref90.dmnd")
        .ifEmpty { error "DIAMOND database not found at ${params.diamond_db_dir}/uniref90.dmnd" }
    
    SPLIT_PROTEINS(PROKKA_ANNOTATION.out.faa)
    
    // Fix: Properly handle the chunk files with genome_type
    chunks_with_type = SPLIT_PROTEINS.out.chunks
        .transpose()
        .map { genome_type, chunk_file -> 
            tuple(genome_type, chunk_file)
        }
        .combine(diamond_db_ch)
    
    DIAMOND_VS_UNIREF90(chunks_with_type)
    
    // Group DIAMOND results by genome type
    diamond_grouped = DIAMOND_VS_UNIREF90.out.diamond_result
        .groupTuple()
    
    MERGE_DIAMOND_RESULTS(diamond_grouped)
    
    // Stage 4: Detect frameshifts from DIAMOND
    DETECT_FRAMESHIFTS_DIAMOND(MERGE_DIAMOND_RESULTS.out.merged_diamond)
    
    // Stage 5: Summarize frameshifts
    SUMMARIZE_FRAMESHIFTS_DIAMOND(DETECT_FRAMESHIFTS_DIAMOND.out.frameshifts)

    // ========================================
    // STAGE 6: Gene-level comparison (CDS-based)
    // ========================================
    
    // Get CDS sequences (nucleotides)
    cds_original = EXTRACT_CDS_SEQUENCES.out.cds_fasta
        .filter { it[0] == 'original' }
        .map { it[1] }
    
    cds_polished = EXTRACT_CDS_SEQUENCES.out.cds_fasta
        .filter { it[0] == 'polished' }
        .map { it[1] }
    
    // Get GFF files
    gff_original = PROKKA_ANNOTATION.out.gff
        .filter { it[0] == 'original' }
        .map { it[1] }
    
    gff_polished = PROKKA_ANNOTATION.out.gff
        .filter { it[0] == 'polished' }
        .map { it[1] }
    
    // Step 1: Add genomic coordinates to CDS FASTA headers
    ADD_GENOMIC_COORDS(
        cds_original.map { tuple('original', it) }
            .mix(cds_polished.map { tuple('polished', it) }),
        gff_original.map { tuple('original', it) }
            .mix(gff_polished.map { tuple('polished', it) })
            .groupTuple()
            .map { type, gffs -> tuple(type, gffs[0]) }
    )
    
    // Split by genome type for database creation
    coords_original = ADD_GENOMIC_COORDS.out.cds_with_coords
        .filter { it[0] == 'original' }
        .map { it[1] }
    
    coords_polished = ADD_GENOMIC_COORDS.out.cds_with_coords
        .filter { it[0] == 'polished' }
        .map { it[1] }
    
    // Step 2: Create BLAST databases
    CREATE_BLAST_DB_ORIGINAL(coords_original, 'original')
    CREATE_BLAST_DB_POLISHED(coords_polished, 'polished')
    
    // Step 3: Run bidirectional BLASTs
    BLAST_ORIGINAL_TO_POLISHED(
        coords_original,
        CREATE_BLAST_DB_POLISHED.out.blast_db
    )
    
    BLAST_POLISHED_TO_ORIGINAL(
        coords_polished,
        CREATE_BLAST_DB_ORIGINAL.out.blast_db
    )
    
    // Step 4: Classify genes (identical vs different)
    CLASSIFY_GENES(
        coords_original,
        coords_polished,
        BLAST_ORIGINAL_TO_POLISHED.out.blast_result,
        BLAST_POLISHED_TO_ORIGINAL.out.blast_result
    )
    
    // Step 5: Detect frameshifts and generate quality report
    DETECT_FRAMESHIFTS_AND_QUALITY_REPORT(
        CLASSIFY_GENES.out.identical_genes,
        CLASSIFY_GENES.out.original_different,
        CLASSIFY_GENES.out.polished_different,
        BLAST_ORIGINAL_TO_POLISHED.out.blast_result,
        BLAST_POLISHED_TO_ORIGINAL.out.blast_result
    )

    // ========================================
    // STAGE 7: Transcriptome-based validation
    // ========================================
    //
    // Requiere que ya existan:
    //  - comparison_results/gene_mapping_analysis/original_different.txt
    //  - comparison_results/gene_mapping_analysis/polished_different.txt
    //  - comparison_results/frameshifts/original_frameshifts.tsv
    //  - comparison_results/frameshifts/polished_frameshifts.tsv
    //  - BAMs de RNA-seq en /scratch/Giovanni.Guzman/nextflow_output
    //
    // El script usa rutas absolutas, así que no necesita canales de entrada.
    Channel
        .value(true)
        .set { dummy_ch }

    TRANSCRIPTOME_VALIDATION(dummy_ch)

    // STAGE 8: Automated genome correction + fusion detection
    //
    // Depende de:
    //  - polished_exact_conflicts.bed
    //  - polished_validation_*.txt
    //  - polished_assembly_8.fasta
    //  - polished.gff de Prokka
    //
    // Todos esos los usa el script con rutas absolutas.

    AUTOMATED_GENOME_CORRECTION(TRANSCRIPTOME_VALIDATION.out.log)
}

// ========================================
// PROCESSES
// ========================================

process PROKKA_ANNOTATION {
    tag "${genome_type}"
    publishDir "${params.outdir}/prokka/${genome_type}", mode: 'copy'
    
    module 'build-env/2020'
    
    cpus 2
    memory '8 GB'
    time '2h'
    
    input:
    tuple val(genome_type), path(genome)
    
    output:
    tuple val(genome_type), path("${genome_type}/${genome_type}.gff"), emit: gff
    tuple val(genome_type), path("${genome_type}/${genome_type}.faa"), emit: faa
    tuple val(genome_type), path("${genome_type}/${genome_type}.gff"), path("${genome_type}/${genome_type}.faa"), emit: gff_faa
    tuple val(genome_type), path("${genome_type}/${genome_type}.txt"), emit: stats
    path "${genome_type}/*"
    
    script:
    """
    source activate prokka_env
    
    echo "=== Annotating ${genome_type} genome with PROKKA ==="
    
    prokka --outdir ${genome_type} \
        --prefix ${genome_type} \
        --cpus ${task.cpus} \
        --kingdom Archaea \
        --genus Lokiarchaea \
        --usegenus \
        --locustag LOKI_${genome_type} \
        --rfam \
        --metagenome \
        --force \
        ${genome}
    
    # Verify output
    if [ ! -s "${genome_type}/${genome_type}.gff" ]; then
        echo "ERROR: PROKKA GFF is empty or missing!"
        exit 1
    fi
    
    if [ ! -s "${genome_type}/${genome_type}.faa" ]; then
        echo "ERROR: PROKKA FAA is empty or missing!"
        exit 1
    fi
    
    echo "✓ PROKKA annotation completed for ${genome_type}"
    echo "  - CDS features: \$(grep -c 'CDS' ${genome_type}/${genome_type}.gff)"
    echo "  - Protein sequences: \$(grep -c '^>' ${genome_type}/${genome_type}.faa)"
    """
}

process EXTRACT_CDS_SEQUENCES {
    tag "${genome_type}"
    publishDir "${params.outdir}/cds_sequences/${genome_type}", mode: 'copy'
    
    module 'build-env/2020'
    module 'bedtools/2.27.1-foss-2018b'
    module 'emboss/6.6.0-foss-2018b'
    
    cpus 2
    memory '4 GB'
    time '1h'
    
    input:
    tuple val(genome_type), path(gff), path(faa)
    
    output:
    tuple val(genome_type), path("${genome_type}_cds.fasta"), emit: cds_fasta
    tuple val(genome_type), path("${genome_type}_cds_translated.faa"), emit: cds_translated
    
    script:
    """
    # Extract CDS coordinates from GFF
    grep -w "CDS" ${gff} | \
        awk 'BEGIN{OFS="\\t"} {
            split(\$9, attrs, ";");
            id = "";
            for (i in attrs) {
                if (attrs[i] ~ /^ID=/) {
                    gsub(/ID=/, "", attrs[i]);
                    id = attrs[i];
                    break;
                }
            }
            print \$1, \$4-1, \$5, id, ".", \$7
        }' > ${genome_type}_cds.bed
    
    # Get genome from GFF (it's at the bottom)
    awk '/^##FASTA/,0' ${gff} | grep -v "^##FASTA" > ${genome_type}_genome.fasta
    
    # Extract sequences
    bedtools getfasta -fi ${genome_type}_genome.fasta \
        -bed ${genome_type}_cds.bed \
        -s -name -fo ${genome_type}_cds.fasta
    
    # Translate to verify
    transeq -sequence ${genome_type}_cds.fasta \
        -outseq ${genome_type}_cds_translated.faa \
        -table ${params.transtable} \
        -trim yes
    
    NUM_CDS=\$(grep -c '^>' ${genome_type}_cds.fasta)
    echo "✓ Extracted \${NUM_CDS} CDS sequences for ${genome_type}"
    """
}

process SPLIT_PROTEINS {
    tag "${genome_type}"
    publishDir "${params.outdir}/diamond_input/${genome_type}", mode: 'copy'
    
    module 'build-env/2020'
    module 'python/3.7.4-gcccore-8.3.0'
    module 'biopython/1.75-foss-2019b-python-3.7.4'
    
    cpus 1
    memory '4 GB'
    time '1h'
    
    input:
    tuple val(genome_type), path(faa)
    
    output:
    tuple val(genome_type), path("${genome_type}_chunk_*.faa"), emit: chunks
    
    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    
    # Target: 50,000 aa per chunk for faster processing
    TARGET_AA = 50_000
    
    current_chunk = 1
    current_aa = 0
    current_seqs = []
    
    print(f"Processing ${genome_type}: ${faa}")
    
    with open("${faa}") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq_len = len(record.seq)
            
            if current_aa > 0 and current_aa + seq_len > TARGET_AA:
                out_file = f"${genome_type}_chunk_{current_chunk:04d}.faa"
                with open(out_file, "w") as out:
                    SeqIO.write(current_seqs, out, "fasta")
                
                print(f"Created {out_file}: {len(current_seqs)} seqs, {current_aa:,} aa")
                
                current_chunk += 1
                current_seqs = []
                current_aa = 0
            
            current_seqs.append(record)
            current_aa += seq_len
    
    # Write last chunk
    if current_seqs:
        out_file = f"${genome_type}_chunk_{current_chunk:04d}.faa"
        with open(out_file, "w") as out:
            SeqIO.write(current_seqs, out, "fasta")
        print(f"Created {out_file}: {len(current_seqs)} seqs, {current_aa:,} aa")
    
    print(f"Total chunks for ${genome_type}: {current_chunk}")
    """
}

process DIAMOND_VS_UNIREF90 {
    tag "${genome_type}_${chunk_file.simpleName}"
    publishDir "${params.outdir}/diamond_results/${genome_type}", mode: 'copy'
    
    stageInMode 'symlink'
    
    cpus 16
    memory '48 GB'
    time '48h'
    maxForks 20
    
    module 'build-env/f2022'
    module 'diamond/2.1.8-gcc-12.2.0'
    
    input:
    tuple val(genome_type), path(chunk_file), path(diamond_db)
    
    output:
    tuple val(genome_type), path("${chunk_file.simpleName}.diamond.tsv"), emit: diamond_result
    
    script:
    """
    set -euo pipefail
    
    [[ ! -s "${chunk_file}" ]] && { echo "ERROR: Empty chunk file"; exit 1; }
    [[ ! -f "${diamond_db}" ]] && { echo "ERROR: DIAMOND DB not found"; exit 2; }
    
    LOCALTMP="/scratch/\${SLURM_JOB_ID}/diamond_tmp"
    mkdir -p \${LOCALTMP}
    
    OUTFILE="${chunk_file.simpleName}.diamond.tsv"
    
    echo "=== DIAMOND search for ${genome_type}: ${chunk_file.simpleName} ==="
    NUM_SEQS=\$(grep -c '^>' "${chunk_file}")
    echo "Sequences: \${NUM_SEQS}"
    echo "Threads: ${task.cpus}"
    
    diamond blastp \
        --query "${chunk_file}" \
        --db "${diamond_db.simpleName}" \
        --out "\${OUTFILE}" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \
        --threads ${task.cpus} \
        --tmpdir \${LOCALTMP} \
        --evalue 1e-5 \
        --max-target-seqs 20 \
        --block-size 4 \
        --index-chunks 1 \
        --id 30 \
        --query-cover 40 \
        --compress 0
    
    rm -rf \${LOCALTMP}
    
    if [ -f "\${OUTFILE}" ]; then
        HITS=\$(wc -l < "\${OUTFILE}")
        echo "✓ Completed: \${HITS} hits found"
    else
        echo "✓ Completed: no hits"
        touch "\${OUTFILE}"
    fi
    """
}

process MERGE_DIAMOND_RESULTS {
    tag "${genome_type}"
    publishDir "${params.outdir}/diamond_merged", mode: 'copy'
    
    cpus 2
    memory '8 GB'
    time '2h'
    
    input:
    tuple val(genome_type), path(diamond_chunks)
    
    output:
    tuple val(genome_type), path("${genome_type}_all.diamond.tsv"), emit: merged_diamond
    
    script:
    """
    echo "=== Merging DIAMOND results for ${genome_type} ==="
    
    # Count input files
    NUM_FILES=\$(ls -1 *.diamond.tsv | wc -l)
    echo "Merging \${NUM_FILES} chunk files"
    
    # Concatenate all results
    cat *.diamond.tsv > ${genome_type}_all.diamond.tsv
    
    TOTAL_HITS=\$(wc -l < ${genome_type}_all.diamond.tsv)
    echo "✓ Total hits for ${genome_type}: \${TOTAL_HITS}"
    """
}

process DETECT_FRAMESHIFTS_DIAMOND {
    tag "${genome_type}"
    publishDir "${params.outdir}/frameshifts", mode: 'copy'
    
    module 'build-env/2020'
    module 'python/3.10.4-gcccore-11.3.0'
    
    cpus 2
    memory '8 GB'
    time '2h'
    
    input:
    tuple val(genome_type), path(diamond_results)
    
    output:
    tuple val(genome_type), path("${genome_type}_frameshifts.tsv"), emit: frameshifts
    tuple val(genome_type), path("${genome_type}_frameshift_detection.log"), emit: log
    
    script:
    """
    python3 ${params.scripts_dir}/detect_frameshifts.py \
        ${diamond_results} \
        ${genome_type}_frameshifts.tsv \
        > ${genome_type}_frameshift_detection.log 2>&1
    
    cat ${genome_type}_frameshift_detection.log
    
    if [ -f "${genome_type}_frameshifts.tsv" ]; then
        NUM_FS=\$(tail -n +2 ${genome_type}_frameshifts.tsv | wc -l)
        echo ""
        echo "✓ Frameshift candidates detected for ${genome_type}: \${NUM_FS}"
    fi
    """
}

process SUMMARIZE_FRAMESHIFTS_DIAMOND {
    tag "${genome_type}"
    publishDir "${params.outdir}/frameshifts", mode: 'copy'
    
    module 'build-env/2020'
    module 'python/3.10.4-gcccore-11.3.0'
    
    cpus 1
    memory '4 GB'
    time '1h'
    
    input:
    tuple val(genome_type), path(frameshifts_tsv)
    
    output:
    tuple val(genome_type), path("${genome_type}_frameshift_summary.txt"), emit: summary
    
    script:
    """
    python3 ${params.scripts_dir}/summarize_frameshifts.py \
        ${frameshifts_tsv} \
        ${genome_type}_frameshift_summary.txt
    
    echo ""
    echo "=== FRAMESHIFT SUMMARY FOR ${genome_type.toUpperCase()} ==="
    cat ${genome_type}_frameshift_summary.txt
    """
}

process ADD_GENOMIC_COORDS {
    tag "${genome_type}"
    publishDir "${params.outdir}/cds_with_coords", mode: 'copy'
    
    module 'build-env/2020'
    module 'python/3.10.4-gcccore-11.3.0'
    module 'biopython/1.75-foss-2019b-python-3.7.4'
    
    cpus 1
    memory '4 GB'
    time '1h'
    
    input:
    tuple val(genome_type), path(cds_fasta)
    tuple val(genome_type), path(gff)
    
    output:
    tuple val(genome_type), path("${genome_type}_cds_with_coords.fasta"), emit: cds_with_coords
    
    script:
    """
    echo "=== Adding genomic coordinates to ${genome_type} CDS ==="
    
    python3 ${params.scripts_dir}/add_genomic_coords_to_fasta.py \
        ${cds_fasta} \
        ${gff} \
        ${genome_type}_cds_with_coords.fasta
    
    NUM_SEQS=\$(grep -c '^>' ${genome_type}_cds_with_coords.fasta)
    echo "✓ Processed \${NUM_SEQS} CDS sequences with coordinates"
    """
}

process CREATE_BLAST_DB_ORIGINAL {
    tag "original_db"
    publishDir "${params.outdir}/blast_dbs", mode: 'copy'
    
    module 'build-env/2020'
    module 'blast+/2.8.1-foss-2018b'
    
    cpus 1
    memory '4 GB'
    time '1h'
    
    input:
    path cds_with_coords
    val genome_type
    
    output:
    path "original_cds_db.*", emit: blast_db
    
    script:
    """
    echo "=== Creating BLAST database for original genome ==="
    
    makeblastdb \
        -in ${cds_with_coords} \
        -dbtype nucl \
        -out original_cds_db \
        -parse_seqids
    
    echo "✓ Original BLAST database created"
    """
}

process CREATE_BLAST_DB_POLISHED {
    tag "polished_db"
    publishDir "${params.outdir}/blast_dbs", mode: 'copy'
    
    module 'build-env/2020'
    module 'blast+/2.8.1-foss-2018b'
    
    cpus 1
    memory '4 GB'
    time '1h'
    
    input:
    path cds_with_coords
    val genome_type
    
    output:
    path "polished_cds_db.*", emit: blast_db
    
    script:
    """
    echo "=== Creating BLAST database for polished genome ==="
    
    makeblastdb \
        -in ${cds_with_coords} \
        -dbtype nucl \
        -out polished_cds_db \
        -parse_seqids
    
    echo "✓ Polished BLAST database created"
    """
}

process BLAST_ORIGINAL_TO_POLISHED {
    tag "original_to_polished"
    publishDir "${params.outdir}/gene_mapping_analysis", mode: 'copy'
    
    module 'build-env/2020'
    module 'blast+/2.8.1-foss-2018b'
    
    cpus 4
    memory '16 GB'
    time '2h'
    
    input:
    path query
    path db_files
    
    output:
    path "original_to_polished.blastn", emit: blast_result
    
    script:
    """
    echo "=== Running BLAST: Original → Polished ==="
    
    blastn \
        -query ${query} \
        -db polished_cds_db \
        -out original_to_polished.blastn \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
        -evalue 1e-10 \
        -max_target_seqs 10 \
        -num_threads ${task.cpus}
    
    NUM_HITS=\$(wc -l < original_to_polished.blastn)
    echo "✓ Found \${NUM_HITS} hits"
    """
}

process BLAST_POLISHED_TO_ORIGINAL {
    tag "polished_to_original"
    publishDir "${params.outdir}/gene_mapping_analysis", mode: 'copy'
    
    module 'build-env/2020'
    module 'blast+/2.8.1-foss-2018b'
    
    cpus 4
    memory '16 GB'
    time '2h'
    
    input:
    path query
    path db_files
    
    output:
    path "polished_to_original.blastn", emit: blast_result
    
    script:
    """
    echo "=== Running BLAST: Polished → Original ==="
    
    blastn \
        -query ${query} \
        -db original_cds_db \
        -out polished_to_original.blastn \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
        -evalue 1e-10 \
        -max_target_seqs 10 \
        -num_threads ${task.cpus}
    
    NUM_HITS=\$(wc -l < polished_to_original.blastn)
    echo "✓ Found \${NUM_HITS} hits"
    """
}

process CLASSIFY_GENES {
    tag "gene_classification"
    publishDir "${params.outdir}/gene_mapping_analysis", mode: 'copy'
    
    module 'build-env/2020'
    module 'python/3.10.4-gcccore-11.3.0'
    
    cpus 2
    memory '8 GB'
    time '1h'
    
    input:
    path original_cds_coords
    path polished_cds_coords
    path blast_orig_to_pol
    path blast_pol_to_orig
    
    output:
    path "identical_genes.txt", emit: identical_genes
    path "original_different.txt", emit: original_different
    path "polished_different.txt", emit: polished_different
    path "classification_summary.log", emit: log
    
    script:
    """
    echo "=== Classifying genes: Identical vs Different ==="
    
    python3 ${params.scripts_dir}/classify_identical_genes.py \
        ${original_cds_coords} \
        ${polished_cds_coords} \
        ${blast_orig_to_pol} \
        ${blast_pol_to_orig} \
        > classification_summary.log 2>&1
    
    cat classification_summary.log
    
    if [ ! -f "identical_genes.txt" ]; then
        echo "ERROR: identical_genes.txt not created"
        exit 1
    fi
    
    echo "✓ Gene classification completed"
    """
}

process DETECT_FRAMESHIFTS_AND_QUALITY_REPORT {
    tag "frameshift_detection"
    publishDir "${params.outdir}/gene_mapping_analysis", mode: 'copy'
    
    module 'build-env/2020'
    module 'python/3.10.4-gcccore-11.3.0'
    
    cpus 2
    memory '16 GB'
    time '1h'
    
    input:
    path identical_genes
    path original_different
    path polished_different
    path blast_orig_to_pol
    path blast_pol_to_orig
    
    output:
    path "category1_frameshifts_corrected.txt", emit: cat1
    path "category2_frameshifts_introduced.txt", emit: cat2
    path "category3_both_frameshifts.txt",       emit: cat3
    path "category4_length_variants.txt",        emit: cat4
    path "category5_lost_in_polished.txt",       emit: cat5
    path "category6_new_in_polished.txt",        emit: cat6
    path "uncategorized_genes.txt",              emit: uncategorized
    path "QUALITY_REPORT.txt",                   emit: quality_report
    path "frameshift_analysis.log",              emit: log
    
    script:
    """
    echo "============================================================================"
    echo "COMPREHENSIVE GENE COMPARISON ANALYSIS"
    echo "============================================================================"
    
    echo "Archivos en el working dir:"
    ls -1
    
    echo ""
    echo "Calling detect_frameshifts.py with explicit arguments:"
    echo "  identical_genes      = ${identical_genes}"
    echo "  original_different   = ${original_different}"
    echo "  polished_different   = ${polished_different}"
    echo "  blast_orig_to_pol    = ${blast_orig_to_pol}"
    echo "  blast_pol_to_orig    = ${blast_pol_to_orig}"
    echo ""

    python3 ${params.scripts_dir}/detect_frameshifts.py \
        ${identical_genes} \
        ${original_different} \
        ${polished_different} \
        ${blast_orig_to_pol} \
        ${blast_pol_to_orig} \
        > frameshift_analysis.log 2>&1
    
    cat frameshift_analysis.log
    
    if [ -f "QUALITY_REPORT.txt" ]; then
        echo ""
        echo "============================================================================"
        echo "QUALITY REPORT"
        echo "============================================================================"
        cat QUALITY_REPORT.txt
    else
        echo "WARNING: QUALITY_REPORT.txt not generated"
    fi
    
    echo ""
    echo "✓ Analysis completed"
    """
}

process TRANSCRIPTOME_VALIDATION {
    tag "transcriptome_validation"
    // No hace falta publishDir porque el script escribe con rutas absolutas

    module 'build-env/2020'
    module 'python/3.6.6-foss-2018b'
    module 'pysam/0.15.1-foss-2018b-python-3.6.6'

    cpus 2
    memory '8 GB'
    time '2h'

    /*
     * No necesita input porque:
     * - Usa rutas absolutas a:
     *   - original_different.txt / polished_different.txt
     *   - original_frameshifts.tsv / polished_frameshifts.tsv
     *   - GFFs de Prokka
     *   - BAMs de RNA-seq
     * Todo eso está hard-coded dentro del script.
     */
    input:
    val dummy

    output:
    path "transcriptome_validation.log", emit: log

    script:
    """
    echo "============================================================================"
    echo "TRANSCRIPTOME VALIDATION OF CONFLICTING REGIONS"
    echo "============================================================================"

    echo "Running comprehensive_transcriptome_validation.py"
    echo "Script path: ${params.scripts_dir}/comprehensive_transcriptome_validation.py"
    echo ""

    python3 ${params.scripts_dir}/comprehensive_transcriptome_validation.py \
        > transcriptome_validation.log 2>&1

    echo ""
    echo "=== Transcriptome validation log ==="
    cat transcriptome_validation.log

    echo ""
    echo "✓ Transcriptome validation completed"
    """
}

process AUTOMATED_GENOME_CORRECTION {
    tag "genome_correction"
    // El script escribe todo en BASE_DIR (/scratch/.../nextflow_output),
    // así que aquí solo guardamos el log.
    publishDir "${params.workdir}", mode: 'copy'

    module 'build-env/2020'
    module 'python/3.10.4-gcccore-11.3.0'
    module 'biopython/1.75-foss-2019b-python-3.7.4'
    // Prokka se usa vía subprocess en el script → necesitamos el entorno:
    // ajusta esta línea a como lo cargas en PROKKA_ANNOTATION si es distinto.
    
    cpus 4
    memory '32 GB'
    time '24h'

    /*
     * Input solo para forzar el orden: queremos que esto vaya
     * DESPUÉS de TRANSCRIPTOME_VALIDATION, porque usa:
     *  - polished_exact_conflicts.bed
     *  - polished_validation_*.txt
     */
    input:
    path transcriptome_log

    output:
    path "genome_correction.log", emit: log

    script:
    """
    echo "============================================================================"
    echo "AUTOMATED GENOME CORRECTION & FUSION DETECTION"
    echo "============================================================================"

    # Activar entorno de Prokka si lo necesitas (igual que en PROKKA_ANNOTATION)
    source activate prokka_env

    echo "Running automated_genome_correction_and_fusion_detection.py"
    echo "Script path: ${params.scripts_dir}/automated_genome_correction_and_fusion_detection.py"
    echo ""

    python3 ${params.scripts_dir}/automated_genome_correction_and_fusion_detection.py \
        > genome_correction.log 2>&1

    echo ""
    echo "=== Genome correction log ==="
    cat genome_correction.log

    echo ""
    echo "✓ Automated genome correction completed"
    """
}
