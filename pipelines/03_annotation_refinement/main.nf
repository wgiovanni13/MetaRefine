#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ========================================
//   PATHS: Input genome and outputs
// ========================================

params.genome = "${projectDir}/polished_assembly_8.fasta"
params.outdir = "${projectDir}/results"

// ========================================
// RESOURCES: External tools and databases
// ========================================

params.gms2_dir = "${projectDir}/resources/gms2_linux_64"
params.hh_uniref = "/groups/berger/lab/Giovanni.Guzman/uniclust30/UniRef30_2023_02"
params.hh_pdb70 = "/groups/berger/lab/Giovanni.Guzman/PDB70_db/pdb70"
params.diamond_db_dir = '/groups/berger/lab/Giovanni.Guzman/DIAMOND_db'
params.diamond_db_name = 'uniref90'

// ========================================
//     SCRIPTS: All helper scripts
// ========================================

params.scripts_dir = "${projectDir}/scripts"

// Category 1: Range categorization

params.categorizer = "${params.scripts_dir}/categorize_by_ranges_unique.py"
params.extract_cat3 = "${params.scripts_dir}/extract_cat3_sequences.sh"

// Translation scripts

params.translate_script = "${params.scripts_dir}/translate_6frames.sh"
params.one_stop_script = "${params.scripts_dir}/filter_one_stop.sh"

// HHpred analysis

params.group_hhr_script = "${params.scripts_dir}/group_hhr_by_coords.py"
params.select_best_script = "${params.scripts_dir}/select_best_query_per_group.py"
params.analyze_singles_script = "${params.scripts_dir}/analyze_singles.py"
params.filter_singles_script = "${params.scripts_dir}/filter_singles.py"
params.merge_ranges_script = "${params.scripts_dir}/merge_final_ranges.py"

// Frame identification

params.identify_frames_script = "${params.scripts_dir}/identify_correct_frames.py"
params.update_frames_script = "${params.scripts_dir}/update_merged_ranges_with_frames.py"
params.extract_unknown_script = "${params.scripts_dir}/extract_unknown_frame_sequences.sh"
params.translate_unknown_script = "${params.scripts_dir}/translate_unknown_frames_6frames.sh"
params.extract_final_script = "${params.scripts_dir}/extract_final_sequences.sh"

// BLAST and frameshift detection

params.blast_script = "${params.scripts_dir}/blast_all_vs_all.sh"
params.detect_frameshifts_blast_script = "${params.scripts_dir}/detect_frameshifts_from_blast.py"
params.translate_final_script = "${params.scripts_dir}/translate_sequences.sh"
params.extract_frame_script = "${params.scripts_dir}/extract_correct_frame.py"

// Category 2: Frameshift detection with DIAMOND

params.detect_frameshifts_script = "${params.scripts_dir}/detect_frameshifts.py"
params.summarize_frameshifts_script = "${params.scripts_dir}/summarize_frameshifts.py"

// Category 3: New functional genes

params.find_new_genes_script = "${params.scripts_dir}/find_new_functional_genes.py"

// R scripts

params.clean_gff_script = "${params.scripts_dir}/clean_gff.R"
params.analyze_gff_script = "${params.scripts_dir}/analyze_gff.R"

// ========================================
//    PARAMETERS: Analysis thresholds
// ========================================

// General
params.transtable = 11
params.promote_thresh = 10.0

// HHpred
params.hh_cpus = 2
params.hh_mem = '64 GB'
params.hh_time = '2d'
params.hhblits_n = 3
params.hhblits_e = '1e-6'
params.hhsearch_e = '1e-6'
params.max_hits = 10

// BLAST
params.blast_evalue = '1e-5'
params.min_blast_identity = 70.0
params.min_blast_coverage = 50.0

// Frameshift detection
params.max_frameshift_gap = 200
params.min_frameshift_query_coverage = 80.0
params.max_frameshift_subject_ratio = 1.2

// Category 3: Tools to compare
params.tools_to_compare = "prokka,prodigal,gms2,glimmer,rast"


// ========================================
//               WORKFLOW
// ========================================

workflow {
    
    // Stage 1: Gene prediction

    genome_ch = Channel.fromPath(params.genome)
    
    GENEMARK(genome_ch)
    GLIMMER(genome_ch)
    PRODIGAL(genome_ch)
    PROKKA(genome_ch)
    
    // Stage 2: Clean and analyze GFFs

    CLEAN_GFF(
        GENEMARK.out.gff,
        GLIMMER.out.gff,
        PRODIGAL.out.fna,
        PROKKA.out.gff,
        Channel.fromPath("${params.outdir}/rast/loki_b35_polished_rast.gff"),
        Channel.fromPath(params.clean_gff_script)
    )
    
    ANALYZE_GFF(
        CLEAN_GFF.out.collect(),
        Channel.fromPath(params.analyze_gff_script)
    )
    
    // Stage 3: Categorize ranges

    CATEGORIZE_RANGES(
        ANALYZE_GFF.out,
        Channel.fromPath(params.categorizer)
    )
    
    // Stage 4: Process Category 3

    EXTRACT_CAT3_SEQUENCES(
        CATEGORIZE_RANGES.out.cat3,
        genome_ch,
        ANALYZE_GFF.out,
        Channel.fromPath(params.extract_cat3)
    )
    
    TRANSLATE_6FRAMES(
        EXTRACT_CAT3_SEQUENCES.out.fasta,
        Channel.fromPath(params.translate_script)
    )
    
    FILTER_ONE_STOP_FRAMES(
        TRANSLATE_6FRAMES.out.faa,
        Channel.fromPath(params.one_stop_script)
    )
    
    // Stage 5: HHpred analysis

    // Validate HHpred databases exist
    if (!file("${params.hh_uniref}_a3m.ffdata").exists()) {
        error "ERROR: HHpred UniRef database not found at: ${params.hh_uniref}"
    }
    if (!file("${params.hh_pdb70}_a3m.ffdata").exists()) {
        error "ERROR: HHpred PDB70 database not found at: ${params.hh_pdb70}"
    }

    SPLIT_FASTA_ONESTOP(
        FILTER_ONE_STOP_FRAMES.out.faa
    )

    HHPRED(
        SPLIT_FASTA_ONESTOP.out.split_files.flatten()
    )

    hhr_collected = HHPRED.out.hhr_files
        .flatten()
        .filter { it.name.endsWith('.hhr') }
        .collect()

    GROUP_HHR_FILES(
        Channel.fromPath(params.group_hhr_script),
        hhr_collected
    )
    
    SELECT_BEST_QUERY(
        GROUP_HHR_FILES.out.groups,
        Channel.fromPath(params.select_best_script)
    )
    
    ANALYZE_SINGLES(
        Channel.fromPath(params.analyze_singles_script),
        GROUP_HHR_FILES.out.singles,
        ANALYZE_GFF.out
    )
    
    FILTER_SINGLES(
        Channel.fromPath(params.filter_singles_script),
        GROUP_HHR_FILES.out.singles,
        ANALYZE_SINGLES.out.detailed,
        GROUP_HHR_FILES.out.groups
    )
    
    // Stage 6: Merge and correct frames

    MERGE_FINAL_RANGES(
        Channel.fromPath(params.merge_ranges_script),
        CATEGORIZE_RANGES.out.cat2,
        SELECT_BEST_QUERY.out.best_queries,
        FILTER_SINGLES.out.all_keep,
        ANALYZE_GFF.out
    )
    
    EXTRACT_UNKNOWN_FRAME_RANGES(
        MERGE_FINAL_RANGES.out.merged_bed,
        genome_ch,
        Channel.fromPath(params.extract_unknown_script)
    )
    
    TRANSLATE_UNKNOWN_FRAMES_6FRAMES(
        EXTRACT_UNKNOWN_FRAME_RANGES.out.fasta,
        Channel.fromPath(params.translate_unknown_script)
    )
    
    IDENTIFY_CORRECT_FRAMES(
        TRANSLATE_UNKNOWN_FRAMES_6FRAMES.out.faa,
        Channel.fromPath(params.identify_frames_script)
    )
    
    UPDATE_MERGED_RANGES_WITH_FRAMES(
        IDENTIFY_CORRECT_FRAMES.out.frames_tsv,
        MERGE_FINAL_RANGES.out.merged_tsv,
        MERGE_FINAL_RANGES.out.merged_bed,
        Channel.fromPath(params.update_frames_script)
    )
    
    // Stage 7: Extract final sequences

    EXTRACT_FINAL_SEQUENCES(
        UPDATE_MERGED_RANGES_WITH_FRAMES.out.updated_bed,
        genome_ch,
        Channel.fromPath(params.extract_final_script)
    )
    
    // Stage 8: Frameshift detection (Method 1: BLAST ALL GENES VS ALL GENES)

    BLAST_ALL_VS_ALL(
        EXTRACT_FINAL_SEQUENCES.out.fasta,
        Channel.fromPath(params.blast_script)
    )
    
    DETECT_FRAMESHIFTS(
        BLAST_ALL_VS_ALL.out.blast_filtered,
        Channel.fromPath(params.detect_frameshifts_blast_script)
    )
    
    // Stage 9: Translate to proteins

    TRANSLATE_FINAL_RANGES(
        EXTRACT_FINAL_SEQUENCES.out.fasta,
        Channel.fromPath(params.translate_final_script)
    )
    
    EXTRACT_FINAL_FRAME(
        TRANSLATE_FINAL_RANGES.out.all_frames,
        Channel.fromPath(params.extract_frame_script)
    )

    // Stage 10: Frameshift detection (Method 2: DIAMOND)
    
    // Validate or prepare DIAMOND database
    
    diamond_db_path = "${params.diamond_db_dir}/uniref90.dmnd"
    
    if (file(diamond_db_path).exists()) {
        println "✓ DIAMOND database found: ${diamond_db_path}"
        diamond_db_ch = Channel.fromPath(diamond_db_path)
    } else {
        println "DIAMOND database not found. It will be downloaded and prepared."
        println "  This is a one-time operation that may take 2-3 hours."
        PREPARE_DIAMOND_DB()
        diamond_db_ch = PREPARE_DIAMOND_DB.out.diamond_db
    }

    SPLIT_PROTEINS_BY_SIZE(
        EXTRACT_FINAL_FRAME.out.proteins
    )
    
    chunks_ch = SPLIT_PROTEINS_BY_SIZE.out.chunks.flatten()
    diamond_db_ch = Channel.fromPath("${params.diamond_db_dir}/uniref90.dmnd")
        .ifEmpty { error "DIAMOND database not found" }
    
    DIAMOND_VS_UNIREF90_OPTIMIZED(
        chunks_ch.combine(diamond_db_ch)
    )
    
    PREPARE_DIAMOND_RESULTS(
        DIAMOND_VS_UNIREF90_OPTIMIZED.out.diamond_result.collect()
    ) 

    DETECT_FRAMESHIFTS_PER_CHUNK(
        DIAMOND_VS_UNIREF90_OPTIMIZED.out.diamond_result
    )
    
    MERGE_FRAMESHIFT_RESULTS(
        DETECT_FRAMESHIFTS_PER_CHUNK.out.frameshift_chunk.collect()
    )
    
    // Stage 11: New functional genes discovery

    tools_ch = Channel.from(params.tools_to_compare.split(',').collect { it.trim() })

    diamond_results_ch = PREPARE_DIAMOND_RESULTS.out.diamond_tsv.collect()
    
    combined_ch = tools_ch
        .combine(ANALYZE_GFF.out)           
        .combine(diamond_results_ch)         
        .map { tool, table, diamond_files ->
            tuple(tool, table, diamond_files)
        }
    
    FIND_NEW_FUNCTIONAL_GENES(combined_ch)
    
    SUMMARIZE_ALL_NEW_GENES(
        FIND_NEW_FUNCTIONAL_GENES.out.new_genes.collect(),
        FIND_NEW_FUNCTIONAL_GENES.out.summary.collect()
    )
    
}

// ========================================
//                PROCESSES
// ========================================

// ========================================
//      STAGE 1: GENE PREDICTION
// ========================================

process GENEMARK {
    publishDir "${params.outdir}/gms2", mode: 'copy'

    module 'build-env/f2022'
    module 'perl/5.38.0-gcccore-13.2.0'

    input:
    path genome

    output:
    path "gms2_output.gff", emit: gff
    path "gms2_output.faa", emit: faa

    script:
    """
    ${params.gms2_dir}/gms2.pl --seq ${genome} \
        --output gms2_output.gff \
        --format gff \
        --faa gms2_output.faa \
        --genome-type archaea
    """
}

process GLIMMER {
    publishDir "${params.outdir}/glimmer3", mode: 'copy'

    module 'build-env/2020'
    module 'samtools/1.10-gcc-8.3.0'

    input:
    path genome

    output:
    path "loki_glimmer3.predict", emit: predict
    path "loki_glimmer3.gff", emit: gff

    script:
    """
    source activate glimmer3_env

    # Detect long ORFs
    long-orfs --no_header --max_olap 30 --trans_table 11 \
        --start_codons atg,gtg,ttg --stop_codons tag,tga,taa \
        --cutoff 1.15 ${genome} annotations.longorfs

    # Extract reads
    extract -t ${genome} annotations.longorfs > extracted_reads.train

    # Build ICM model
    build-icm --depth 7 --width 12 --period 3 --trans_table 11 \
        --stop_codons tag,tga,taa --reverse loki.icm < extracted_reads.train

    # Gene prediction
    glimmer3 -o30 -g110 -t30 -r ${genome} loki.icm loki_glimmer3

    # Convert to GFF
    awk '{
        if (\$1 ~ /^>/) {
            contig=\$1; gsub(">", "", contig);
        } else {
            start = \$2; end = \$3; strand = (\$4 ~ /-/) ? "-" : "+";
            print contig "\\tglimmer3\\tCDS\\t" start "\\t" end "\\t.\\t" strand "\\t.\\tID=" \$1 ";"
        }
    }' loki_glimmer3.predict > loki_glimmer3.gff
    """
}

process PRODIGAL {
    publishDir "${params.outdir}/prodigal", mode: 'copy'

    module 'build-env/2020'
    module 'prodigal/2.6.3-gcccore-7.3.0'

    input:
    path genome

    output:
    path "prodigal_output.gff3", emit: gff
    path "prodigal_output.faa", emit: faa
    path "prodigal_output.fna", emit: fna

    script:
    """
    prodigal -i ${genome} \
        -o prodigal_output.gff3 \
        -a prodigal_output.faa \
        -d prodigal_output.fna \
        -p single -g 11
    """
}

process PROKKA {
    publishDir "${params.outdir}/prokka", mode: 'copy'

    module 'build-env/2020'

    input:
    path genome

    output:
    path "Lokiarchaea/*"
    path "Lokiarchaea/Lokiarchaea.gff", emit: gff

    script:
    """
    source activate prokka_env
    
    prokka --outdir Lokiarchaea \
        --prefix Lokiarchaea \
        --cpus 1 \
        --kingdom Archaea \
        --genus Lokiarchaea \
        --usegenus \
        --locustag LOKI \
        --rfam \
        --metagenome \
        --force ${genome}
    
    if [ ! -s "Lokiarchaea/Lokiarchaea.gff" ]; then
        echo "ERROR: PROKKA GFF is empty or missing!"
        ls -lah Lokiarchaea/
        exit 1
    fi
    
    echo "✓ PROKKA GFF: \$(wc -l < Lokiarchaea/Lokiarchaea.gff) lines"
    """
}

// ========================================
//   STAGE 2: CLEAN AND ANALYZE GFFs
// ========================================

process CLEAN_GFF {
    publishDir "${params.outdir}/cleaned_gff", mode: 'copy'

    module 'build-env/f2022'
    module 'r/4.4.1-gfbf-2023b'

    input:
    path genemark_gff
    path glimmer_gff
    path prodigal_fna
    path prokka_gff
    path rast_gff
    path clean_script

    output:
    path "cleaned_*.gff"

    script:
    """
    Rscript ${clean_script} genemark ${genemark_gff}
    Rscript ${clean_script} glimmer ${glimmer_gff}
    Rscript ${clean_script} prodigal ${prodigal_fna}
    Rscript ${clean_script} prokka ${prokka_gff}
    Rscript ${clean_script} rast ${rast_gff}
    """
}

process ANALYZE_GFF {
    publishDir "${params.outdir}/analysis", mode: 'copy'
    
    memory '32 GB'

    module 'build-env/f2022'
    module 'r/4.4.1-gfbf-2023b'

    input:
    path cleaned_files
    path analyze_script

    output:
    path "final_overlap_results_combined.tsv"

    script:
    """
    Rscript ${analyze_script}
    """
}

// ========================================
//      STAGE 3: CATEGORIZE RANGES
// ========================================

process CATEGORIZE_RANGES {
    publishDir "${params.outdir}/categorized_ranges", mode: 'copy'

    module 'build-env/f2022'
    module 'python/3.11.5-gcccore-13.2.0'

    input:
    path combined_tsv
    path categorizer_script

    output:
    path "cat1_single_tool.tsv", emit: cat1
    path "cat2_full_consensus_5.tsv", emit: cat2
    path "cat3_partial_overlaps.tsv", emit: cat3
    path "cat4_other_full_or_uncategorized.tsv", emit: cat4
    path "category_summary.tsv", emit: summary
    path "all_unique_ranges.tsv", emit: all_unique

    script:
    """
    python3 ${categorizer_script} \
        --input ${combined_tsv} \
        --thresh ${params.promote_thresh} \
        --outdir .
    """
}

// ========================================
//      STAGE 4: PROCESS CATEGORY 3
// ========================================

process EXTRACT_CAT3_SEQUENCES {
    publishDir "${params.outdir}/cat3_extracted_sequences", mode: 'copy'

    module 'build-env/f2022'
    module 'bedtools/2.31.1-gcc-13.2.0'

    input:
    path cat3_tsv
    path genome_fa
    path combined_tsv
    path extract_script

    output:
    path "cat3_ranges.bed", emit: bed
    path "cat3_ranges_extracted.fasta", emit: fasta

    script:
    """
    bash ${extract_script} \
        ${genome_fa} \
        ${cat3_tsv} \
        cat3_ranges_extracted.fasta \
        ${combined_tsv}
    """
}

process TRANSLATE_6FRAMES {
    publishDir "${params.outdir}/cat3_translated", mode: 'copy'

    module 'build-env/2020'
    module 'emboss/6.6.0-foss-2018b'

    input:
    path fasta_in
    path translate_script

    output:
    path "cat3_ranges_extracted.6frames.archaea.faa", emit: faa

    script:
    """
    bash ${translate_script} \
        ${fasta_in} \
        cat3_ranges_extracted.6frames.archaea.faa \
        ${params.transtable}
    """
}

process FILTER_ONE_STOP_FRAMES {
    publishDir "${params.outdir}/cat3_onestop", mode: 'copy'

    module 'build-env/f2022'

    input:
    path faa6
    path filter_script

    output:
    path "ranges_with_exactly1_terminal_stop.txt", emit: list
    path "cat3_ranges_extracted.6frames.archaea.oneStopStripped.faa", emit: faa

    script:
    """
    bash ${filter_script} \
        ${faa6} \
        ranges_with_exactly1_terminal_stop.txt \
        cat3_ranges_extracted.6frames.archaea.oneStopStripped.faa
    """
}

// ========================================
//       STAGE 5: HHPRED ANALYSIS
// ========================================

process SPLIT_FASTA_ONESTOP {
    publishDir "${params.outdir}/hhpred_results/split", mode: 'copy'

    module 'build-env/2020'
    module 'seqkit/0.13.2'

    input:
    path faa_in

    output:
    path "split/*.faa", emit: split_files

    script:
    """
    mkdir -p split
    seqkit split -i -O split -s 1 ${faa_in}
    """
}

process HHPRED {
    publishDir "${params.outdir}/hhpred_results", mode: 'copy'

    cpus params.hh_cpus
    memory params.hh_mem
    time params.hh_time

    module 'build-env/2020'

    input:
    path fa

    output:
    path "a3m/*.a3m", optional: true, emit: a3m_files
    path "hhm/*.hhm", optional: true, emit: hhm_files
    path "hhr/*.hhr", optional: true, emit: hhr_files

    script:
    """
    set -euo pipefail
    
    source activate hhsuite_env 2>/dev/null || true
    
    [[ ! -s "${fa}" ]] && { echo "ERROR: Empty input file"; exit 1; }
    
    echo "=== HHPRED on \$(hostname) ==="
    echo "Input: ${fa}"
    echo "Start: \$(date)"
    
    # Usar DBs directamente (sin copia)
    DB_UNIREF="${params.hh_uniref}"
    DB_PDB70="${params.hh_pdb70}"
    
    # Validar que existan
    [ ! -f "\${DB_UNIREF}_a3m.ffdata" ] && {
        echo "ERROR: UniRef DB not found: \${DB_UNIREF}_a3m.ffdata"
        exit 2
    }
    
    [ ! -f "\${DB_PDB70}_a3m.ffdata" ] && {
        echo "ERROR: PDB70 DB not found: \${DB_PDB70}_a3m.ffdata"
        exit 2
    }
    
    # Preparar directorios
    mkdir -p a3m hhm hhr
    
    # Extraer ID de la secuencia
    HEADER=\$(grep '^>' "${fa}" | head -1 | sed 's/^>//' | cut -d' ' -f1)
    BASE=\$(echo "\${HEADER}" | tr '|/: ' '____')
    
    echo "Processing: \${HEADER}"
    echo "Base name: \${BASE}"
    
    START=\$(date +%s)
    
    # hhblits con optimizaciones de I/O
    echo "[1/3] Running hhblits (using ${params.hhblits_n} iterations)..."
    hhblits \
        -M first \
        -i "${fa}" \
        -d "\${DB_UNIREF}" \
        -oa3m "a3m/\${BASE}.a3m" \
        -n ${params.hhblits_n} \
        -e ${params.hhblits_e} \
        -cpu ${task.cpus} \
        -v 0 \
        2>&1 | grep -E "^(Done|Error|WARNING)" || true
    
    if [ ! -s "a3m/\${BASE}.a3m" ]; then
        echo "ERROR: hhblits failed or produced empty output"
        exit 3
    fi
    
    A3M_SEQS=\$(grep -c "^>" "a3m/\${BASE}.a3m" || echo 0)
    echo "  ✓ Generated MSA with \${A3M_SEQS} sequences"
    
    # hhmake
    echo "[2/3] Running hhmake..."
    hhmake \
        -v 0 \
        -i "a3m/\${BASE}.a3m" \
        -o "hhm/\${BASE}.hhm" \
        2>&1
    
    if [ ! -s "hhm/\${BASE}.hhm" ]; then
        echo "ERROR: hhmake failed or produced empty output"
        exit 4
    fi
    echo "  ✓ HMM profile created"
    
    # hhsearch con optimizaciones
    echo "[3/3] Running hhsearch (top 10 hits)..."
    hhsearch \
        -i "hhm/\${BASE}.hhm" \
        -d "\${DB_PDB70}" \
        -o "hhr/\${BASE}.hhr" \
        -cpu ${task.cpus} \
        -E ${params.hhsearch_e} \
        -p 20 \
        -Z 10 -z 1 \
        -b 10 -B 10 \
        -v 0 \
        2>&1 | grep -E "^(Done|Error|WARNING)" || true
    
    if [ ! -s "hhr/\${BASE}.hhr" ]; then
        echo "ERROR: hhsearch failed or produced empty output"
        exit 5
    fi
    
    HITS=\$(grep -c "^No " "hhr/\${BASE}.hhr" || echo 0)
    echo "  ✓ Found \${HITS} structural hits"
    
    END=\$(date +%s)
    ELAPSED=\$((END - START))
    echo ""
    echo "✓ Completed in \${ELAPSED}s (\$((ELAPSED / 60))m \$((ELAPSED % 60))s)"
    echo "End: \$(date)"
    """
}

process GROUP_HHR_FILES {
    publishDir "${params.outdir}/hhpred_grouped", mode: 'copy'

    module 'build-env/f2022'
    module 'python/3.11.5-gcccore-13.2.0'

    input:
    path group_script
    path 'hhr_files/*.hhr'

    output:
    path "groups.tsv", emit: groups
    path "singles.tsv", emit: singles
    path "groups_by_size_summary.tsv", emit: summary
    path "skipped_unexpected_format.txt", optional: true, emit: skipped

    script:
    """
    set -e
    
    NUM_HHR=\$(ls hhr_files/*.hhr 2>/dev/null | wc -l || echo 0)
    echo "=== GROUP_HHR_FILES Debug ==="
    echo "Number of HHR files received: \${NUM_HHR}"
    echo "Files:"
    ls -lh hhr_files/ || echo "Directory empty or not found"
    echo "============================="
    
    if [ "\${NUM_HHR}" -eq 0 ]; then
        echo "ERROR: No HHR files received from HHPRED process!"
        echo "This usually means HHPRED didn't produce any output."
        exit 1
    fi
    
    python3 ${group_script} \
        --hhrdir hhr_files \
        --outdir .
    
    echo "✓ Grouping completed successfully"
    """
}

process SELECT_BEST_QUERY {
    publishDir "${params.outdir}/hhpred_best_queries", mode: 'copy'
    
    module 'build-env/f2022'
    module 'python/3.11.5-gcccore-13.2.0'
    
    input:
    path groups_tsv
    path select_script
    
    output:
    path "best_query_per_group.tsv", emit: best_queries
    
    script:
    """
    python3 ${select_script} \
        --groups ${groups_tsv} \
        --hhrdir ${params.outdir}/hhpred_results/hhr \
        --outdir . \
        --max-hits ${params.max_hits}
    """
}

process ANALYZE_SINGLES {
    publishDir "${params.outdir}/hhpred_grouped", mode: 'copy'

    module 'build-env/2020'
    module 'python/3.6.6-foss-2018b'
    module 'pandas/0.25.3-foss-2018b-python-3.6.6'

    input:
    path analyze_script
    path singles_tsv
    path annotation_tsv

    output:
    path "singles_detailed_analysis.tsv", emit: detailed
    path "singles_analysis_summary.tsv", emit: summary

    script:
    """
    python3 ${analyze_script} \
        --singles ${singles_tsv} \
        --annotation ${annotation_tsv} \
        --outdir .
    """
}

process FILTER_SINGLES {
    publishDir "${params.outdir}/hhpred_filtered", mode: 'copy'

    module 'build-env/2020'
    module 'python/3.6.6-foss-2018b'
    module 'pandas/0.25.3-foss-2018b-python-3.6.6'

    input:
    path filter_script
    path singles_tsv
    path detailed_tsv
    path groups_tsv

    output:
    path "keep_five_tools.tsv", optional: true, emit: five_tools
    path "keep_overlap_with_five_tools.tsv", optional: true, emit: overlap_five
    path "keep_different_strand_overlaps.tsv", optional: true, emit: diff_strand
    path "keep_same_strand_low_overlap.tsv", optional: true, emit: same_strand_low
    path "keep_no_overlaps.tsv", optional: true, emit: no_overlaps
    path "keep_by_priority_resolution.tsv", optional: true, emit: priority_winners
    path "all_ranges_to_keep.tsv", emit: all_keep
    path "removed_covered_by_groups.tsv", optional: true, emit: removed_covered
    path "removed_asymmetric_overlap.tsv", optional: true, emit: removed_asymmetric
    path "removed_by_priority.tsv", optional: true, emit: removed_priority
    path "uncategorized_ranges.tsv", optional: true, emit: uncategorized
    path "filtering_summary.txt", emit: summary

    script:
    """
    python3 ${filter_script} \
        --singles ${singles_tsv} \
        --detailed ${detailed_tsv} \
        --groups ${groups_tsv} \
        --coverage-threshold 90.0 \
        --overlap-threshold 30.0 \
        --asymmetric-high 80.0 \
        --asymmetric-low 45.0 \
        --outdir .
    """
}

// ========================================
//   STAGE 6: MERGE AND CORRECT FRAMES
// ========================================

process MERGE_FINAL_RANGES {
    publishDir "${params.outdir}/final_merged_ranges", mode: 'copy'
    
    module 'build-env/2020'
    module 'python/3.6.6-foss-2018b'
    module 'pandas/0.25.3-foss-2018b-python-3.6.6'
    
    input:
    path merge_script
    path cat2_tsv
    path best_query_tsv
    path filtered_tsv
    path annotation_tsv
    
    output:
    path "merged_final_ranges.tsv", emit: merged_tsv
    path "merged_final_ranges.bed", emit: merged_bed
    path "merge_summary.txt", emit: summary
    
    script:
    """
    python3 ${merge_script} \
        --cat2 ${cat2_tsv} \
        --best-query ${best_query_tsv} \
        --filtered ${filtered_tsv} \
        --annotation ${annotation_tsv} \
        --chrom-name "CP104013.1" \
        --output merged_final_ranges.tsv \
        --bed-output merged_final_ranges.bed \
        > merge_summary.txt 2>&1
    """
}

process EXTRACT_UNKNOWN_FRAME_RANGES {
    publishDir "${params.outdir}/unknown_frame_sequences", mode: 'copy'
    
    module 'build-env/f2022'
    module 'bedtools/2.31.1-gcc-13.2.0'
    
    input:
    path merged_bed
    path genome_fa
    path extract_script
    
    output:
    path "unknown_frame_ranges.bed", emit: bed
    path "unknown_frame_ranges.fasta", emit: fasta
    path "unknown_frame_stats.txt", emit: stats
    
    script:
    """
    bash ${extract_script} \
        ${genome_fa} \
        ${merged_bed} \
        unknown_frame_ranges.fasta
    """
}

process TRANSLATE_UNKNOWN_FRAMES_6FRAMES {
    publishDir "${params.outdir}/unknown_frame_translated", mode: 'copy'
    
    module 'build-env/2020'
    module 'emboss/6.6.0-foss-2018b'
    
    input:
    path fasta_in
    path translate_script
    
    output:
    path "unknown_frame_ranges.6frames.archaea.faa", emit: faa
    path "translation_stats.txt", emit: stats
    
    script:
    """
    bash ${translate_script} \
        ${fasta_in} \
        unknown_frame_ranges.6frames.archaea.faa \
        ${params.transtable}
    """
}

process IDENTIFY_CORRECT_FRAMES {
    publishDir "${params.outdir}/frame_identification", mode: 'copy'
    
    module 'build-env/2020'
    module 'python/3.6.6-foss-2018b'
    module 'pandas/0.25.3-foss-2018b-python-3.6.6'
    
    input:
    path faa_6frames
    path identify_script
    
    output:
    path "identified_frames.tsv", emit: frames_tsv
    path "frame_identification.log", emit: log
    
    script:
    """
    python3 ${identify_script} \
        --faa ${faa_6frames} \
        --output identified_frames.tsv \
        > frame_identification.log 2>&1
    
    cat frame_identification.log
    """
}

process UPDATE_MERGED_RANGES_WITH_FRAMES {
    publishDir "${params.outdir}/final_merged_ranges_updated", mode: 'copy'
    
    module 'build-env/2020'
    module 'python/3.6.6-foss-2018b'
    module 'pandas/0.25.3-foss-2018b-python-3.6.6'
    
    input:
    path frames_tsv
    path merged_tsv
    path merged_bed
    path update_script
    
    output:
    path "merged_final_ranges_updated.tsv", emit: updated_tsv
    path "merged_final_ranges_updated.bed", emit: updated_bed
    path "update_summary.txt", emit: summary
    
    script:
    """
    python3 ${update_script} \
        --frames ${frames_tsv} \
        --input-tsv ${merged_tsv} \
        --input-bed ${merged_bed} \
        --output-tsv merged_final_ranges_updated.tsv \
        --output-bed merged_final_ranges_updated.bed \
        --chrom-name "CP104013.1" \
        > update_summary.txt 2>&1
    
    cat update_summary.txt
    """
}

// ========================================
//   STAGE 7: EXTRACT FINAL SEQUENCES
// ========================================

process EXTRACT_FINAL_SEQUENCES {
    publishDir "${params.outdir}/final_sequences", mode: 'copy'
    
    module 'build-env/f2022'
    module 'bedtools/2.31.1-gcc-13.2.0'
    
    input:
    path updated_bed
    path genome_fa
    path extract_script
    
    output:
    path "final_ranges.fasta", emit: fasta
    path "final_ranges_corrected.bed", emit: bed_corrected
    path "extraction_summary.txt", emit: summary
    
    script:
    """
    bash ${extract_script} \
        ${genome_fa} \
        ${updated_bed} \
        final_ranges.fasta \
        > extraction_summary.txt 2>&1
    
    cat extraction_summary.txt
    """
}

// ========================================
// STAGE 8: FRAMESHIFT DETECTION (BLAST)
// ========================================

process BLAST_ALL_VS_ALL {
    publishDir "${params.outdir}/blast_all_vs_all", mode: 'copy'
    
    module 'build-env/2020'
    module 'blast+/2.8.1-foss-2018b'
    
    cpus 4
    memory '4 GB'
    time '1h'
    
    input:
    path fasta
    path blast_script
    
    output:
    path "final_ranges_blast.blastn", emit: blast_full
    path "final_ranges_blast.blastn.filtered", emit: blast_filtered
    path "final_ranges_blast_db.*", emit: blast_db
    path "blast_summary.txt", emit: summary
    
    script:
    """
    bash ${blast_script} \
        ${fasta} \
        final_ranges_blast \
        ${params.blast_evalue} \
        > blast_summary.txt 2>&1
    
    cat blast_summary.txt
    """
}

process DETECT_FRAMESHIFTS {
    publishDir "${params.outdir}/frameshift_detection", mode: 'copy'
    
    module 'build-env/f2022'
    module 'python/3.11.5-gcccore-13.2.0'
    
    input:
    path blast_filtered
    path detect_script
    
    output:
    path "frameshift_candidates.txt", emit: detailed
    path "frameshift_candidates.tsv", emit: tsv
    path "frameshift_detection.log", emit: log
    
    script:
    """
    python3 ${detect_script} \
        --blast ${blast_filtered} \
        --output frameshift_candidates.txt \
        --output-tsv frameshift_candidates.tsv \
        --min-identity ${params.min_blast_identity} \
        --max-gap ${params.max_frameshift_gap} \
        --min-query-coverage ${params.min_frameshift_query_coverage} \
        --max-subject-ratio ${params.max_frameshift_subject_ratio} \
        > frameshift_detection.log 2>&1
    
    cat frameshift_detection.log
    """
}

// ========================================
//    STAGE 9: TRANSLATE TO PROTEINS
// ========================================

process TRANSLATE_FINAL_RANGES {
    publishDir "${params.outdir}/translation_final", mode: 'copy'
    
    module 'build-env/2020'
    module 'emboss/6.6.0-foss-2018b'
    
    input:
    path nucleotides_fasta
    path translate_script
    
    output:
    path "final_ranges_6frames.faa", emit: all_frames
    path "translation_final.log", emit: log
    
    script:
    """
    bash ${translate_script} \
        ${nucleotides_fasta} \
        final_ranges_6frames.faa \
        ${params.transtable} \
        > translation_final.log 2>&1
    
    cat translation_final.log
    """
}

process EXTRACT_FINAL_FRAME {
    publishDir "${params.outdir}/final_proteins", mode: 'copy'
    
    module 'build-env/f2022'
    module 'python/3.11.5-gcccore-13.2.0'
    
    input:
    path all_frames_fasta
    path extract_script
    
    output:
    path "final_proteins.faa", emit: proteins
    path "frame_selection_stats.tsv", emit: stats
    path "rejected_sequences.txt", emit: rejected
    path "frame_extraction.log", emit: log
    
    script:
    """
    python3 ${extract_script} \
        --input ${all_frames_fasta} \
        --output final_proteins.faa \
        --stats frame_selection_stats.tsv \
        --rejected rejected_sequences.txt \
        > frame_extraction.log 2>&1
    
    cat frame_extraction.log
    """
}

// ========================================
// STAGE 10: FRAMESHIFT DETECTION (DIAMOND)
// ========================================

process PREPARE_DIAMOND_DB {
    publishDir "${params.diamond_db_dir}", mode: 'copy', saveAs: { filename -> 
        filename.equals('uniref90.dmnd') ? filename : null 
    }

    module 'build-env/f2022'
    module 'diamond/2.1.8-gcc-12.2.0'

    cpus 2
    memory '16 GB'
    time '8h'

    output:
    path "uniref90.dmnd", emit: diamond_db
    path "diamond_db_prep.log", emit: log

    when:
    !file("${params.diamond_db_dir}/uniref90.dmnd").exists()

    script:
    """
    echo "=== Preparing DIAMOND UniRef90 Database ===" | tee diamond_db_prep.log
    echo "Start time: \$(date)" | tee -a diamond_db_prep.log
    echo "" | tee -a diamond_db_prep.log

    # Download UniRef90 FASTA
    echo "Downloading UniRef90 database..." | tee -a diamond_db_prep.log
    wget -c ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz 2>&1 | tee -a diamond_db_prep.log

    # Verify download
    if [ ! -f uniref90.fasta.gz ]; then
        echo "ERROR: Download failed" | tee -a diamond_db_prep.log
        exit 1
    fi

    # Convert to DIAMOND format
    echo "" | tee -a diamond_db_prep.log
    echo "Converting to DIAMOND format (this may take 1-2 hours)..." | tee -a diamond_db_prep.log
    diamond makedb \
        --in uniref90.fasta.gz \
        -d uniref90 \
        --threads ${task.cpus} \
        2>&1 | tee -a diamond_db_prep.log

    # Verify DB creation
    if [ ! -f uniref90.dmnd ]; then
        echo "ERROR: DIAMOND database creation failed" | tee -a diamond_db_prep.log
        exit 1
    fi

    # Clean up
    echo "" | tee -a diamond_db_prep.log
    echo "Cleaning up temporary files..." | tee -a diamond_db_prep.log
    rm uniref90.fasta.gz

    echo "" | tee -a diamond_db_prep.log
    echo "End time: \$(date)" | tee -a diamond_db_prep.log
    echo "Database size: \$(ls -lh uniref90.dmnd | awk '{print \$5}')" | tee -a diamond_db_prep.log
    echo "" | tee -a diamond_db_prep.log
    echo "✓ UniRef90 DIAMOND database ready at: uniref90.dmnd" | tee -a diamond_db_prep.log

    cat diamond_db_prep.log
    """
}

process SPLIT_PROTEINS_BY_SIZE {
    publishDir "${params.outdir}/diamond_vs_uniref90/chunks", mode: 'copy'

    module 'build-env/2020'
    module 'python/3.10.4-gcccore-11.3.0'

    input:
    path faa_in

    output:
    path "chunk_*.faa", emit: chunks

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    
    # Target: Number of aa per chunk
    TARGET_AA = 100_000
    
    current_chunk = 1
    current_aa = 0
    current_seqs = []
    
    print(f"Processing input file: ${faa_in}")
    
    with open("${faa_in}") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq_len = len(record.seq)
            
            # If adding this sequence exceeds the limit, write chunk
            if current_aa > 0 and current_aa + seq_len > TARGET_AA:
                out_file = f"chunk_{current_chunk:04d}.faa"
                with open(out_file, "w") as out:
                    SeqIO.write(current_seqs, out, "fasta")
                
                print(f"Created chunk {current_chunk}: {len(current_seqs)} seqs, {current_aa:,} aa")
                
                current_chunk += 1
                current_seqs = []
                current_aa = 0
            
            current_seqs.append(record)
            current_aa += seq_len
    
    # Write last chunk
    if current_seqs:
        out_file = f"chunk_{current_chunk:04d}.faa"
        with open(out_file, "w") as out:
            SeqIO.write(current_seqs, out, "fasta")
        print(f"Created chunk {current_chunk}: {len(current_seqs)} seqs, {current_aa:,} aa")
    
    print(f"\\nTotal chunks created: {current_chunk}")
    """
}

process DIAMOND_VS_UNIREF90_OPTIMIZED {
    publishDir "${params.outdir}/diamond_vs_uniref90/results", mode: 'copy'
    
    stageInMode 'symlink'

    cpus 16
    memory '48 GB'
    time '48h'
    maxForks 10

    module 'build-env/f2022'
    module 'diamond/2.1.8-gcc-12.2.0'

    input:
    tuple path(chunk_faa), path(diamond_db)

    output:
    path "*.diamond.gz", emit: diamond_result

    tag "${chunk_faa.baseName}"

    script:
    """
    set -euo pipefail

    # Validate inputs
    [[ ! -s "${chunk_faa}" ]] && { echo "ERROR: Empty chunk file"; exit 1; }
    [[ ! -f "${diamond_db}" ]] && { echo "ERROR: DIAMOND DB not found"; exit 2; }

    # Tmp local on node
    LOCALTMP="/scratch/\${SLURM_JOB_ID}/diamond_tmp"
    mkdir -p \${LOCALTMP}

    # Output name
    OUTFILE="${chunk_faa.baseName}.diamond"

    echo "=== Processing chunk: ${chunk_faa.baseName} ==="
    NUM_SEQS=\$(grep -c '^>' "${chunk_faa}")
    TOTAL_AA=\$(grep -v '^>' "${chunk_faa}" | tr -d '\\n' | wc -c)
    echo "Sequences: \${NUM_SEQS}"
    echo "Total amino acids: \${TOTAL_AA}"
    echo "Threads: ${task.cpus}"
    echo "Memory: ${task.memory}"
    echo ""

    # DIAMOND with optimized parameters
    diamond blastp \
        --query "${chunk_faa}" \
        --db "${diamond_db.baseName}" \
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
        --compress 1

    # Compress output if not already compressed
    if [ -f "\${OUTFILE}" ] && [ ! -f "\${OUTFILE}.gz" ]; then
        gzip "\${OUTFILE}"
    fi

    # Clean tmp local
    rm -rf \${LOCALTMP}

    # Statistics
    if [ -f "\${OUTFILE}.gz" ]; then
        HITS=\$(zcat "\${OUTFILE}.gz" | wc -l)
        echo ""
        echo "✓ Completed: \${HITS} hits found"
    else
        echo ""
        echo "✓ Completed: no hits"
        touch "\${OUTFILE}.gz"
    fi
    """
}

process PREPARE_DIAMOND_RESULTS {
    publishDir "${params.outdir}/diamond_vs_uniref90/results", mode: 'copy', pattern: "*.diamond.tsv"

    memory '4 GB'
    cpus 2
    time '2h'

    input:
    path diamond_gz_files

    output:
    path "*.diamond.tsv", emit: diamond_tsv

    script:
    """
    set -euo pipefail

    echo "=== Decompressing DIAMOND results ==="
    
    for gz_file in *.diamond.gz; do
        if [ -s "\${gz_file}" ]; then
            base=\$(basename "\${gz_file}" .gz)
            echo "Decompressing \${gz_file}..."
            gunzip -c "\${gz_file}" > "\${base}"
            
            LINES=\$(wc -l < "\${base}")
            echo "  \${base}: \${LINES} hits"
        else
            echo "WARNING: Empty file \${gz_file}"
        fi
    done

    TOTAL_FILES=\$(ls -1 *.diamond.tsv 2>/dev/null | wc -l)
    echo ""
    echo "✓ Prepared \${TOTAL_FILES} DIAMOND result files"
    """
}

process DETECT_FRAMESHIFTS_PER_CHUNK {
    publishDir "${params.outdir}/frameshifts/chunks", mode: 'copy', pattern: "*.tsv"

    memory '4 GB'
    cpus 1
    time '1h'

    module 'build-env/2020'
    module 'python/3.10.4-gcccore-11.3.0'

    input:
    path diamond_chunk

    output:
    path "${diamond_chunk.baseName}.frameshifts.tsv", emit: frameshift_chunk

    tag "${diamond_chunk.baseName}"

    script:
    """
    set -euo pipefail
    
    echo "=== Detecting frameshifts in ${diamond_chunk.baseName} ==="
    
    python3 ${params.detect_frameshifts_script} \
        ${diamond_chunk} \
        ${diamond_chunk.baseName}.frameshifts.tsv
    
    if [ -f "${diamond_chunk.baseName}.frameshifts.tsv" ]; then
        FRAMESHIFTS=\$(tail -n +2 ${diamond_chunk.baseName}.frameshifts.tsv | wc -l)
        echo "✓ Found \${FRAMESHIFTS} frameshift candidates"
    fi
    """
}

process MERGE_FRAMESHIFT_RESULTS {
    publishDir "${params.outdir}/frameshifts", mode: 'copy'

    memory '8 GB'
    cpus 2
    time '2h'

    module 'build-env/2020'
    module 'python/3.10.4-gcccore-11.3.0'

    input:
    path frameshift_chunks

    output:
    path "all_frameshifts.tsv", emit: frameshifts_merged
    path "frameshift_summary.txt", emit: summary

    script:
    """
    set -euo pipefail

    echo "=== Merging frameshift results ==="
    
    # Get first file
    FIRST_FILE=\$(ls *.frameshifts.tsv | head -n 1)
    
    # Get header from first file
    head -n 1 \${FIRST_FILE} > all_frameshifts.tsv

    # Append all data (skip headers)
    for file in *.frameshifts.tsv; do
        tail -n +2 "\${file}" >> all_frameshifts.tsv
    done

    NUM_FILES=\$(ls *.frameshifts.tsv | wc -l)
    echo "✓ Merged \${NUM_FILES} chunk files"
    
    # Generate summary
    python3 ${params.summarize_frameshifts_script} \
        all_frameshifts.tsv \
        frameshift_summary.txt
    
    echo ""
    echo "=== SUMMARY ==="
    cat frameshift_summary.txt
    """
}

// ========================================
// STAGE 11: NEW FUNCTIONAL GENES DISCOVERY
// ========================================

process FIND_NEW_FUNCTIONAL_GENES {
    tag "${tool_name}"
    publishDir "${params.outdir}/new_functional_genes", mode: 'copy'

    memory '8 GB'
    cpus 2
    time '2h'

    module 'build-env/2020'
    module 'python/3.10.4-gcccore-11.3.0'

    input:
    tuple val(tool_name), path(merged_table), path('diamond_results/*.diamond.tsv')

    output:
    path "new_genes_not_in_${tool_name}.tsv", emit: new_genes
    path "new_genes_not_in_${tool_name}.summary.txt", emit: summary

    script:
    """
    set -euo pipefail

    echo "=== Finding new functional genes not in ${tool_name.toUpperCase()} ==="
    
    # Debug: verify DIAMOND results are present
    echo "DEBUG: DIAMOND results directory contents:"
    ls -lh diamond_results/ || echo "WARNING: No diamond_results directory!"
    echo ""
    echo "DEBUG: Counting DIAMOND chunk files:"
    CHUNK_COUNT=\$(ls -1 diamond_results/*.diamond.tsv 2>/dev/null | wc -l)
    echo "Found \${CHUNK_COUNT} DIAMOND chunk files"
    
    if [ "\${CHUNK_COUNT}" -eq 0 ]; then
        echo "ERROR: No DIAMOND result files found!"
        ls -lah diamond_results/
        exit 1
    fi
    echo ""

    python3 ${params.find_new_genes_script} \
        ${merged_table} \
        diamond_results \
        new_genes_not_in_${tool_name}.tsv \
        ${tool_name} 2> new_genes_not_in_${tool_name}.summary.txt

    if [ -f "new_genes_not_in_${tool_name}.tsv" ]; then
        NEW_GENES=\$(tail -n +2 "new_genes_not_in_${tool_name}.tsv" | wc -l || true)
        echo ""
        echo "=== RESULTS for ${tool_name.toUpperCase()} ==="
        echo "New functional genes found: \${NEW_GENES}"

        if [ "\${NEW_GENES}" -gt 0 ]; then
            echo ""
            echo "Top 5 genes:"
            head -n 6 "new_genes_not_in_${tool_name}.tsv" | tail -n 5 | \
              awk -F'\\t' '{printf "  %s | %s | Identity: %s%% | %s\\n", \$1, \$4, \$6, \$5}' | \
              cut -c1-100
        fi
    fi

    echo "✓ Analysis complete for ${tool_name}"
    """
}

process SUMMARIZE_ALL_NEW_GENES {
    publishDir "${params.outdir}/new_functional_genes", mode: 'copy'

    memory '4 GB'
    cpus 1
    time '30m'

    module 'build-env/2020'
    module 'python/3.10.4-gcccore-11.3.0'

    input:
    path tsv_files
    path summary_files

    output:
    path "combined_new_genes_summary.txt", emit: combined_summary
    path "all_tools_comparison.tsv", emit: comparison

    script:
    """
    set -euo pipefail
    shopt -s nullglob

    echo "=== Generating combined summary ===" > combined_new_genes_summary.txt
    echo "" >> combined_new_genes_summary.txt

    # Debug: listing received files
    echo "DEBUG: TSV files received:" >&2
    ls -lh new_genes_not_in_*.tsv 2>/dev/null || echo "No TSV files found" >&2
    echo "" >&2
    echo "DEBUG: Summary files received:" >&2
    ls -lh new_genes_not_in_*.summary.txt 2>/dev/null || echo "No summary files found" >&2
    echo "" >&2

    # Counting by tool
    for file in new_genes_not_in_*.tsv; do
        [ -f "\${file}" ] || continue
        tool=\$(basename "\${file}" .tsv | sed 's/new_genes_not_in_//')
        count=\$(tail -n +2 "\${file}" 2>/dev/null | wc -l || echo "0")
        printf "%-15s : %5d new functional genes\\n" "\${tool^^}" "\${count}" >> combined_new_genes_summary.txt
    done

    echo "" >> combined_new_genes_summary.txt
    echo "=== Detailed summaries ===" >> combined_new_genes_summary.txt
    echo "" >> combined_new_genes_summary.txt

    # Append individual summaries
    for file in new_genes_not_in_*.summary.txt; do
        [ -f "\${file}" ] || continue
        tool=\$(basename "\${file}" .summary.txt | sed 's/new_genes_not_in_//')
        echo "===== \${tool^^} =====" >> combined_new_genes_summary.txt
        cat "\${file}" 2>/dev/null >> combined_new_genes_summary.txt || echo "Error reading \${file}" >> combined_new_genes_summary.txt
        echo "" >> combined_new_genes_summary.txt
    done

    # Comparative table
    echo -e "tool\tnew_functional_genes\tavg_identity\tavg_coverage" > all_tools_comparison.tsv

    for file in new_genes_not_in_*.tsv; do
        [ -f "\${file}" ] || continue
        tool=\$(basename "\${file}" .tsv | sed 's/new_genes_not_in_//')
        count=\$(tail -n +2 "\${file}" 2>/dev/null | wc -l || echo "0")

        if [ "\${count}" -gt 0 ]; then
            avg_identity=\$(tail -n +2 "\${file}" | awk -F'\\t' '{sum+=\$7; n++} END {if(n>0) printf "%.1f", sum/n; else print "0"}')
            avg_coverage=\$(tail -n +2 "\${file}" | awk -F'\\t' '{sum+=\$8; n++} END {if(n>0) printf "%.1f", sum/n; else print "0"}')
        else
            avg_identity="0.0"
            avg_coverage="0.0"
        fi

        echo -e "$tool\t$count\t$avg_identity\t$avg_coverage" >> all_tools_comparison.tsv
    done

    echo ""
    echo "=== COMPARISON TABLE ==="
    if [ -s all_tools_comparison.tsv ]; then
        column -t all_tools_comparison.tsv 2>/dev/null || cat all_tools_comparison.tsv
    else
        echo "WARNING: all_tools_comparison.tsv is empty!"
    fi

    echo ""
    echo "=== COMBINED SUMMARY ==="
    cat combined_new_genes_summary.txt
    
    echo ""
    echo "✓ Combined summary generated successfully"
    """
}
