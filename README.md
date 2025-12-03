<p align="center">
  <img src="assets/metarefine_logo.png" width="300">
</p>

<h1 align="center">MetaRefine: A Modular Framework for Prokaryotic Genome Refinement</h1>

<p align="center">
<a href="https://github.com/wgiovanni13/MetaRefine/issues"><img src="https://img.shields.io/github/issues/wgiovanni13/MetaRefine"></a>
<a href="https://github.com/wgiovanni13/MetaRefine"><img src="https://img.shields.io/github/last-commit/wgiovanni13/MetaRefine"></a>
<a href="LICENSE"><img src="https://img.shields.io/github/license/wgiovanni13/MetaRefine"></a>
<img src="https://img.shields.io/badge/Nextflow-DSL2-blue">
<img src="https://img.shields.io/badge/Category-Genomics-success">
</p>

---

# Table of Contents

- [About the project](#about-the-project)
- [Biological & Conceptual Overview](#biological-and-conceptual-overview)
- [Lokiarchaeon as a case of study](#lokiarchaeon-as-a-case-of-study)
- [Pipeline architecture](#pipeline-architecture)
- [When should users run each pipeline?](#when-should-users-run-each-pipeline)
- [Dependencies & External Resources](#dependencies-and-external-resources)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Contact](#contact)
- [Acknowledgements](#acknowledgements)

---

## 1. About the project

*MetaRefine is a modular Nextflow-based suite for **Prokaryotic Genome Refinement** that integrates:

- ** Hybrid long and short redas metagenome assembly and iterative polishing (Pipeline 1) **
- ** Metatranscriptome assembly and based evidence for codign sequences (Pipeline 2) **
- ** Multi-tool structural and functional annotation refinement (Pipeline 3) **
- ** Tool validation in Lokiarchaeota genome: Polishing and Transcriptome validation (Pipeline 4) **

Although the framework is demonstrated on a **Lokiarchaeota B-35** genome, it is designed to be **generalizable to prokaryotic genomes** (Bacteria and Archaea) with a genetic code 11 and genetic code 34 to be incorporated soon, and applicable to:

- Reference Genomes Only (Usage of Pipeline 3) 
- Datasets where either ** DNA and RNA ** data are available (Usage of Pipeline 2 or 3)

---

## 2. Biological & Conceptual Overview

High-throughput sequencing has enabled reconstruction of prokaryotic genomes directly from complex environments. However:

- Assemblies can contain **local errors** (indels or mis-joins) that cause **frameshifts**.
- Gene callers often **disagree* on gene boundaries and even presence/absence of genes.
- **Metatranscriptomic evidence** is rarely used systematically to refine genome-level annotation. 

---

## 3. Lokiarchaeon as a case of study

It is used to demonstrate:

- How iterative polishing improves the genome compared to the published reference.
- How multi-tool structural and functional annotation helps to identify frameshifts, new functional genes that can be missed by using only one annotation tool, and to correct gene boundaries.
- How metatranscriptome helps **validate or challenge** predicted frameshifts and ultimately correct them. 

---

## 4. Pipeline architecture

### 4.1. Pipeline Architecture for users 

```
┌─────────────────────────────────────────────────────────────┐
│                  METAREFINE WORKFLOW SUITE                  │
└─────────────────────────────────────────────────────────────┘

Pipeline 1: Metagenome Assembly & Polishing
============================================
Raw Nanopore Reads → QC (NanoStat/NanoPlot)
           ↓
    Trimming (NanoFilt)
           ↓
    Assembly (MetaFlye)
           ↓
    Extract Longest Contig (Depending on the number of chromosomes)
           ↓
    Polish with NextPolish (Illumina short reads)
           ↓
    Iterative Pilon Polishing → Final Assembly
                                       ↓
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Pipeline 2: Metatranscriptome Assembly
=======================================
RNA-seq Paired-End Reads → QC (FastQC/MultiQC)
           ↓
    Trimming (fastp)
           ↓
    Align to Reference (Bowtie2)
           ↓
    Filter Properly Paired Reads (Samtools)
           ↓
    Merge BAMs → Extract FASTQ
           ↓
    De novo Assembly (rnaSPAdes) → Transcripts
           ↓
    Map to Genome (minimap2)
           ↓
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Pipeline 3: Annotation Refinement, Frameshift Detection and Gene Discovery  
==========================================================================
Polished Assembly → Multi-tool Annotation:
                    - Prokka
                    - Prodigal (-g 11 / -g 34 implemented soon)
                    - GeneMarkS-2
                    - Glimmer3
                    - RAST
           ↓
    Clean & Analyze GFFs
           ↓
    Range Categorization:
           ↓
    HHsuite-based structural support (UniClust30 + PDB70)
           ↓
    Start Codon Correction
           ↓
    Frameshift Detection (2 methods):
    1. BLAST all-vs-all (within genome)
    2. DIAMOND vs UniRef90 (functional validation)
           ↓
    New Functional Gene Discovery
           ↓
    Final Annotation (.gff)+ Report
```
### 4.2. Pipeline Architecture for Validation & Improvement Analysisi in Lokiarchaeon B-35

```
┌─────────────────────────────────────────────────────────────┐
│               VALIDATION WORKFLOW in LOKI B-35              │
└─────────────────────────────────────────────────────────────┘

Pipeline 4: Validation & Improvement Analysis
==========================================================================

Layer A: Genome Comparison

Prokka Annotation in Reference Genome & Polished Genome
                 ↓
Run bi-directional BLAST to map CDS between both genomes  
                 ↓
     Frameshift comparison (corrected)       

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Layer B: Transcriptome Based Validation

Integration Frameshift results & DIAMOND frameshift calls
                 ↓
Confirmation of conflict regions in transcriptome
                 ↓
 Fusion detection & gene improvement 
                 ↓
   Automated Genome correction 

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Layer C: Validation of Genome Correction

Deletion of conflict region & Reannotation | Coverage-based improvement (short reads)

---
```

## 5. When should users run each pipeline?

MetaRefine is modular: users do **not** need to run all pipelines. Typical use-cases:

### 5.1. Only a finished genome (no reads available)

You have a prokaryotic genome FASTA, but no raw sequencing data.

- **Recommended:** run **Pipeline 3** (Annotation Refinement, Frameshift Detection and Gene Discovery).
- **Not required:** Pipeline 1 (metagenome assembly) and Pipeline 2 (metatranscriptome) can be skipped.
- Goal: refine gene models, detect frameshifts, identify new functional genes supported by homology and structure (DIAMOND + HHsuite).


### 5.2. Metagenomic long + short reads

You have **Nanopore** and **Illumina short reads** from a metagenomic sample:

- **Recommended sequence:**
  1. Run **Pipeline 1** to:
     - assemble the genome from metagenomic long reads
     - extract the target contig(s)
     - polish with short reads
  2. Take the polished assembly and run **Pipeline 3** for refined annotation.
- Goal: obtain a high-quality, polished genome and a refined annotation starting from raw metagenomic reads.


### 5.3. Metatranscriptomic reads with a reference / polished genome

You have **RNA-seq short reads** (metatranscriptome) and a reference genome (either from literature or from Pipeline 1).

- **Recommended sequence:**
  1. Run **Pipeline 2** to:
     - QC + trim RNA-seq reads
     - map them to the reference genome,
     - assemble transcripts and map them back to the genome.
  2. Combine **Pipeline 2 outputs** with **Pipeline 3 results** to:
     - validate frameshifts and gene models
     - identify transcriptionally supported conflicts
  3. Run **Pipeline 4**, if you want a full polishing + transcriptome-based correction analysis.
- Goal: use expression evidence to support or challenge your annotation and identify candidate corrections.


### 5.4. Full MetaRefine workflow

You have:

- Metagenomic long + short reads
- Metatranscriptomic reads
- Interest in **quantitatively showing** how polishing and transcriptome-based corrections improve annotation

- **Recommended full workflow:**
  1. **Pipeline 1** – assemble and polish the genome.
  2. **Pipeline 2** – process and assemble the metatranscriptome; map transcripts to the genome.
  3. **Pipeline 3** – perform multi-tool annotation refinement and gene discovery on the polished genome.
  4. **Pipeline 4** – perform:
     - original vs polished vs corrected genome comparison,
     - transcriptome-based validation of frameshifts and conflicts,
     - reannotation and coverage-based improvement analysis (short reads).

---

## 6. Dependencies & External Resources

### 6.1 Core Dependencies

- Nextflow (≥ 21.x, DSL2 enabled)
- Java 8 or newer
- Python 3 (with pandas and biopython)
- R
- Conda/mamba or HPC modules

### 6.2 Bioinformatics tools

**Metagenome assembly & polishing (Pipeline 1)**

- sra-tools (prefetch, fasterq-dump)
- NanoStat, NanoPlot, NanoFilt (Nanopack)
- FastQC, MultiQC
- fastp
- flye / metaflye
- bwa, samtools
- NextPolish
- Pilon

**Metatranscriptome Assembly (Pipeline 2)**

- FastQC, MultiQC
- fastp
- bowtie2
- samtools
- rnaSPAdes
- minimap2

**Annotation refinement (Pipeline 3)**
- Prokka
- Prodigal (genetic code 11; 34 coming soon)
- GeneMarkS-2
- Glimmer3
- RAST (as precomputed input)
- bedtools
- EMBOSS (transeq)
- blast+
- DIAMOND
- HHsuite (hhblits, hhmake and hhsearch)

**Validation & improvement (Pipeline 4)**

Same tools as above, plus:
pysam (for transcriptome validation scripts)
additional Python utilities for genome correction and coverage analysis

### 6.3 External Databases

- UniRef90 DIAMOND database (~45 GB)
- UniClust30 (~86 GB)
- PDB70 (~56 GB)

---

## 7. Installation

```bash

# Clone repository (needs to be finished)

git clone https://github.com/YourUsername/METAREFINE.git
cd METAREFINE

# Install dependencies

conda env create -f envs/pipeline1.yml
conda env create -f envs/pipeline2.yml
conda env create -f envs/pipeline3.yml
conda env create -f envs/pipeline4.yml

```

## 8. Quick Start: Running each pipeline

### Run Individual Pipelines

#### Pipeline 1: Metagenome Assembly & Polishing

```bash
cd metagenome

nextflow run main.nf \
  --longReadID SRR22273710 \
  --shortReadIDs 'SRR22273278,SRR22273277,SRR22273276' \
  --outDir results \
  -profile <your_profile>
```

#### Pipeline 2: Metatranscriptome Assembly

```bash
cd metatranscriptome

nextflow run main.nf \
  --sraIDs 'SRRXXXXXX,SRRYYYYYY,SRRZZZZZZ' \
  --genomeDir /path/to/reference_genome_dir \
  --refGenome /path/to/reference_genome.fasta \
  --outDir /scratch/your_user/nextflow_pipeline_metatranscriptome \
  -profile <your_profile>

```

#### Pipeline 3: Annotation Refinement

```bash
cd lokib35

nextflow run main.nf \
  --genome polished_assembly_8.fasta \
  --outdir results \
  --gms2_dir /path/to/gms2 \
  --hh_uniref /path/to/UniRef30_db \
  --hh_pdb70 /path/to/pdb70_db \
  --diamond_db_dir /path/to/DIAMOND_db \
  --diamond_db_name uniref90 \
  --scripts_dir ./scripts \
  -profile <your_profile>
```

#### Pipeline 4: Validation & Improvement Analysis

```bash
cd tool_validation

nextflow run compare_genomes.nf \
  --workdir /scratch/<your_user>/nextflow_output \
  --outdir /scratch/<your_user>/nextflow_output/comparison_results \
  --genome_original /scratch/.../CP104013.1.fasta \
  --genome_polished /scratch/.../polished_assembly_8.fasta \
  --diamond_db_dir /groups/.../DIAMOND_db \
  --scripts_dir /scratch/.../nextflow_output/scripts \
  -profile <your_profile>
```
---

## 9. Contact

- **Email:** wagner.guzman@gmi.oeaw.ac.at

---

## 10. Acknowledgements

- **Lab:** Berger Lab and IT Resources

---
