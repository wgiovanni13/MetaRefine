#!/usr/bin/env bash

#SBATCH --job-name=translate_6frames
#SBATCH --partition=c
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=translate6_%j.out
#SBATCH --error=translate6_%j.err
#SBATCH --mail-user=wagner.guzman@gmi.oeaw.ac.at
#SBATCH --mail-type=FAIL

# Translate sequences to 6 reading frames using EMBOSS transeq
# Args: Input FASTA, Output FASTA, Genetic Code Table

IN="$1"                 # unknown_frame_ranges.fasta
OUT="$2"                # unknown_frame_ranges.6frames.archaea.faa
TABLE="${3:-11}"        # Genetic Code: 11 (Bacteria/Archaea)

# Needed modules (already loaded by Nextflow process)
# ml build-env/2020
# ml emboss/6.6.0-foss-2018b

# Translate to 6 frames
transeq -sequence "$IN" -outseq "$OUT" -frame 6 -table "$TABLE" -auto

# Generate statistics
in_n=$(grep -c '^>' "$IN" || true)
out_n=$(grep -c '^>' "$OUT" || true)
expected=$((in_n * 6))

echo "Translation completed" > translation_stats.txt
echo "Input sequences        : $in_n" >> translation_stats.txt
echo "Translated sequences   : $out_n" >> translation_stats.txt
echo "Expected (6 frames)    : $expected" >> translation_stats.txt
echo "Output                 : $OUT" >> translation_stats.txt

# Print to stdout as well
cat translation_stats.txt

