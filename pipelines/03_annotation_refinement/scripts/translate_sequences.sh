#!/usr/bin/env bash

#SBATCH --job-name=translate_sequences
#SBATCH --partition=c
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=translate_sequences_%j.out
#SBATCH --error=translate_sequences_%j.err
#SBATCH --mail-user=wagner.guzman@gmi.oeaw.ac.at
#SBATCH --mail-type=FAIL

# Translate nucleotide sequences to 6 frames using EMBOSS transeq

IN="$1"                 # Input FASTA (nucleotides)
OUT="$2"                # Output FASTA (6 frames)
TABLE="${3:-11}"        # Genetic Code: 11 (Archaea/Bacteria)

# Load modules

ml build-env/2020
ml emboss/6.6.0-foss-2018b

echo "=== Translating sequences to 6 frames ==="
echo "Input:  $IN"
echo "Output: $OUT"
echo "Table:  $TABLE (Archaea/Bacteria)"

transeq -sequence "$IN" -outseq "$OUT" -frame 6 -table "$TABLE" -auto

# Statistics
in_n=$(grep -c '^>' "$IN" || true)
out_n=$(grep -c '^>' "$OUT" || true)

echo ""
echo "=== Translation Summary ==="
echo "Input sequences:      $in_n"
echo "Translated sequences: $out_n (expected: $((in_n * 6)))"
echo "Output saved to:      $OUT"

