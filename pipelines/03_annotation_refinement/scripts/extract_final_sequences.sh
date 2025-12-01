#!/usr/bin/env bash

#SBATCH --job-name=extract_final_sequences
#SBATCH --partition=c
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=extract_final_sequences_%j.out
#SBATCH --error=extract_final_sequences_%j.err
#SBATCH --mail-user=wagner.guzman@gmi.oeaw.ac.at
#SBATCH --mail-type=FAIL

# Extract sequences from final merged ranges BED file
# Args: Genome (.fasta), BED file, Output (.fasta)

GENOME="$1"
BED="$2"
OUT="$3"

# Extract the actual contig name from FASTA

CONTIG=$(awk '/^>/{sub(/^>/,""); print $1; exit}' "$GENOME")
echo "Detected contig name: $CONTIG"

# Replace chromosome name in BED with actual contig name

awk -v OFS='\t' -v CONTIG="$CONTIG" '{
    $1 = CONTIG; print
}' "$BED" > final_ranges_corrected.bed

# Count ranges

n_ranges=$(wc -l < final_ranges_corrected.bed)
echo "Total ranges to extract: $n_ranges"

# Count by strand

echo ""
echo "By strand:"
cut -f6 final_ranges_corrected.bed | sort | uniq -c

# Extract sequences using bedtools

bedtools getfasta \
    -fi "$GENOME" \
    -bed final_ranges_corrected.bed \
    -name \
    -s \
    -fo "$OUT"

# Summary

echo ""
echo "Sequences extracted: $(grep -c '^>' "$OUT")"
echo "Output: $OUT"

