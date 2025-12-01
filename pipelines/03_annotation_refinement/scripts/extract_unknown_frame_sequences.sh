#!/usr/bin/env bash

#SBATCH --job-name=extract_unknown_frame_sequences
#SBATCH --partition=c
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=extract_unknown_frame_sequences_%j.out
#SBATCH --error=extract_unknown_frame_sequences_%j.err
#SBATCH --mail-user=wagner.guzman@gmi.oeaw.ac.at
#SBATCH --mail-type=FAIL

# Extract sequences for ranges with unknown frame
# Args: Genome (.fasta), Merged BED file, Output (.fasta)

GENOME="$1"
MERGED_BED="$2"
OUT="$3"

# Extract the actual contig name from FASTA
CONTIG=$(awk '/^>/{sub(/^>/,""); print $1; exit}' "$GENOME")
echo "Detected contig name: $CONTIG"

# Filter ranges with unknown frame (contain '_unk')
grep '_unk' "$MERGED_BED" > unknown_frame_ranges_original.bed

# Replace chromosome name with actual contig name
awk -v OFS='\t' -v CONTIG="$CONTIG" '{
    $1 = CONTIG; print
}' unknown_frame_ranges_original.bed > unknown_frame_ranges.bed

# Count ranges
n_ranges=$(wc -l < unknown_frame_ranges.bed)
echo "Total ranges with unknown frame: $n_ranges" > unknown_frame_stats.txt

# Count by strand
echo "" >> unknown_frame_stats.txt
echo "By strand:" >> unknown_frame_stats.txt
cut -f6 unknown_frame_ranges.bed | sort | uniq -c >> unknown_frame_stats.txt

# Extract sequences using bedtools
bedtools getfasta \
    -fi "$GENOME" \
    -bed unknown_frame_ranges.bed \
    -name \
    -s \
    -fo "$OUT"

# Summary
echo "" >> unknown_frame_stats.txt
echo "Sequences extracted:" >> unknown_frame_stats.txt
grep -c '^>' "$OUT" >> unknown_frame_stats.txt

echo "Output: $OUT"
cat unknown_frame_stats.txt
