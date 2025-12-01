#!/usr/bin/env bash

#SBATCH --job-name=blast_all_vs_all
#SBATCH --partition=c
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=blast_all_vs_all_%j.out
#SBATCH --error=blast_all_vs_all_%j.err
#SBATCH --mail-user=wagner.guzman@gmi.oeaw.ac.at
#SBATCH --mail-type=FAIL

# Perform BLAST all-vs-all on final sequences
# Args: Input FASTA, Output prefix, E-value threshold

FASTA="$1"
OUT_PREFIX="$2"
EVALUE="${3:-1e-5}"

echo "Creating BLAST database..."
makeblastdb \
    -in "$FASTA" \
    -dbtype nucl \
    -out "${OUT_PREFIX}_db" \
    -title "Final_Ranges_DB"

echo ""
echo "Running BLASTn all-vs-all..."
blastn \
    -query "$FASTA" \
    -db "${OUT_PREFIX}_db" \
    -out "${OUT_PREFIX}.blastn" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue "$EVALUE" \
    -num_threads 4 \
    -task blastn

# Filter out self-hits

echo ""
echo "Filtering self-hits..."
awk '$1 != $2' "${OUT_PREFIX}.blastn" > "${OUT_PREFIX}.blastn.filtered"

# Summary statistics

total_hits=$(wc -l < "${OUT_PREFIX}.blastn")
filtered_hits=$(wc -l < "${OUT_PREFIX}.blastn.filtered")
self_hits=$((total_hits - filtered_hits))

echo ""
echo "=== BLAST SUMMARY ==="
echo "Total hits: $total_hits"
echo "Self-hits removed: $self_hits"
echo "Filtered hits: $filtered_hits"
echo ""
echo "Output files:"
echo "  Full BLAST: ${OUT_PREFIX}.blastn"
echo "  Filtered (no self-hits): ${OUT_PREFIX}.blastn.filtered"

# Top hits preview

echo ""
echo "Top 10 hits (by bitscore):"
sort -k12,12nr "${OUT_PREFIX}.blastn.filtered" | head -10 | \
    awk 'BEGIN{print "Query\tSubject\tPident\tLength\tBitscore"} 
         {print $1"\t"$2"\t"$3"\t"$4"\t"$12}'

