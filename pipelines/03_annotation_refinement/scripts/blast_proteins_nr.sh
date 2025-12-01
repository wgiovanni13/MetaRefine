#!/usr/bin/env bash

#SBATCH --job-name=blast_proteins_nr
#SBATCH --partition=m
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=translate6_%j.out
#SBATCH --error=translate6_%j.err
#SBATCH --mail-user=wagner.guzman@gmi.oeaw.ac.at
#SBATCH --mail-type=FAIL

# BLASTp against NR database with top 10 hits per query

QUERY="$1"        # final_proteins.faa
OUTFILE="$2"      # output file
NR_DB="$3"        # path to NR database
THREADS="${4:-8}"
EVALUE="${5:-1e-5}"

echo "=== BLASTp against NR database ==="
echo "Query: $QUERY"
echo "Output: $OUTFILE"
echo "Database: $NR_DB"
echo "Threads: $THREADS"
echo "E-value: $EVALUE"
echo ""

# Count queries

NUM_QUERIES=$(grep -c '^>' "$QUERY")
echo "Number of queries: $NUM_QUERIES"
echo ""

# Run BLASTp

echo "Running BLASTp..."
blastp \
    -query "$QUERY" \
    -db "$NR_DB" \
    -out "$OUTFILE" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle" \
    -evalue "$EVALUE" \
    -num_threads "$THREADS" \
    -max_target_seqs 10 \
    -max_hsps 1

echo ""
echo "BLASTp completed!"
echo "Total hits: $(wc -l < "$OUTFILE")"
echo "Queries with hits: $(cut -f1 "$OUTFILE" | sort -u | wc -l)"
echo ""

# Summary statistics

echo "=== Summary per query ==="
awk '{print $1}' "$OUTFILE" | sort | uniq -c | sort -rn | head -20

