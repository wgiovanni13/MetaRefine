#!/usr/bin/env bash

#SBATCH --job-name=extract_cat3
#SBATCH --partition=c
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --output=extract_cat3_%j.out
#SBATCH --error=extract_cat3_%j.err
#SBATCH --mail-user=wagner.guzman@gmi.oeaw.ac.at
#SBATCH --mail-type=FAIL

# Args: Genome (.fasta), Cat3 (.tsv of 3rd category), Output (.fasta) and initial tsv to extract strand

GENOME="$1"
CAT3="$2"
OUT="$3"
COMBINED="$4"

# Modules needed

ml build-env/f2022
ml bedtools/2.31.1-gcc-13.2.0

# Contig (first FASTA header for bed format)

CONTIG=$(awk '/^>/{sub(/^>/,""); print $1; exit}' "$GENOME")

# BED (0 based start position) with header "Start-End"

awk -v OFS='\t' -v C="$CONTIG" '
  NR>1 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ {
    s=$1; e=$2; if (s>e){t=s; s=e; e=t}
    print C, s-1, e, s "-" e
  }' "$CAT3" \
| sort -t $'\t' -k1,1 -k2,2n -k3,3n -u > cat3_ranges.bed

# Extract FASTA (bedtools appends "::contig:start-end" after -name)

RAW="${OUT}.raw"
bedtools getfasta -fi "$GENOME" -bed cat3_ranges.bed -fo "$RAW" -name

# Build Start-End header with strand map from combined file

awk -F'\t' 'NR>1{
  if($1!="NA" && $2!="NA" && $7~/^[+-]$/) print $1"\t"$2"\t"$7;
  if($3!="NA" && $4!="NA" && $8~/^[+-]$/) print $3"\t"$4"\t"$8;
}' "$COMBINED" | sort -u > ranges_strand.map

# Rename headers: from headers from bed file >Start-End::... to >Start-End_(strand)

awk -v MAP="ranges_strand.map" '
BEGIN{
  FS="\t";
  while((getline<MAP)>0){ key=$1"-"$2; st[key]=$3 }
}

# Headers

/^>/{
  h=$0; sub(/^>/,"",h);
  split(h,a,"::"); r=a[1]; s=st[r];
  if(s=="") print ">" r;
  else      print ">" r "_(" s ")";
  next
}

# Sequences

{ print }
' "$RAW" > "$OUT"

# Summary

echo "Contig: $CONTIG"
echo "BED lines: $(wc -l < cat3_ranges.bed)"
echo "FASTA seqs: $(grep -c '^>' "$OUT")"
echo "Output: $OUT"

