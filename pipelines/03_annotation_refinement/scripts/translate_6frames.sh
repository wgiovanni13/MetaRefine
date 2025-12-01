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

# Routes for Inputs and outputs

IN="$1"                 # cat3_ranges_extracted.fasta
OUT="$2"                # cat3_ranges_extracted.6frames.archaea.faa
TABLE="${3:-11}"        # Genetic Code: 11 (Bacteria/Archaea)

# Needed modules

ml build-env/2020
ml emboss/6.6.0-foss-2018b

transeq -sequence "$IN" -outseq "$OUT" -frame 6 -table "$TABLE" -auto

in_n=$(grep -c '^>' "$IN" || true)
out_n=$(grep -c '^>' "$OUT" || true)
echo "Input sequences        : $in_n"
echo "Translated sequences   : $out_n (expected ~ 6x)"
echo "Output                 : $OUT"
