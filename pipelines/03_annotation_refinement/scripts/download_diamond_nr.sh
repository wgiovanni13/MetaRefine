#!/usr/bin/env bash

#SBATCH --job-name=download_diamond_nr
#SBATCH --partition=c
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=download_diamond_nr_%j.out
#SBATCH --error=download_diamond_nr_%j.err
#SBATCH --mail-user=wagner.guzman@gmi.oeaw.ac.at
#SBATCH --mail-type=FAIL

# Loading modules

ml build-env/f2022
ml diamond/2.1.8-gcc-12.2.0

# Create directories

OUTPUT_DIR="/groups/berger/lab/Giovanni.Guzman/DIAMOND_db"
mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

echo "Downloading pre-formatted DIAMOND NR database..."

# Download NR database in FASTA format

wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

# Convert into diamond

diamond makedb --in nr.gz -d nr --threads 2

echo "Done!"
