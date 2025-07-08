#!/bin/bash

#SBATCH --job-name=sourmash
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mail-user megan.j@wustl.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


source /ref/sahlab/software/miniforge3/bin/activate
conda activate sourmash


# Define directories
dir=cluster_metagenomes

# Create a directory for the visualizations
mkdir -p "$dir"/visualization

# Generate visualizations using sourmash plot
sourmash plot "$dir"/distance_matrix.csv --output-dir "$dir"/visualization


conda deactivate
conda deactivate

