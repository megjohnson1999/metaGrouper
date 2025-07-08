#!/bin/bash

#SBATCH --job-name=sourmash
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --mail-user megan.j@wustl.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


source /ref/sahlab/software/miniforge3/bin/activate
conda activate sourmash

# Define directories and the metadata CSV file
INPUT_DIR="test_0415.out/results/output/host_removed"
OUTPUT_DIR="cluster_metagenomes"
SIG_DIR="${OUTPUT_DIR}/signatures"
TEMP_DIR="${OUTPUT_DIR}/combined_reads"
METADATA_FILE="metadata.csv"

# Create output directories
mkdir -p ${SIG_DIR} ${TEMP_DIR}

# Read the metadata CSV file and store sample metadata in an associative array
declare -A sample_metadata

# Read metadata CSV (assuming it has 'sample_name' and 'metadata' columns)
tail -n +2 ${METADATA_FILE} | while IFS=',' read -r sample_name metadata; do
    # Store metadata for each sample (assuming metadata is in the second column)
    sample_metadata["${sample_name}"]="${metadata}"
done

# Process all paired-end samples in the input directory
# Find all R1 files and extract sample names
find ${INPUT_DIR} -name "*_R1.fastq*" | while read r1_file; do
    # Extract base filename without path and R1 suffix
    base_name=$(basename ${r1_file})
    sample_name=${base_name/_R1.fastq*/}

    # Try to match the "GEMM_*" part of the name from sourmash to the metadata
    matched_metadata=""

    # Extract the GEMM part from the sourmash sample name (e.g., "GEMM_010_18M")
    gemm_part=$(echo ${sample_name} | grep -o 'GEMM_[0-9]*_[0-9]*[A-Za-z]*')

    # Remove the "_hr" suffix if it exists
    gemm_part=$(echo ${gemm_part} | sed 's/_hr$//')

    # Look for a match in the metadata
    for key in "${!sample_metadata[@]}"; do
        if [[ "${key}" == *"${gemm_part}"* ]]; then
            matched_metadata="${sample_metadata[${key}]}"
            break
        fi
    done

    if [ -z "$matched_metadata" ]; then
        echo "No matching metadata for sample: ${sample_name}. Skipping."
        continue
    fi

    # Construct R2 filename
    r2_file=${r1_file/_R1.fastq/_R2.fastq}
    echo "Processing sample: ${sample_name} with metadata: ${matched_metadata}"

    # Combine paired reads
    cat ${r1_file} ${r2_file} > ${TEMP_DIR}/${sample_name}_combined.fastq

    # Compute signature with explicit name parameter and metadata
    sourmash compute -k 31 --scaled 1000 --track-abundance \
        --name "${sample_name}_${matched_metadata}" \
        ${TEMP_DIR}/${sample_name}_combined.fastq \
        -o ${SIG_DIR}/${sample_name}.sig

    # Remove combined file to save space
    rm ${TEMP_DIR}/${sample_name}_combined.fastq
done

# Compare all signatures
echo "Comparing all samples..."
sourmash compare ${SIG_DIR}/*.sig -o ${OUTPUT_DIR}/distance_matrix.csv

# Create a directory for the visualizations
mkdir -p ${OUTPUT_DIR}/visualization

# Generate visualizations using sourmash plot
sourmash plot ${OUTPUT_DIR}/distance_matrix.csv --output-dir ${OUTPUT_DIR}/visualization

conda deactivate
conda deactivate

