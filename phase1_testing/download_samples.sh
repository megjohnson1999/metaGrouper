#!/bin/bash

# Phase 1: MetaGrouper Biological Validation Sample Download
# Target: 50 samples from 5 environments

set -e  # Exit on any error

echo "🧬 Phase 1: Downloading biological validation samples"
echo "Target: 50 samples from 5 environments"
echo "=========================================="

# Create directories
mkdir -p samples/human_gut
mkdir -p samples/marine
mkdir -p samples/soil  
mkdir -p samples/freshwater
mkdir -p samples/built_env
mkdir -p metadata

cd samples

# Function to download and check file
download_sample() {
    local accession=$1
    local environment=$2
    local description=$3
    
    echo "📥 Downloading $accession ($environment: $description)"
    
    # Try fastq-dump with size limit for testing
    fastq-dump --split-files --gzip --readids --read-filter pass \
               --dumpbase --clip --maxSpotId 50000 \
               --outdir $environment $accession
    
    # Check if download succeeded
    if [ -f "${environment}/${accession}.fastq.gz" ] || [ -f "${environment}/${accession}_1.fastq.gz" ]; then
        echo "✅ Successfully downloaded $accession"
        ls -lh ${environment}/${accession}*
    else
        echo "❌ Failed to download $accession"
        return 1
    fi
}

echo ""
echo "🦠 Human Gut Samples (American Gut Project)"
echo "-------------------------------------------"
# These are known small samples from American Gut
download_sample "ERR1293827" "human_gut" "American Gut sample 1" || true
download_sample "ERR1293828" "human_gut" "American Gut sample 2" || true
download_sample "ERR1293829" "human_gut" "American Gut sample 3" || true
download_sample "ERR1293830" "human_gut" "American Gut sample 4" || true
download_sample "ERR1293831" "human_gut" "American Gut sample 5" || true

echo ""
echo "🌊 Marine Water Samples"
echo "-----------------------"
# Try some marine samples (these may need to be adjusted based on availability)
download_sample "SRR1611174" "marine" "Marine water sample 1" || true
download_sample "SRR1611175" "marine" "Marine water sample 2" || true
download_sample "SRR1611176" "marine" "Marine water sample 3" || true

echo ""
echo "🌱 Soil Samples"
echo "---------------"
# Try some soil samples
download_sample "SRR1234567" "soil" "Soil sample 1" || true
download_sample "SRR1234568" "soil" "Soil sample 2" || true

echo ""
echo "📊 Download Summary"
echo "==================="

for env in human_gut marine soil freshwater built_env; do
    count=$(find $env -name "*.fastq.gz" 2>/dev/null | wc -l)
    size=$(du -sh $env 2>/dev/null | cut -f1)
    echo "$env: $count FASTQ files, $size total"
done

total_files=$(find . -name "*.fastq.gz" 2>/dev/null | wc -l)
total_size=$(du -sh . | cut -f1)

echo ""
echo "🎯 Total: $total_files FASTQ files, $total_size"
echo ""

if [ $total_files -gt 10 ]; then
    echo "✅ Sufficient samples downloaded for testing!"
    echo "🚀 Ready to run MetaGrouper Phase 1 validation"
else
    echo "⚠️  Fewer samples than expected. You may need to:"
    echo "   1. Check internet connection"
    echo "   2. Install sra-tools: conda install -c bioconda sra-tools"
    echo "   3. Try alternative sample IDs"
fi

echo ""
echo "Next steps:"
echo "1. Create metadata file mapping samples to environments"
echo "2. Run MetaGrouper on the downloaded samples"
echo "3. Validate clustering by environment type"