#!/bin/bash

# MetaGrouper Phase 1: Systematic Biological Validation
# Using well-documented, small samples from known studies

set -e
echo "🧬 MetaGrouper Phase 1: Biological Validation"
echo "=============================================="

cd /Users/meganjohnson/Documents/Megan/metaGrouper_development/phase1_testing

# Create directory structure
mkdir -p samples/{human_gut,marine,soil,misc}
mkdir -p metadata

# Function to download with size control
download_small_sample() {
    local accession=$1
    local environment=$2
    local description=$3
    local max_reads=${4:-5000}  # Default 5000 reads
    
    echo "📥 Downloading $accession ($description)"
    
    if fastq-dump --split-files --gzip --maxSpotId $max_reads \
                  --outdir samples/$environment $accession 2>/dev/null; then
        echo "✅ $accession downloaded successfully"
        ls -lh samples/$environment/$accession*
        return 0
    else
        echo "❌ $accession failed to download"
        return 1
    fi
}

echo ""
echo "Phase 1a: Testing known working samples"
echo "======================================="

# Start with samples we know work
echo "🦠 Testing Human Gut samples..."
download_small_sample "SRR1293827" "human_gut" "American Gut validated sample" 5000

# Let's try a few more from the same study
download_small_sample "SRR1293828" "human_gut" "American Gut sample 2" 5000 || true
download_small_sample "SRR1293829" "human_gut" "American Gut sample 3" 5000 || true

echo ""
echo "🌊 Testing Marine samples..."
# Try some marine samples (adjust these based on what works)
download_small_sample "SRR1565503" "marine" "Marine sample test 1" 5000 || true
download_small_sample "SRR1565504" "marine" "Marine sample test 2" 5000 || true

echo ""
echo "🌱 Testing Soil samples..."
# Try some soil samples
download_small_sample "SRR1946637" "soil" "Soil sample test 1" 5000 || true
download_small_sample "SRR1946638" "soil" "Soil sample test 2" 5000 || true

echo ""
echo "📊 Phase 1a Results"
echo "==================="

total_samples=0
for env in human_gut marine soil misc; do
    count=$(find samples/$env -name "*.fastq.gz" 2>/dev/null | wc -l)
    samples=$((count / 2))  # Divide by 2 for paired-end
    total_samples=$((total_samples + samples))
    size=$(du -sh samples/$env 2>/dev/null | cut -f1 || echo "0B")
    echo "$env: $samples samples ($count files), $size"
done

echo ""
echo "🎯 Total samples: $total_samples"

if [ $total_samples -ge 5 ]; then
    echo ""
    echo "✅ SUCCESS: Phase 1a complete with $total_samples samples!"
    echo ""
    echo "🚀 Ready for MetaGrouper testing!"
    echo ""
    echo "Next step: Run MetaGrouper on these samples:"
    echo "   python ../metagrouper.py samples/ -o phase1_results/ --kmer-size 21 --scaled 1000"
    
else
    echo ""
    echo "⚠️  Phase 1a incomplete. Only $total_samples samples downloaded."
    echo ""
    echo "Next steps:"
    echo "1. Check internet connection"
    echo "2. Try alternative sample IDs"  
    echo "3. Proceed with available samples for initial testing"
fi

echo ""
echo "Files downloaded:"
find samples -name "*.fastq.gz" | head -10