# MetaGrouper Interactive Visualizations Guide

## Overview
MetaGrouper generates multiple interactive HTML files to provide different perspectives on your metagenomic data. Each file serves a specific purpose and offers unique insights.

## 📊 Interactive Files Generated

### 1. `interactive_enhanced.html` 🌟 **[RECOMMENDED]**
**Purpose**: Comprehensive multi-method data exploration  
**Features**:
- ✅ **Method Switching**: Toggle between PCA, t-SNE, and UMAP
- ✅ **Dynamic Coloring**: Color points by any metadata variable
- ✅ **Fixed Dropdowns**: Working color switching (no gray points)
- ✅ **Color Legends**: Automatic colorbar with variable names
- ✅ **Rich Hover Info**: Detailed sample information on hover
- ✅ **Built-in Zoom/Pan**: Native Plotly interactivity

**Best for**: 
- Comprehensive data exploration
- Finding patterns across different projection methods
- Understanding sample relationships and metadata effects

**Usage**: This is the most advanced visualization - **use this one first!**

---

### 2. `interactive_pca.html`
**Purpose**: Simple PCA-only exploration  
**Features**:
- PCA projection with metadata coloring dropdown
- Hover tooltips with sample details
- Basic zoom/pan functionality

**Best for**:
- Quick PCA analysis
- Users who prefer simpler interface
- When you only need PCA perspective

**Note**: May have older color dropdown implementation

---

### 3. `interactive_heatmap.html`
**Purpose**: Sample-to-sample similarity analysis  
**Features**:
- Interactive distance/similarity matrix
- Hover details showing pairwise relationships
- Zoom/pan to explore sample clusters
- Color-coded similarity values

**Best for**:
- Understanding which samples are most/least similar
- Identifying outliers or tight clusters
- Quality control and sample relationship verification

---

### 4. `interactive_dashboard.html`
**Purpose**: Combined overview  
**Features**:
- Side-by-side PCA and heatmap views
- Synchronized data exploration
- Split-screen layout for comparison

**Best for**:
- Getting both perspectives simultaneously
- Presentations and reports
- Quick overview of results

---

## 🎨 Color Legend Features

### Automatic Colorbars
- **Numerical Variables**: Continuous color scale with colorbar legend
- **Categorical Variables**: Discrete colors with proper mapping
- **Dynamic Titles**: Legend titles update when switching variables

### How to Use Colors
1. **Select Variable**: Use the "Color by" dropdown
2. **View Legend**: Colorbar appears on the right side
3. **Interpret Values**: 
   - Continuous scale for numerical data (Age, Values)
   - Discrete colors for categories (Treatment, Patient ID)

---

## 🚀 Usage Recommendations

### For First-Time Users
1. **Start with**: `interactive_enhanced.html`
2. **Try switching**: Between PCA, t-SNE, and UMAP
3. **Explore coloring**: Switch between different metadata variables
4. **Check relationships**: Use the heatmap for detailed sample comparisons

### For Different Analysis Goals

**🔍 Pattern Discovery**
- Use `interactive_enhanced.html`
- Try different projection methods
- Color by various metadata variables

**📏 Sample Similarity**
- Use `interactive_heatmap.html`
- Hover over cells for exact distances
- Look for cluster patterns

**📊 Quick Overview**
- Use `interactive_dashboard.html`
- Get both perspectives at once

**🎯 Simple Analysis**
- Use `interactive_pca.html`
- Focus on PCA results only

---

## 🔧 Technical Notes

### Browser Compatibility
- Works in all modern browsers (Chrome, Firefox, Safari, Edge)
- Self-contained HTML files (no internet connection needed)
- Files are ~4-5MB each (includes all data and libraries)

### Performance
- Optimized for datasets up to 100-200 samples
- Larger datasets may have slower interactions
- Use `--max-reads` parameter for testing with large datasets

### Sharing
- HTML files are completely portable
- Can be shared via email, cloud storage, or web hosting
- No dependencies or additional software required

---

## 🐛 Troubleshooting

### Color Dropdown Not Working
- **Problem**: Points turn gray when switching colors
- **Solution**: Use `interactive_enhanced.html` (has the fix)
- **Status**: Fixed in latest version

### Large File Sizes
- **Cause**: Self-contained HTML with embedded data
- **Normal Size**: 4-5MB for typical datasets
- **Reduction**: Use `--max-reads` for smaller test files

### Slow Performance
- **Cause**: Large datasets (>100 samples)
- **Solutions**: 
  - Use smaller k-mer sizes (`--kmer-size 15`)
  - Limit reads (`--max-reads 1000`)
  - Use faster computer/browser

---

## 📈 Best Practices

### Data Exploration Workflow
1. **Start broad**: Use enhanced HTML with all methods
2. **Identify patterns**: Look for clusters and outliers
3. **Test hypotheses**: Color by relevant metadata
4. **Verify relationships**: Cross-check with heatmap
5. **Document findings**: Save interesting views/discoveries

### Metadata Coloring Tips
- **Categorical variables**: Look for clear separation between groups
- **Numerical variables**: Look for gradients and correlation patterns
- **Temporal data**: Check if time-related variables show progression
- **Treatment effects**: Compare before/after or control/treatment

---

## 🎯 Quick Start Commands

```bash
# Generate all interactive visualizations
python metagrouper.py data/ -m metadata.csv --interactive -o results/

# Generate only interactive (skip static plots)
python metagrouper.py data/ -m metadata.csv --interactive-only -o results/

# Quick test with enhanced visualizations
python metagrouper.py data/ --interactive --kmer-size 15 --max-reads 500
```

The enhanced HTML file provides the most comprehensive and reliable interactive experience for exploring your metagenomic data! 🚀