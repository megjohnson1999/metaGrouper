# MetaGrouper Recommended Metadata Variables

## 🎯 Quick Start: Manual Variable Selection

For immediate analysis, use these biologically relevant variables:

```bash
python metagrouper.py /path/to/fastq/files \
    --metadata metadata_with_samples_final.csv \
    --variables celiacs_group case_control Sex "Delivery Mode" HLA "Dx Status" month Country GEMM PID "Onset Dx" disease_onset_details \
    --output focused_analysis
```

## 📊 Variable Categories

### ✅ **Highly Recommended (Primary Biological Variables)**
- `celiacs_group` - Disease classification
- `case_control` - Case/control status  
- `Sex` - Demographic factor
- `Delivery Mode` - Birth method (major microbiome factor!)
- `HLA` - Genetic markers
- `Dx Status` - Diagnosis status

### 🔬 **Recommended (Secondary Variables)**
- `month` - Temporal progression/age
- `Country` - Geographic origin
- `GEMM` - Patient/study identifier (for grouping)
- `PID` - Patient ID (for grouping)
- `Onset Dx` - Disease onset timing
- `disease_onset_details` - Disease specifics
- `Classification` - Sample classification

### ❌ **Exclude (Technical/Administrative)**
- `Unnamed: 0` - Row numbers
- `Plate`, `Well` - Lab processing coordinates  
- `Name`, `Golay.Barcode` - Technical identifiers
- `Collaborator.ID` - Administrative ID
- `Sample.Name.External.ID` - Technical ID
- `Sample_ID` - Used for matching only
- `PIDu`, `PIDcc` - Derived/duplicate variables
- `Study ID` - Administrative (unless needed for multi-study analysis)

## 🧬 **Variable Selection Strategy**

### **For Celiac Disease Research:**
Focus on disease-related variables:
```bash
--variables celiacs_group case_control "Dx Status" "Onset Dx" HLA Sex month
```

### **For Microbiome Development:**
Focus on developmental factors:
```bash
--variables "Delivery Mode" Sex month celiacs_group Country
```

### **For Patient Grouping:**
Include patient identifiers for co-assembly:
```bash
--variables GEMM PID celiacs_group case_control Sex
```

## 📈 **Expected Results**

With proper variable selection, you should see:
- **Significant PERMANOVA results** for disease-related variables
- **Meaningful clustering** based on biological factors
- **Assembly recommendations** grouped by patient or disease status
- **Clean visualizations** without technical noise

## 🚀 **Auto-Filtering (Coming Soon)**

MetaGrouper will soon support automatic variable filtering:
```bash
--auto-filter-variables  # Automatically detect biological variables
--exclude-variables Plate Well "Golay.Barcode"  # Manual exclusions
```