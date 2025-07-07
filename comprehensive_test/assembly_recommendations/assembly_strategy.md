# MetaGrouper Assembly Strategy Recommendations

## Overall Strategy: Grouped

**Confidence Score:** 0.41

**Rationale:** Mixed strategy: 8 co-assembly groups covering 24/24 samples.

**Primary Grouping Criterion:** patient_id

## Strategy Summary

- **Total Assemblies:** 8
- **Co-assemblies:** 8
- **Individual Assemblies:** 0
- **Samples in Co-assembly:** 24 (100.0%)

## Recommended Assembly Groups

### Group 1: patient_id_P001

- **Samples:** sample_001, sample_002, sample_003
- **Grouping Criterion:** patient_id
- **Criterion Value:** P001
- **Average Distance:** 0.938
- **Confidence Score:** 0.41

**Expected Benefits:**
- Biologically meaningful grouping by patient_id
- Reduced inter-sample contamination
- Better representation of group-specific features

**Expected Challenges:**
- May miss cross-group shared sequences
- Groups based on patient_id may have variable quality
- High within-group diversity (avg dist: 0.938)

### Group 2: patient_id_P002

- **Samples:** sample_004, sample_005, sample_006
- **Grouping Criterion:** patient_id
- **Criterion Value:** P002
- **Average Distance:** 0.924
- **Confidence Score:** 0.41

**Expected Benefits:**
- Biologically meaningful grouping by patient_id
- Reduced inter-sample contamination
- Better representation of group-specific features

**Expected Challenges:**
- May miss cross-group shared sequences
- Groups based on patient_id may have variable quality
- High within-group diversity (avg dist: 0.924)

### Group 3: patient_id_P003

- **Samples:** sample_008, sample_007, sample_009
- **Grouping Criterion:** patient_id
- **Criterion Value:** P003
- **Average Distance:** 0.939
- **Confidence Score:** 0.41

**Expected Benefits:**
- Biologically meaningful grouping by patient_id
- Reduced inter-sample contamination
- Better representation of group-specific features

**Expected Challenges:**
- May miss cross-group shared sequences
- Groups based on patient_id may have variable quality
- High within-group diversity (avg dist: 0.939)

### Group 4: patient_id_P004

- **Samples:** sample_010, sample_011, sample_012
- **Grouping Criterion:** patient_id
- **Criterion Value:** P004
- **Average Distance:** 0.939
- **Confidence Score:** 0.41

**Expected Benefits:**
- Biologically meaningful grouping by patient_id
- Reduced inter-sample contamination
- Better representation of group-specific features

**Expected Challenges:**
- May miss cross-group shared sequences
- Groups based on patient_id may have variable quality
- High within-group diversity (avg dist: 0.939)

### Group 5: patient_id_P005

- **Samples:** sample_014, sample_013, sample_015
- **Grouping Criterion:** patient_id
- **Criterion Value:** P005
- **Average Distance:** 0.937
- **Confidence Score:** 0.41

**Expected Benefits:**
- Biologically meaningful grouping by patient_id
- Reduced inter-sample contamination
- Better representation of group-specific features

**Expected Challenges:**
- May miss cross-group shared sequences
- Groups based on patient_id may have variable quality
- High within-group diversity (avg dist: 0.937)

### Group 6: patient_id_P006

- **Samples:** sample_016, sample_018, sample_017
- **Grouping Criterion:** patient_id
- **Criterion Value:** P006
- **Average Distance:** 0.938
- **Confidence Score:** 0.41

**Expected Benefits:**
- Biologically meaningful grouping by patient_id
- Reduced inter-sample contamination
- Better representation of group-specific features

**Expected Challenges:**
- May miss cross-group shared sequences
- Groups based on patient_id may have variable quality
- High within-group diversity (avg dist: 0.938)

### Group 7: patient_id_P007

- **Samples:** sample_021, sample_019, sample_020
- **Grouping Criterion:** patient_id
- **Criterion Value:** P007
- **Average Distance:** 0.854
- **Confidence Score:** 0.41

**Expected Benefits:**
- Biologically meaningful grouping by patient_id
- Reduced inter-sample contamination
- Better representation of group-specific features

**Expected Challenges:**
- May miss cross-group shared sequences
- Groups based on patient_id may have variable quality
- High within-group diversity (avg dist: 0.854)

### Group 8: patient_id_P008

- **Samples:** sample_022, sample_023, sample_024
- **Grouping Criterion:** patient_id
- **Criterion Value:** P008
- **Average Distance:** 0.869
- **Confidence Score:** 0.41

**Expected Benefits:**
- Biologically meaningful grouping by patient_id
- Reduced inter-sample contamination
- Better representation of group-specific features

**Expected Challenges:**
- May miss cross-group shared sequences
- Groups based on patient_id may have variable quality
- High within-group diversity (avg dist: 0.869)

## Assembly Commands

### MEGAHIT

**patient_id_P001:**
```bash
megahit -r sample_001.fastq,sample_002.fastq,sample_003.fastq -o patient_id_P001_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99
```

**patient_id_P002:**
```bash
megahit -r sample_004.fastq,sample_005.fastq,sample_006.fastq -o patient_id_P002_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99
```

**patient_id_P003:**
```bash
megahit -r sample_008.fastq,sample_007.fastq,sample_009.fastq -o patient_id_P003_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99
```

**patient_id_P004:**
```bash
megahit -r sample_010.fastq,sample_011.fastq,sample_012.fastq -o patient_id_P004_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99
```

**patient_id_P005:**
```bash
megahit -r sample_014.fastq,sample_013.fastq,sample_015.fastq -o patient_id_P005_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99
```

**patient_id_P006:**
```bash
megahit -r sample_016.fastq,sample_018.fastq,sample_017.fastq -o patient_id_P006_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99
```

**patient_id_P007:**
```bash
megahit -r sample_021.fastq,sample_019.fastq,sample_020.fastq -o patient_id_P007_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99
```

**patient_id_P008:**
```bash
megahit -r sample_022.fastq,sample_023.fastq,sample_024.fastq -o patient_id_P008_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99
```

### SPADES

**patient_id_P001:**
```bash
python -c "import shutil; out=open('patient_id_P001_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_001.fastq" "sample_002.fastq" "sample_003.fastq"]]; out.close()"
```

```bash
spades.py --meta -s patient_id_P001_combined.fastq -o patient_id_P001_spades_assembly
```

**patient_id_P002:**
```bash
python -c "import shutil; out=open('patient_id_P002_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_004.fastq" "sample_005.fastq" "sample_006.fastq"]]; out.close()"
```

```bash
spades.py --meta -s patient_id_P002_combined.fastq -o patient_id_P002_spades_assembly
```

**patient_id_P003:**
```bash
python -c "import shutil; out=open('patient_id_P003_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_008.fastq" "sample_007.fastq" "sample_009.fastq"]]; out.close()"
```

```bash
spades.py --meta -s patient_id_P003_combined.fastq -o patient_id_P003_spades_assembly
```

**patient_id_P004:**
```bash
python -c "import shutil; out=open('patient_id_P004_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_010.fastq" "sample_011.fastq" "sample_012.fastq"]]; out.close()"
```

```bash
spades.py --meta -s patient_id_P004_combined.fastq -o patient_id_P004_spades_assembly
```

**patient_id_P005:**
```bash
python -c "import shutil; out=open('patient_id_P005_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_014.fastq" "sample_013.fastq" "sample_015.fastq"]]; out.close()"
```

```bash
spades.py --meta -s patient_id_P005_combined.fastq -o patient_id_P005_spades_assembly
```

**patient_id_P006:**
```bash
python -c "import shutil; out=open('patient_id_P006_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_016.fastq" "sample_018.fastq" "sample_017.fastq"]]; out.close()"
```

```bash
spades.py --meta -s patient_id_P006_combined.fastq -o patient_id_P006_spades_assembly
```

**patient_id_P007:**
```bash
python -c "import shutil; out=open('patient_id_P007_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_021.fastq" "sample_019.fastq" "sample_020.fastq"]]; out.close()"
```

```bash
spades.py --meta -s patient_id_P007_combined.fastq -o patient_id_P007_spades_assembly
```

**patient_id_P008:**
```bash
python -c "import shutil; out=open('patient_id_P008_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_022.fastq" "sample_023.fastq" "sample_024.fastq"]]; out.close()"
```

```bash
spades.py --meta -s patient_id_P008_combined.fastq -o patient_id_P008_spades_assembly
```

