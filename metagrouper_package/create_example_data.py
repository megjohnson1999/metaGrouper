#!/usr/bin/env python3
"""
Create example dataset for MetaGrouper testing.

Generates realistic viral metagenomic samples with associated metadata
for testing all three phases of MetaGrouper.
"""

import random
import tempfile
from pathlib import Path
import pandas as pd
import numpy as np


def generate_viral_sequences(n_sequences: int = 100, seq_length: int = 150) -> list:
    """Generate realistic viral-like sequences."""
    
    # Common viral k-mer patterns (simplified)
    viral_patterns = [
        "ATCGATCGATCG",  # Simple repeat
        "GCTAGCTAGCTA",  # Another repeat
        "AAATTTAAATTT",  # AT-rich region
        "CCCGGGCCCGGG",  # GC-rich region
        "TGATCAGATCAG",  # Mixed pattern
    ]
    
    sequences = []
    
    for i in range(n_sequences):
        # Choose a base pattern
        base_pattern = random.choice(viral_patterns)
        
        # Extend to desired length
        sequence = ""
        while len(sequence) < seq_length:
            if random.random() < 0.7:  # 70% from pattern
                sequence += base_pattern
            else:  # 30% random
                sequence += random.choice("ATCG")
        
        # Trim to exact length
        sequence = sequence[:seq_length]
        
        # Add some mutations
        seq_list = list(sequence)
        mutation_rate = random.uniform(0.01, 0.05)  # 1-5% mutations
        for j in range(len(seq_list)):
            if random.random() < mutation_rate:
                seq_list[j] = random.choice("ATCG")
        
        sequences.append(''.join(seq_list))
    
    return sequences


def create_sample_fastq(sample_name: str, patient_group: int, timepoint: str, 
                       output_dir: Path, n_reads: int = 1000) -> str:
    """Create a FASTQ file for one sample."""
    
    # Generate sequences with some similarity based on patient group
    sequences = []
    
    # Group-specific patterns to create similarity
    if patient_group == 1:
        # Group 1: More AT-rich
        base_sequences = generate_viral_sequences(n_reads // 4, 150)
        for seq in base_sequences:
            # Create variations
            for _ in range(4):
                seq_list = list(seq)
                # Add AT bias
                for i in range(len(seq_list)):
                    if random.random() < 0.1:
                        seq_list[i] = random.choice("AT")
                sequences.append(''.join(seq_list))
    
    elif patient_group == 2:
        # Group 2: More GC-rich
        base_sequences = generate_viral_sequences(n_reads // 4, 150)
        for seq in base_sequences:
            for _ in range(4):
                seq_list = list(seq)
                # Add GC bias
                for i in range(len(seq_list)):
                    if random.random() < 0.1:
                        seq_list[i] = random.choice("GC")
                sequences.append(''.join(seq_list))
    
    else:
        # Group 3: Mixed/diverse
        sequences = generate_viral_sequences(n_reads, 150)
    
    # Add timepoint-specific variations
    if timepoint == "week4" or timepoint == "week8":
        # Later timepoints have more diversity (viral evolution)
        for i in range(len(sequences)):
            seq_list = list(sequences[i])
            evolution_rate = 0.02 if timepoint == "week4" else 0.04
            for j in range(len(seq_list)):
                if random.random() < evolution_rate:
                    seq_list[j] = random.choice("ATCG")
            sequences[i] = ''.join(seq_list)
    
    # Write FASTQ file
    fastq_file = output_dir / f"{sample_name}.fastq"
    with open(fastq_file, 'w') as f:
        for i, seq in enumerate(sequences):
            # FASTQ format
            f.write(f"@{sample_name}_read_{i+1}\n")
            f.write(f"{seq}\n")
            f.write(f"+\n")
            f.write(f"{'I' * len(seq)}\n")  # High quality scores
    
    return str(fastq_file)


def create_metadata(samples_info: list, output_dir: Path) -> str:
    """Create realistic metadata file."""
    
    metadata_records = []
    
    for sample_name, patient_id, timepoint, treatment, location in samples_info:
        record = {
            'sample_id': sample_name,
            'patient_id': patient_id,
            'timepoint': timepoint,
            'treatment': treatment,
            'location': location,
            'age': random.randint(25, 75),
            'gender': random.choice(['Male', 'Female']),
            'bmi': round(random.uniform(18.5, 35.0), 1),
            'viral_load': random.randint(1000, 100000),
            'collection_date': f"2023-{random.randint(1, 12):02d}-{random.randint(1, 28):02d}"
        }
        metadata_records.append(record)
    
    metadata_df = pd.DataFrame(metadata_records)
    metadata_file = output_dir / "realistic_metadata.csv"
    metadata_df.to_csv(metadata_file, index=False)
    
    return str(metadata_file)


def main():
    """Create complete example dataset."""
    
    print("🧬 Creating MetaGrouper Example Dataset")
    print("=" * 50)
    
    output_dir = Path("example_data")
    output_dir.mkdir(exist_ok=True)
    
    # Define study design
    patients = ['P001', 'P002', 'P003', 'P004', 'P005', 'P006', 'P007', 'P008']
    timepoints = ['baseline', 'week4', 'week8']
    treatments = ['control', 'treatment_A', 'treatment_B']
    locations = ['site_A', 'site_B', 'site_C']
    
    samples_info = []
    sample_count = 0
    
    # Create longitudinal study design
    for patient in patients:
        # Assign patient to groups
        patient_num = int(patient[1:])
        patient_group = ((patient_num - 1) // 3) + 1  # Groups 1, 2, 3
        
        # Assign treatment and location
        treatment = treatments[patient_num % len(treatments)]
        location = locations[patient_num % len(locations)]
        
        # Create samples for each timepoint
        for timepoint in timepoints:
            sample_count += 1
            sample_name = f"sample_{sample_count:03d}"
            
            # Create FASTQ file
            print(f"Creating {sample_name} ({patient}, {timepoint})...")
            fastq_path = create_sample_fastq(
                sample_name, patient_group, timepoint, output_dir
            )
            
            samples_info.append((sample_name, patient, timepoint, treatment, location))
    
    # Create metadata file
    print("Creating metadata file...")
    metadata_path = create_metadata(samples_info, output_dir)
    
    # Create README for the example data
    readme_content = f"""# MetaGrouper Example Dataset

This synthetic dataset contains {len(samples_info)} viral metagenomic samples for testing MetaGrouper.

## Study Design

- **Patients:** {len(patients)} individuals
- **Timepoints:** {', '.join(timepoints)}
- **Treatments:** {', '.join(treatments)}
- **Locations:** {', '.join(locations)}

## Data Structure

### Samples
- Synthetic viral-like sequences (~1000 reads per sample)
- Group-specific patterns to test similarity detection
- Temporal evolution to test longitudinal analysis

### Metadata
- Patient demographics (age, gender, BMI)
- Clinical variables (viral load, collection date)
- Study design variables (treatment, location, timepoint)

## Expected Results

### Phase 1: K-mer Analysis
- Samples should cluster by patient group
- Some temporal patterns within patients

### Phase 2: Metadata Analysis
- **Significant variables:** patient_id, timepoint, treatment
- **Less significant:** demographic variables
- **Clustering:** Should identify patient groups

### Phase 3: Assembly Strategy
- **Recommended strategy:** Grouped assembly
- **Groups:** Based on patient similarity
- **Confidence:** Medium to high

## Usage

```bash
# Basic analysis
python metagrouper.py example_data/ -o test_results/

# Full analysis with metadata
python metagrouper.py example_data/ -m example_data/realistic_metadata.csv -o full_results/

# Memory-efficient for testing
python metagrouper.py example_data/ -o results/ --use-sketching --sketch-size 500
```

## Files

- `sample_*.fastq` - FASTQ files for each sample
- `realistic_metadata.csv` - Sample metadata
- `README.md` - This file

Generated with MetaGrouper example data generator.
"""
    
    with open(output_dir / "README.md", 'w') as f:
        f.write(readme_content)
    
    print(f"\n✅ Example dataset created!")
    print(f"📁 Location: {output_dir.absolute()}")
    print(f"📊 Samples: {len(samples_info)}")
    print(f"🧬 Patients: {len(patients)}")
    print(f"📋 Metadata: {metadata_path}")
    print(f"\n🚀 Test with:")
    print(f"python metagrouper.py example_data/ -m example_data/realistic_metadata.csv -o test_results/")


if __name__ == "__main__":
    main()