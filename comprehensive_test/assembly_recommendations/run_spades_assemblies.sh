#!/bin/bash

# MetaGrouper SPADES Assembly Commands
# Generated assembly strategy: grouped

# Group: patient_id_P001
python -c "import shutil; out=open('patient_id_P001_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_001.fastq" "sample_002.fastq" "sample_003.fastq"]]; out.close()"
spades.py --meta -s patient_id_P001_combined.fastq -o patient_id_P001_spades_assembly

# Group: patient_id_P002
python -c "import shutil; out=open('patient_id_P002_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_004.fastq" "sample_005.fastq" "sample_006.fastq"]]; out.close()"
spades.py --meta -s patient_id_P002_combined.fastq -o patient_id_P002_spades_assembly

# Group: patient_id_P003
python -c "import shutil; out=open('patient_id_P003_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_008.fastq" "sample_007.fastq" "sample_009.fastq"]]; out.close()"
spades.py --meta -s patient_id_P003_combined.fastq -o patient_id_P003_spades_assembly

# Group: patient_id_P004
python -c "import shutil; out=open('patient_id_P004_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_010.fastq" "sample_011.fastq" "sample_012.fastq"]]; out.close()"
spades.py --meta -s patient_id_P004_combined.fastq -o patient_id_P004_spades_assembly

# Group: patient_id_P005
python -c "import shutil; out=open('patient_id_P005_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_014.fastq" "sample_013.fastq" "sample_015.fastq"]]; out.close()"
spades.py --meta -s patient_id_P005_combined.fastq -o patient_id_P005_spades_assembly

# Group: patient_id_P006
python -c "import shutil; out=open('patient_id_P006_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_016.fastq" "sample_018.fastq" "sample_017.fastq"]]; out.close()"
spades.py --meta -s patient_id_P006_combined.fastq -o patient_id_P006_spades_assembly

# Group: patient_id_P007
python -c "import shutil; out=open('patient_id_P007_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_021.fastq" "sample_019.fastq" "sample_020.fastq"]]; out.close()"
spades.py --meta -s patient_id_P007_combined.fastq -o patient_id_P007_spades_assembly

# Group: patient_id_P008
python -c "import shutil; out=open('patient_id_P008_combined.fastq', 'wb'); [shutil.copyfileobj(open(f, 'rb'), out) for f in ["sample_022.fastq" "sample_023.fastq" "sample_024.fastq"]]; out.close()"
spades.py --meta -s patient_id_P008_combined.fastq -o patient_id_P008_spades_assembly

