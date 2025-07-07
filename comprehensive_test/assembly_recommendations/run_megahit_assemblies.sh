#!/bin/bash

# MetaGrouper MEGAHIT Assembly Commands
# Generated assembly strategy: grouped

# Group: patient_id_P001
megahit -r sample_001.fastq,sample_002.fastq,sample_003.fastq -o patient_id_P001_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99

# Group: patient_id_P002
megahit -r sample_004.fastq,sample_005.fastq,sample_006.fastq -o patient_id_P002_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99

# Group: patient_id_P003
megahit -r sample_008.fastq,sample_007.fastq,sample_009.fastq -o patient_id_P003_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99

# Group: patient_id_P004
megahit -r sample_010.fastq,sample_011.fastq,sample_012.fastq -o patient_id_P004_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99

# Group: patient_id_P005
megahit -r sample_014.fastq,sample_013.fastq,sample_015.fastq -o patient_id_P005_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99

# Group: patient_id_P006
megahit -r sample_016.fastq,sample_018.fastq,sample_017.fastq -o patient_id_P006_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99

# Group: patient_id_P007
megahit -r sample_021.fastq,sample_019.fastq,sample_020.fastq -o patient_id_P007_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99

# Group: patient_id_P008
megahit -r sample_022.fastq,sample_023.fastq,sample_024.fastq -o patient_id_P008_coassembly --min-contig-len 500 --k-list 21,29,39,59,79,99

