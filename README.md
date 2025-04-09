# DNA Barcode Generator with Primer3 and BLAST

### Author: Anthony Glaude from Labo-MAB, Pre. Brunet

## Overview

This script generates DNA barcodes by creating primers with **Primer3** and assembling them into barcodes with randomly generated inserts. It performs a series of validation checks to ensure that the barcodes meet stringent criteria, including checks for Hamming distance, internal repeats, palindromes, and BLAST verification. The barcodes are BLASTed against a viral genome database and each other to avoid undesired similarities.

## Features

### Primer Creation with Primer3
- **Size**: 18–22 nucleotides
- **Melting Temperature (Tm)**: 57–63 °C 
- **GC Content**: 45–55%
- **BLAST Verification**: Primers are verified against a viral genome database to ensure no undesired similarity.
  - **E-value Threshold**: 0.1 (even slight similarity is rejected for primer specificity)

### Barcode Assembly
- A **200nt barcode** is assembled by concatenating the forward and reverse primers with a randomly generated insert in between.

### Barcode Validation Checks
- **Hamming Distance**: At least 140 out of 200 nucleotides must differ between any two accepted barcodes.
- **K-mer Check**: Using a sliding window (10-mers), no identical fragments should be present across different barcodes.
- **Internal Repeats**: The barcode must not contain repeated k-mer sequences within itself.
- **Partial Palindromes**: Sequences with semi-palindromic structures are rejected (e.g., a 6-nt fragment with a minimum loop of 3 nt).

### Final BLAST Verification
- The assembled barcodes are BLASTed against each other to ensure there is no undesired similarity.
- 
## Database Requirements

This script relies on a **BLAST database** of viral genome sequences. By default, the database is expected to be in the `genome_db/viral_sequences/` directory. If you wish to use a custom genome database, follow the steps below to create and set up the database:

### Creating the BLAST Database

1. **Download your genome sequences**: You will need a FASTA file containing the genome sequences (e.g., `genome_for_db.fasta`). These sequences can be viral or any other type depending on your needs. You can obtain genome sequences from various sources, such as [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/).

2. **Create the BLAST database**:
   To create the BLAST database from your FASTA file, use the `makeblastdb` command. This will index your genome sequences for use in BLAST searches.

   Run the following command in your terminal:

   ```bash
   makeblastdb -in genome_for_db.fasta -dbtype nucl -out genome_db/viral_sequences
## Configuration Parameters

Below are the main configuration settings you can modify in the script to customize the primer generation and barcode assembly process:

### General Configuration
```python
# === Configurations ===
N_TOTAL = 200  # Length of the final barcode (can be modified for different lengths)
INSERT_LENGTH = 160  # Length of the random insert between primers
N_PAIRS = 2  # Number of primer pairs to generate
GC_TARGET = 50  # Target GC content percentage for primer sequences
HAMMING_MIN = 140  # Minimum Hamming distance between barcodes
KMER_SIZE = 10  # Size of k-mers to check for duplicates in barcodes
output_dir = "work_repertory"  # Directory for saving results (change to desired output location)
BLAST_DB = os.path.join(output_dir, "genome_db", "viral_sequences")  # Path to the BLAST database
E_VALUE_CUTOFF = 0.1  # E-value threshold for BLAST search (controls similarity rejection)

# === Primer Generation Parameters ===
PRIMER_SIZE = 20  # (Fw & Rv) Primer size in nucleotides
PRIMER_MIN_SIZE = 18  # Minimum primer size (can be adjusted)
PRIMER_MAX_SIZE = 22  # Maximum primer size (can be adjusted)
PRIMER_OPT_TM = 60.0  # Optimal melting temperature (Tm) for primers (in °C)
PRIMER_MIN_TM = 57.0  # Minimum allowed melting temperature (Tm) for primers
PRIMER_MAX_TM = 63.0  # Maximum allowed melting temperature (Tm) for primers
PRIMER_NUM_RETURN = 1  # Number of primer pairs to return per primer design step
```

## Requirements

- Python 3.x
- `primer3-py` (for primer design)
- `Biopython` (for sequence manipulation and BLAST parsing)
- `BLAST+` command-line tools (for BLAST searches)
- `shutil` and `os` libraries (for file management)

### Installation

You can install the required Python dependencies using `pip` or `conda install`

```bash
pip install primer3-py biopython
```

Ensure that the **BLAST+** command-line tools are installed on your system. You can download them from [NCBI BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

Make sure that the BLAST tools are accessible from your system's PATH.

## Usage

1. **Run the Script**: Execute the script to generate primers, assemble barcodes, and perform the necessary checks. The following files will be generated:
   - `final_barcodes.fasta`: Contains the generated barcodes.
   - `primers.tsv`: Contains primer details (forward and reverse sequences, Tm, GC%, etc.).

2. **Temporary Files**: The script will generate temporary BLAST files (e.g., `temp_query.fasta`, `temp_result.xml`) but will clean them up automatically after processing.

3. **Customization**: Adjust the following parameters to fit your experimental needs:
   - Primer size, melting temperature, GC content, and BLAST E-value threshold.
   - Barcode length (currently set to 200 nt).
   - Internal repeat, palindrome, and Hamming distance check criteria.

4. **Clean-up**: The script automatically deletes unnecessary temporary files to maintain a clean working directory. (Only `result.xml` is kept)

## Example

```bash
python barcode_generator.py
```

This will run the script and generate the `final_barcodes.fasta` and `primers.tsv` files in the working directory.

