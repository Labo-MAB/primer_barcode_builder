import os
import random
import subprocess
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
import primer3
import numpy as np
import shutil

# === Configurations ===
N_TOTAL = 200
INSERT_LENGTH = 160
N_PAIRS = 2
GC_TARGET = 50
HAMMING_MIN = 140
KMER_SIZE = 10
output_dir = "work_repertory"
BLAST_DB = os.path.join(output_dir, "genome_db", "viral_sequences")
E_VALUE_CUTOFF = 0.1

# === Primer Generation Parameters ===
PRIMER_SIZE = 20 # (Fw & Rv)
PRIMER_MIN_SIZE = 18
PRIMER_MAX_SIZE = 22
PRIMER_OPT_TM = 60.0
PRIMER_MIN_TM = 57.0
PRIMER_MAX_TM = 63.0
PRIMER_NUM_RETURN = 1

# === Utility Functions ===
def _reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.
    """
    return str(Seq(seq).reverse_complement())

def _gc_content(seq):
    """
    Calculates the GC content percentage of a DNA sequence.
    """
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)

def _generate_seq_guided(length, gc_target):
    """
    Generates a random DNA sequence of specified length, 
    targeting a specified GC content percentage.
    """
    seq = []
    p_gc = gc_target / 100
    for i in range(length):
        allowed = list("ACGT")
        if i >= 3 and seq[-1] == seq[-2] == seq[-3]:
            allowed.remove(seq[-1])
        weights = [p_gc/2 if n in "GC" else (1-p_gc)/2 for n in allowed]
        nucleotide = random.choices(allowed, weights=weights, k=1)[0]
        seq.append(nucleotide)
    return "".join(seq)

def _hamming_distance(a, b):
    """
    Computes the Hamming distance between two sequences of equal length.
    Returns np.inf if sequences have different lengths.
    """
    if len(a) != len(b):
        return np.inf
    return sum(x != y for x, y in zip(a, b))

def _generate_kmers(seq, k):
    """
    Generates a set of k-mers from a given DNA sequence.
    """
    return {seq[i:i+k] for i in range(len(seq)-k+1)}

def _blast_sequence(seq, db=BLAST_DB, evalue=E_VALUE_CUTOFF):
    """
    Performs a BLAST search for a given DNA sequence against a specified database.
    Returns True if a significant alignment is found.
    """
    if not os.path.exists(f"{db}.nin"):
        raise FileNotFoundError(f"BLAST database not found: {db}")
    
    with open("temp_query.fasta", "w") as f:
        f.write(f">query\n{seq}\n")

    subprocess.run([ 
        "blastn", "-query", "temp_query.fasta", "-db", db,
        "-evalue", str(evalue), "-outfmt", "5", "-out", "temp_result.xml"
    ], check=True)

    with open("temp_result.xml") as result:
        blast_records = list(NCBIXML.parse(result))
        for record in blast_records:
            for align in record.alignments:
                for hsp in align.hsps:
                    print("   Rejected by BLAST on external DB:")
                    print(f"  Sequence: {seq}")
                    print(f"  Aligned with: {align.hit_def}")
                    print(f"  E-value: {hsp.expect:.2e}")
                    print(f"  Score: {hsp.score}")
                    return True 
    return False

def _has_internal_repeats(seq, repeat_length=10):
    """
    Returns True if a subsequence of length 'repeat_length' appears more than once in 'seq' (including overlaps).
    """
    counts = {}
    for i in range(len(seq) - repeat_length + 1):
        subseq = seq[i:i+repeat_length]
        counts[subseq] = counts.get(subseq, 0) + 1
        if counts[subseq] > 1:
            return True
    return False

def _has_partial_palindrome(seq, stem_length=6, loop_size=3):
    """
    Returns True if a hairpin-like structure is found in 'seq'.
    It searches for a segment of 'stem_length' (stem1) followed, after a loop of at least 'loop_size' nucleotides,
    by a segment that is the reverse complement of stem1.
    """
    seq_len = len(seq)
    for i in range(seq_len - (2 * stem_length + loop_size) + 1):
        stem1 = seq[i : i + stem_length]
        for j in range(i + stem_length + loop_size, seq_len - stem_length + 1):
            stem2 = seq[j : j + stem_length]
            if stem1 == _reverse_complement(stem2):
                return True
    return False

# ===================================================
# Similarity check between accepted primers via BLAST
# ===================================================
def blast_against_accepted(primer, accepted_primers, evalue=E_VALUE_CUTOFF):
    """
    Blast the given primer against a temporary database built from accepted_primers.
    Returns True if an alignment (other than self) is found.
    """
    if not accepted_primers:
        return False

    temp_db_fasta = "temp_accepted_primers.fasta"
    # Writing the accepted primers to a temporary FASTA file
    with open(temp_db_fasta, "w") as f:
        for i, seq in enumerate(accepted_primers):
            f.write(f">accepted_{i}\n{seq}\n")
    
    # Creating the BLAST database from the temporary file
    makeblastdb_command = ["makeblastdb", "-in", temp_db_fasta, "-dbtype", "nucl", "-out", "temp_accepted_db"]
    subprocess.run(makeblastdb_command, check=True)
    
    # Writing the candidate primer to a temporary FASTA file
    with open("temp_query.fasta", "w") as f:
        f.write(f">query\n{primer}\n")
    
    # Running BLAST on the primer against the accepted primers DB
    blastn_command = ["blastn", "-query", "temp_query.fasta", "-db", "temp_accepted_db", "-evalue", str(evalue), "-outfmt", "5", "-out", "temp_result.xml"]
    subprocess.run(blastn_command, check=True)
    
    with open("temp_result.xml") as result:
        blast_records = list(NCBIXML.parse(result))
        for record in blast_records:
            for align in record.alignments:
                for hsp in align.hsps:
                    print("   Rejected by BLAST between primers:")
                    print(f"  Candidate primer: {primer}")
                    print(f"  Aligned with: {align.hit_def}")
                    print(f"  E-value: {hsp.expect:.2e}")
                    print(f"  Score: {hsp.score}")
                    return True
    return False

# ===================================================
# 1. Primer Generation with Additional Check Between Primers
# ===================================================
def generate_primers(n):
    """
    Generates 'n' pairs of primers with additional checks such as GC content, BLAST against external databases, 
    and complementarity to previously accepted primers.
    """
    primer_pairs = []
    primer_infos = []
    accepted_primers = []  # List to store all accepted primers (forward and reverse)
    while len(primer_pairs) < n:
        seq_template = _generate_seq_guided(INSERT_LENGTH, GC_TARGET)
        result = primer3.bindings.design_primers({
            'SEQUENCE_TEMPLATE': seq_template
        }, {
            'PRIMER_OPT_SIZE': PRIMER_SIZE,
            'PRIMER_MIN_SIZE': PRIMER_MIN_SIZE,
            'PRIMER_MAX_SIZE': PRIMER_MAX_SIZE,
            'PRIMER_OPT_TM': PRIMER_OPT_TM,
            'PRIMER_MIN_TM': PRIMER_MIN_TM,
            'PRIMER_MAX_TM': PRIMER_MAX_TM,
            'PRIMER_NUM_RETURN': PRIMER_NUM_RETURN
        })
        fwd = result.get('PRIMER_LEFT_0_SEQUENCE')
        rev = result.get('PRIMER_RIGHT_0_SEQUENCE')
        tm_fwd = result.get('PRIMER_LEFT_0_TM')
        tm_rev = result.get('PRIMER_RIGHT_0_TM')
        if not fwd or not rev:
            continue
        # Check GC content between 45% and 55%
        if not (45 <= _gc_content(fwd) <= 55) or not (45 <= _gc_content(rev) <= 55):
            continue
        # BLAST check on viral genome for each primer
        if _blast_sequence(fwd) or _blast_sequence(rev):
            continue

        # --- Checking complementarity between accepted primers ---
        if blast_against_accepted(fwd, accepted_primers) or blast_against_accepted(rev, accepted_primers):
            continue
        
        # If the candidate passes all tests, add to the list
        primer_pairs.append((fwd, rev))
        primer_infos.append((
            f"primer_{len(primer_pairs)-1}",
            fwd,
            tm_fwd,
            rev,
            tm_rev,
            len(fwd),
            _gc_content(fwd),
            len(rev),
            _gc_content(rev)
        ))
        # Add accepted primers for future checks
        accepted_primers.append(fwd)
        accepted_primers.append(rev)
    return primer_pairs, primer_infos


# ===================================================
# 2. Barcode Assembly
# ===================================================

def assemble_barcodes(primer_pairs):
    """
    Assembles barcodes by adding random inserts between the forward and reverse primers.
    """
    full_barcodes = []
    for fwd, rev in primer_pairs:
        insert = _generate_seq_guided(N_TOTAL - len(fwd) - len(rev), GC_TARGET)
        full_seq = fwd + insert + _reverse_complement(rev)
        full_barcodes.append(full_seq)
    return full_barcodes

# ======================================================
# BLAST final on generated barcodes (second verification)
# ======================================================

def move_and_clean_blast_files(output_dir):
    """
    Removes unnecessary files in the target directory and moves BLAST files.
    This function cleans up extra files (like `results.xml` and `barcode_db`) 
    and moves the BLAST output files to the target directory.
    """
    # List of files to remove from the output directory
    files_to_delete = [
        "barcode_db.ndb",         # BLAST database file
        "barcode_db.not",         # BLAST database file
        "barcode_db.ntf",         # BLAST database file
        "barcode_db.nto",         # BLAST database file
        "temp_query.fasta",       # Temporary query file
        "temp_result.xml",        # Temporary BLAST result file
    ]
        
    for file in files_to_delete:
        file_path = os.path.join(output_dir, file)  # Files are in `output_dir`
        try:
            if os.path.exists(file_path): 
                os.remove(file_path)
        except FileNotFoundError:
            pass


def run_blast_commands(output_dir):
    """
    Creates a BLAST database from the FASTA file and runs BLAST to verify the generated barcodes.
    """
    fasta_path = os.path.join(output_dir, "final_barcodes.fasta")
    db_prefix = os.path.join(output_dir, "barcode_db")
    results_path = os.path.join(output_dir, "results.xml")
    
    try:
        # Create BLAST database
        makeblastdb_command = ["makeblastdb", "-in", fasta_path, "-dbtype", "nucl", "-out", db_prefix]
        subprocess.run(makeblastdb_command, check=True)
        
        # Run BLAST
        blastn_command = ["blastn", "-query", fasta_path, "-db", db_prefix, "-out", results_path, "-evalue", "1e-1", "-outfmt", "5"]
        subprocess.run(blastn_command, check=True)
        print("Second BLAST verification completed.")
    finally:
        temp_files = [
            "temp_query.fasta",
            "temp_result.xml",
            f"{db_prefix}.nin",
            f"{db_prefix}.nsq",
            f"{db_prefix}.nhr"
        ]
        for file in temp_files:
            try:
                os.remove(file)
            except FileNotFoundError:
                pass



def parse_blast_results(output_dir):
    """
    Parses and displays the results of the BLAST search from the generated XML file.
    """
    results_path = os.path.join(output_dir, "results.xml")
    with open(results_path) as result_handle:
        blast_records = list(NCBIXML.parse(result_handle))
        print(f"Number of BLAST records: {len(blast_records)}")
        for blast_record in blast_records:
            query_id = blast_record.query
            non_self_alignments = []
            for alignment in blast_record.alignments:
                if query_id in alignment.hit_def:
                    continue
                for hsp in alignment.hsps:
                    non_self_alignments.append((alignment.hit_def, hsp.expect, hsp.bits))
            if non_self_alignments:
                for hit in non_self_alignments:
                    print(f"Query: {query_id}")
                    print(f"Subject: {hit[0]}")
                    print(f"E-value: {hit[1]}")
                    print(f"Bitscore: {hit[2]}\n")
            else:
                print(f"{query_id}: no significant alignment.")

# ===================================================
# MAIN
# ===================================================
def main():
    """
    Main function to generate primers, assemble barcodes, check for internal repeats and palindromes,
    and perform BLAST verification on the final set of barcodes.
    """
    final_barcodes = []
    final_infos = []
    kmers_list = []
    attempts = 0

    while len(final_barcodes) < N_PAIRS and attempts < N_PAIRS * 200:
        primers, primer_infos = generate_primers(1)
        barcodes = assemble_barcodes(primers)
        bc = barcodes[0]
        info = primer_infos[0]

        # Hamming distance check between barcodes
        if any(_hamming_distance(bc, other) < HAMMING_MIN for other in final_barcodes):
            attempts += 1
            continue

        # K-mer check between barcodes
        kmers = _generate_kmers(bc, KMER_SIZE)
        if any(kmers & existing for existing in kmers_list):
            attempts += 1
            continue

        # Internal repeats check
        if _has_internal_repeats(bc, repeat_length=KMER_SIZE):
            attempts += 1
            continue

        # Partial palindrome (hairpin) check
        if _has_partial_palindrome(bc, stem_length=6, loop_size=3):
            attempts += 1
            continue

        # Final BLAST check on individual barcode
        if _blast_sequence(bc, db=BLAST_DB, evalue=1e-1):
            attempts += 1
            continue

        # Update primer ID to be globally unique
        info = list(info)
        info[0] = f"primer_{len(final_barcodes)}"
        info = tuple(info)

        final_barcodes.append(bc)
        final_infos.append(info)
        kmers_list.append(kmers)
        attempts += 1
        print(f"Attempt {attempts}: {len(final_barcodes)} valid barcodes")

    save_primers_to_table(final_infos[:N_PAIRS], output_dir)
    save_barcodes_to_fasta(final_barcodes[:N_PAIRS], output_dir)
    print(f"\n{len(final_barcodes[:N_PAIRS])} final barcodes saved in '{output_dir}'!")

    run_blast_commands(output_dir)
    parse_blast_results(output_dir)
    move_and_clean_blast_files(output_dir)

# ===================================================
# Export functions
# ===================================================

def save_primers_to_table(primer_infos, output_dir, filename="primers.tsv"):
    """
    Saves the primer information (forward and reverse sequences, Tm, lengths, GC%) to a TSV file.
    """
    os.makedirs(output_dir, exist_ok=True)
    filepath = os.path.join(output_dir, filename)
    with open(filepath, mode='w') as file:
        file.write("primer\tforward (sequence)\tTm (fw)\treverse (sequence)\tTm (rv)\tlen_fw\tGC%_fw\tlen_rv\tGC%_rv\n")
        for primer_id, fwd, tm_fwd, rev, tm_rev, len_fw, gc_fw, len_rv, gc_rv in primer_infos:
            file.write(f"{primer_id}\t{fwd}\t{tm_fwd:.2f}\t{rev}\t{tm_rev:.2f}\t{len_fw}\t{gc_fw:.2f}\t{len_rv}\t{gc_rv:.2f}\n")

def save_barcodes_to_fasta(barcodes, output_dir, filename="final_barcodes.fasta"):
    """
    Saves the generated barcodes to a FASTA file, each barcode with a unique identifier.
    """
    os.makedirs(output_dir, exist_ok=True)
    filepath = os.path.join(output_dir, filename)
    with open(filepath, "w") as f:
        for i, seq in enumerate(barcodes):
            f.write(f">barcode_{i}\n{seq}\n")

if __name__ == '__main__':
    main()
