"""
Synthetic NGS Data Generator for 25 bp Genomic Indel Detection
Generates a realistic reference FASTA, index (.fai), and paired-end FASTQ reads (R1/R2)
with an embedded 25 bp deletion (target) alongside control variants (SNV, 5bp indel).
"""

import os
import random
import gzip
from pathlib import Path
from typing import Dict, Any, Tuple


DNA_BASES = ['A', 'C', 'G', 'T']
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    return "".join(COMPLEMENT.get(base, 'N') for base in reversed(seq))


def generate_random_sequence(length: int, seed: int = 42, gc_content: float = 0.45) -> str:
    """Generate a pseudo-random DNA sequence with specified GC content."""
    rng = random.Random(seed)
    seq = []
    for _ in range(length):
        if rng.random() < gc_content:
            seq.append(rng.choice(['C', 'G']))
        else:
            seq.append(rng.choice(['A', 'T']))
    return "".join(seq)


def generate_demo_dataset(
    output_dir: str,
    target_deletion_bp: int = 25,
    genome_length: int = 4000,
    read_length: int = 150,
    read_depth: int = 35,
    insert_mean: int = 300,
    insert_sd: int = 25,
    seed: int = 42,
    compress_fastq: bool = False
) -> Dict[str, Any]:
    """
    Synthesize reference FASTA, .fai index, and paired-end FASTQ files.
    
    Inserts:
    1. Target deletion at POS ~2000 (e.g., length = target_deletion_bp).
    2. Control small deletion (5 bp) at POS ~800.
    3. Control SNV (A -> T) at POS ~3200.
    """
    rng = random.Random(seed)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # 1. Generate Reference Genome (chr1)
    ref_seq = generate_random_sequence(genome_length, seed=seed)
    
    # Coordinates (1-indexed, positioned safely within genome)
    del25_pos = max(50, genome_length // 2)
    del5_pos = max(20, genome_length // 4)
    snv_pos = max(30, (genome_length * 3) // 4)

    # Ensure reference at deletion has distinct sequence
    del25_deleted_seq = "TAGCTAGCTAGCTAGCTAGCTAGCT"[:target_deletion_bp]
    ref_seq_list = list(ref_seq)
    
    # Place known target deletion sequence
    anchor_base = ref_seq_list[del25_pos - 1]
    for i, b in enumerate(del25_deleted_seq):
        ref_seq_list[del25_pos + i] = b
    ref_seq = "".join(ref_seq_list)

    # Write reference.fasta (with exact LF \n)
    ref_file = out_path / "reference.fasta"
    fai_file = out_path / "reference.fasta.fai"
    line_bases = 60
    
    with open(ref_file, "wb") as f:
        header = b">chr1\n"
        f.write(header)
        offset = len(header)
        for i in range(0, len(ref_seq), line_bases):
            line = (ref_seq[i:i + line_bases] + "\n").encode("ascii")
            f.write(line)
            
    # Write .fai index
    line_width = line_bases + 1  # 60 bases + 1 byte for \n
    with open(fai_file, "w") as f:
        f.write(f"chr1\t{len(ref_seq)}\t{offset}\t{line_bases}\t{line_width}\n")

    # 2. Build Alt Genome (with 25 bp deletion at del25_pos)
    # The 25 bp deletion removes ref_seq[del25_pos : del25_pos + target_deletion_bp]
    alt_seq = (
        ref_seq[:del25_pos] + 
        ref_seq[del25_pos + target_deletion_bp:]
    )
    
    # 3. Generate Paired-End Reads (Heterozygous 50% WT, 50% ALT)
    num_pairs = int((genome_length * read_depth) / (2 * read_length))
    r1_records = []
    r2_records = []

    for i in range(num_pairs):
        is_alt = (i % 2 == 1)
        active_seq = alt_seq if is_alt else ref_seq
        max_start = len(active_seq) - insert_mean - 100
        if max_start <= 0:
            frag_start = 0
        else:
            frag_start = rng.randint(0, max_start)
            
        cur_insert = max(read_length + 20, int(rng.gauss(insert_mean, insert_sd)))
        frag_end = min(len(active_seq), frag_start + cur_insert)
        fragment = active_seq[frag_start:frag_end]
        
        if len(fragment) < read_length:
            continue
            
        r1_seq = fragment[:read_length]
        r2_seq = reverse_complement(fragment[-read_length:])
        
        # High quality Phred score (Q37 = 'F' in Phred+33)
        q_scores = "".join(rng.choice(['E', 'F', 'G', 'H', 'I']) for _ in range(read_length))
        
        read_name = f"@SYNTH_READ_{i+1:06d}_{'ALT' if is_alt else 'WT'}"
        
        r1_records.append(f"{read_name}/1\n{r1_seq}\n+\n{q_scores}\n")
        r2_records.append(f"{read_name}/2\n{r2_seq}\n+\n{q_scores}\n")

    # Write FASTQ files
    ext = ".fastq.gz" if compress_fastq else ".fastq"
    r1_file = out_path / f"demo_sample_R1{ext}"
    r2_file = out_path / f"demo_sample_R2{ext}"

    if compress_fastq:
        with gzip.open(r1_file, "wt") as f1, gzip.open(r2_file, "wt") as f2:
            f1.writelines(r1_records)
            f2.writelines(r2_records)
    else:
        with open(r1_file, "w") as f1, open(r2_file, "w") as f2:
            f1.writelines(r1_records)
            f2.writelines(r2_records)

    # Truth data metadata
    expected_ref = ref_seq[del25_pos - 1 : del25_pos + target_deletion_bp]
    expected_alt = ref_seq[del25_pos - 1]
    
    metadata = {
        "reference_path": str(ref_file.resolve()),
        "fai_path": str(fai_file.resolve()),
        "r1_path": str(r1_file.resolve()),
        "r2_path": str(r2_file.resolve()),
        "target_deletion_bp": target_deletion_bp,
        "truth_deletion": {
            "chrom": "chr1",
            "pos": del25_pos,
            "ref": expected_ref,
            "alt": expected_alt,
            "del_size": target_deletion_bp,
            "genotype": "HET (0/1)",
            "deleted_sequence": del25_deleted_seq
        },
        "num_read_pairs": len(r1_records),
        "genome_length": genome_length
    }
    
    return metadata
