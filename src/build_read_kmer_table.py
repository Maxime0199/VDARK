#!/usr/bin/env python3
# =============================================================================
# build_read_kmer_table.py
# Associates each tumour-specific read with the tumour-specific k-mers it
# contains. Outputs a TSV used for read clustering in the R pipeline.
# Working directory: /srv/home/mlef0011/VDARK/
# =============================================================================

import csv
import pyfastx

K   = 31
WD  = "/srv/home/mlef0011/VDARK"

KMER_FILE   = f"{WD}/rawdata/kmer/tumour_specific.txt"
OUTPUT_FILE = f"{WD}/rawdata/reads/read_kmer_association.tsv"
READS = [
    f"{WD}/rawdata/reads/tumour_R1_tumour_specific.fq",
    f"{WD}/rawdata/reads/tumour_R2_tumour_specific.fq",
]


def rev_comp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

def canonical(kmer: str) -> str:
    rc = rev_comp(kmer)
    return kmer if kmer < rc else rc

def get_kmers(seq: str, k: int = K) -> list[str]:
    return [canonical(seq[i:i+k]) for i in range(len(seq) - k + 1)]


# Load tumour-specific k-mers (canonical form)
print("[INFO] Loading tumour-specific k-mers...")
with open(KMER_FILE) as f:
    sig_kmers = {canonical(line.split()[0]) for line in f if line.strip()}
print(f"[INFO] {len(sig_kmers):,} k-mers loaded")

# Build read-kmer association table
print("[INFO] Building read-kmer table...")
n_rows = 0
with open(OUTPUT_FILE, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["read_ID", "kmer"])

    for fq_file in READS:
        for name, seq, _ in pyfastx.Fastq(fq_file, build_index=False):
            hits = {k for k in get_kmers(seq) if k in sig_kmers}
            for kmer in hits:
                writer.writerow([name, kmer])
                n_rows += 1

print(f"[INFO] Done — {n_rows:,} rows written to {OUTPUT_FILE}")