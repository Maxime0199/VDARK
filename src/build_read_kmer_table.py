#!/usr/bin/env python3
# =============================================================================
# build_read_kmer_table.py
# Associates each tumour-specific read with the tumour-specific k-mers it
# contains. Outputs a TSV used for read clustering in the R pipeline.
<<<<<<< HEAD
#
# NOTE: read_ID is suffixed with "/R1" or "/R2" to keep mates as distinct
# nodes in the downstream read-clustering graph. Without this, R1 and R2
# of the same fragment share the same QNAME and collapse into a single
# node, which can silently bridge two unrelated loci (e.g. if R1 carries
# a mutation at locus A and R2 carries a mutation at locus B) bypassing
# the MIN_SHARED_KMERS edge-weight filter entirely.
# =============================================================================
import csv
import argparse
import pyfastx
# ── Arguments ────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Build read-kmer association table")
parser.add_argument(
    "-k", "--kmer_size",
    type=int,
    default=31,
    help="k-mer size (default: 31)"
)
parser.add_argument(
    "-f", "--file_kmer",
    required=True,
    help="Input k-mer file"
)
parser.add_argument(
    "-o", "--output_file",
    required=True,
    help="Output file"
)
args = parser.parse_args()
K = args.kmer_size
KMER_FILE = args.file_kmer
OUTPUT_FILE = args.output_file
# ── Files ────────────────────────────────────────────────────────────────────
# Order matters: index 0 = R1, index 1 = R2 (used below to derive the mate tag)
READS = [
    f"rawdata/reads/tumour_R1_tumour_specific_chr21_k{K}.fq",
    f"rawdata/reads/tumour_R2_tumour_specific_chr21_k{K}.fq",
]
MATE_TAGS = ["R1", "R2"]
# ── Sequence utilities ───────────────────────────────────────────────────────
def rev_comp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
def canonical(kmer: str) -> str:
    rc = rev_comp(kmer)
    return kmer if kmer < rc else rc
def get_kmers(seq: str, k: int) -> list[str]:
    return [canonical(seq[i:i+k]) for i in range(len(seq) - k + 1)]
# ── Load tumour-specific k-mers ──────────────────────────────────────────────
print("[INFO] Loading tumour-specific k-mers...")
with open(KMER_FILE) as f:
    sig_kmers = {canonical(line.split()[0]) for line in f if line.strip()}
print(f"[INFO] {len(sig_kmers):,} k-mers loaded")
# ── Build read-kmer association table ────────────────────────────────────────
print("[INFO] Building read-kmer table...")
n_rows = 0
with open(OUTPUT_FILE, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["read_ID", "kmer"])
    for fq_file, mate_tag in zip(READS, MATE_TAGS):
        for name, seq, _ in pyfastx.Fastq(fq_file, build_index=False):
            hits = {k for k in get_kmers(seq, K) if k in sig_kmers}
            read_id = f"{name}/{mate_tag}"
            for kmer in hits:
                writer.writerow([read_id, kmer])
                n_rows += 1
=======
# =============================================================================

import csv
import argparse
import pyfastx


# ── Arguments ────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Build read-kmer association table")
parser.add_argument("-k", "--kmer_size", type=int, required=True,
                    help="k-mer size (e.g. 31)")
parser.add_argument("--wd", default="/srv/home/mlef0011/VDARK",
                    help="Working directory")

args = parser.parse_args()

K  = args.kmer_size
WD = args.wd


# ── Files ────────────────────────────────────────────────────────────────────
KMER_FILE   = f"{WD}/rawdata/kmer/tumour_specific_k{K}.txt"
OUTPUT_FILE = f"{WD}/rawdata/reads/read_kmer_association_k{K}.tsv"

READS = [
    f"{WD}/rawdata/reads/tumour_R1_tumour_specific_k{K}.fq",
    f"{WD}/rawdata/reads/tumour_R2_tumour_specific_k{K}.fq",
]


# ── Sequence utilities ───────────────────────────────────────────────────────
def rev_comp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def canonical(kmer: str) -> str:
    rc = rev_comp(kmer)
    return kmer if kmer < rc else rc


def get_kmers(seq: str, k: int) -> list[str]:
    return [canonical(seq[i:i+k]) for i in range(len(seq) - k + 1)]


# ── Load tumour-specific k-mers ──────────────────────────────────────────────
print("[INFO] Loading tumour-specific k-mers...")

with open(KMER_FILE) as f:
    sig_kmers = {canonical(line.split()[0]) for line in f if line.strip()}

print(f"[INFO] {len(sig_kmers):,} k-mers loaded")


# ── Build read-kmer association table ────────────────────────────────────────
print("[INFO] Building read-kmer table...")

n_rows = 0

with open(OUTPUT_FILE, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["read_ID", "kmer"])

    for fq_file in READS:
        for name, seq, _ in pyfastx.Fastq(fq_file, build_index=False):
            hits = {k for k in get_kmers(seq, K) if k in sig_kmers}

            for kmer in hits:
                writer.writerow([name, kmer])
                n_rows += 1

>>>>>>> 1d6bdf9a4bf4bd777b4dc1b9ca42945bc1278d5e
print(f"[INFO] Done — {n_rows:,} rows written to {OUTPUT_FILE}")