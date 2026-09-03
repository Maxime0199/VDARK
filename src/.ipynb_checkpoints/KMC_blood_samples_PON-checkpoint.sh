#!/bin/bash
set -euo pipefail

WD="/srv/home/mlef0011/VDARK"
TSV="/srv/home/mlef0011/VDARK/rawdata/PON/bam_paths_TCGA.tsv"
OUTDIR="/srv/home/mlef0011/VDARK/rawdata/PON/kmer"
BINDIR="/srv/home/mlef0011/VDARK/rawdata/PON/kmer"
TMPBASE="/srv/home/mlef0011/VDARK/rawdata/PON/tmp"
UNIONDIR="/srv/home/mlef0011/VDARK/rawdata/PON/kmer/union"
KMC="$WD/software/kmc/bin/kmc"
KMC_TOOLS="$WD/software/kmc/bin/kmc_tools"
THREADS=8

# ── Per-sample k-mer count threshold ──────────────────────────────────────
# -ci passed to KMC at counting time: minimum number of occurrences of a
# k-mer WITHIN a single BAM for it to be retained in that sample's DB.
# Default = 1 (present at least once). Override with -c <int>.
CI_SAMPLE=2

# ── Recurrence threshold across samples ───────────────────────────────────
# Minimum number of samples (after union) in which a k-mer must appear.
# Default = 2. Override with -r <int>.
CI_RECURRENCE=2

usage() {
    echo "Usage: $0 [-c CI_SAMPLE] [-r CI_RECURRENCE]" >&2
    echo "  -c  min k-mer occurrences WITHIN a sample to keep it (default: 1)" >&2
    echo "  -r  min number of SAMPLES a k-mer must appear in (default: 2)" >&2
    exit 1
}

while getopts "c:r:h" opt; do
    case "$opt" in
        c) CI_SAMPLE="$OPTARG" ;;
        r) CI_RECURRENCE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

echo "[$(date +%T)] Per-sample threshold: -ci${CI_SAMPLE}  |  Recurrence threshold: >=${CI_RECURRENCE} samples"

mkdir -p "$OUTDIR" "$BINDIR" "$TMPBASE" "$UNIONDIR"

tail -n +2 "$TSV" | while IFS=$'\t' read -r project normal_path tumour_path; do

    id=$(basename "$(dirname "$normal_path")")

    # ── File names now encode CI_SAMPLE so different thresholds never collide ─
    out_prefix="${OUTDIR}/${id}_k31_ci${CI_SAMPLE}"
    bin_prefix="${BINDIR}/${id}_k31_ci${CI_SAMPLE}_bin"
    tmpdir="${TMPBASE}/${id}"

    # ── Skip entirely if the compact (binarized) DB already exists ────────
    if [[ -f "${bin_prefix}.kmc_pre" && -f "${bin_prefix}.kmc_suf" ]]; then
        echo "[$(date +%T)] SKIP ${id} (${project}) — already compacted at ci${CI_SAMPLE}"
        continue
    fi

    # ── Count k-mers only if the raw DB doesn't already exist ─────────────
    if [[ -f "${out_prefix}.kmc_pre" && -f "${out_prefix}.kmc_suf" ]]; then
        echo "[$(date +%T)] REUSE ${id} (${project}) — raw k-mer DB already counted at ci${CI_SAMPLE}"
    else
        mkdir -p "$tmpdir"
        echo "[$(date +%T)] KMC on ${id} (${project}) -ci${CI_SAMPLE}"
        "$KMC" -k31 -fbam -ci"${CI_SAMPLE}" -t"${THREADS}" "$normal_path" "$out_prefix" "$tmpdir"
        #rm -rf "$tmpdir"
    fi

    # ── Compact (binarize): drop counters, keep only k-mer presence ───────
    echo "[$(date +%T)] Compacting ${id}"
    "$KMC_TOOLS" transform "$out_prefix" compact "$bin_prefix"

    # ── Keep only the compact file: remove the raw counted DB ─────────────
    #rm -f "${out_prefix}.kmc_pre" "${out_prefix}.kmc_suf"
done

# =============================================================================
# Union of all binarized DBs (kmc_tools complex) + recurrence filter
# =============================================================================

union_final="${UNIONDIR}/PON_union_ci${CI_SAMPLE}_recurrence_ci${CI_RECURRENCE}"
union_out="${UNIONDIR}/union_out_ci${CI_SAMPLE}"
op_file="${UNIONDIR}/union_op_ci${CI_SAMPLE}.txt"

if [[ -f "${union_final}.kmc_pre" && -f "${union_final}.kmc_suf" ]]; then
    echo "[$(date +%T)] SKIP union — recurrence-filtered PON already exists (ci${CI_SAMPLE}/rec${CI_RECURRENCE})"
else
    # ── Collect all binarized DBs present on disk for THIS ci threshold ────
    bin_prefixes=()
    while IFS= read -r -d '' f; do
        bin_prefixes+=("${f%.kmc_pre}")
    done < <(find "$BINDIR" -maxdepth 1 -name "*_k31_ci${CI_SAMPLE}_bin.kmc_pre" -print0)

    n=${#bin_prefixes[@]}
    if (( n < 2 )); then
        echo "[$(date +%T)] ERROR: need >=2 binarized DBs for a union, found ${n} (ci${CI_SAMPLE})" >&2
        exit 1
    fi

    echo "[$(date +%T)] Building union of ${n} binarized DB(s) at ci${CI_SAMPLE}"

    # ── Write kmc_tools complex operation file ──────────────────────────────
    {
        echo "INPUT:"
        for i in "${!bin_prefixes[@]}"; do
            echo "s$((i+1)) = ${bin_prefixes[$i]}"
        done
        echo "OUTPUT:"
        expr=$(seq 1 "$n" | sed 's/^/s/' | paste -sd+ -)
        echo "${union_out} = ${expr}"
    } > "$op_file"

    echo "[$(date +%T)] Running kmc_tools complex"
    "$KMC_TOOLS" complex "$op_file"

    # ── Recurrence filter: keep k-mers present in >=CI_RECURRENCE samples ───
    # + dump to txt in the same command
    echo "[$(date +%T)] Applying recurrence filter (-ci${CI_RECURRENCE}) + dumping txt"
    "$KMC_TOOLS" transform "$union_out" -ci"${CI_RECURRENCE}" \
        reduce "$union_final" \
        dump "${union_final}.txt"
fi

echo "[$(date +%T)] Done. Final PON DB: ${union_final}  (txt: ${union_final}.txt)"