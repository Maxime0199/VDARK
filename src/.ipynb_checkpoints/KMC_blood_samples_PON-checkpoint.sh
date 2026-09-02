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

mkdir -p "$OUTDIR" "$BINDIR" "$TMPBASE" "$UNIONDIR"

tail -n +2 "$TSV" | while IFS=$'\t' read -r project normal_path tumour_path; do
    id=$(basename "$(dirname "$normal_path")")
    out_prefix="${OUTDIR}/${id}_k31"
    bin_prefix="${BINDIR}/${id}_k31_bin"
    tmpdir="${TMPBASE}/${id}"

    # ── Skip entirely if the compact (binarized) DB already exists ────────
    if [[ -f "${bin_prefix}.kmc_pre" && -f "${bin_prefix}.kmc_suf" ]]; then
        echo "[$(date +%T)] SKIP ${id} (${project}) — already compacted"
        continue
    fi

    # ── Count k-mers only if the raw DB doesn't already exist ─────────────
    if [[ -f "${out_prefix}.kmc_pre" && -f "${out_prefix}.kmc_suf" ]]; then
        echo "[$(date +%T)] REUSE ${id} (${project}) — raw k-mer DB already counted"
    else
        mkdir -p "$tmpdir"
        echo "[$(date +%T)] KMC on ${id} (${project})"
        "$KMC" -k31 -fbam -ci1 -t"${THREADS}" "$normal_path" "$out_prefix" "$tmpdir"
        rm -rf "$tmpdir"
    fi

    # ── Compact (binarize): drop counters, keep only k-mer presence ───────
    echo "[$(date +%T)] Compacting ${id}"
    "$KMC_TOOLS" transform "$out_prefix" compact "$bin_prefix"

    # ── Keep only the compact file: remove the raw counted DB ─────────────
    rm -f "${out_prefix}.kmc_pre" "${out_prefix}.kmc_suf"
done

# =============================================================================
# Union of all binarized DBs (kmc_tools complex) + recurrence filter (>=2)
# =============================================================================

union_raw="${UNIONDIR}/PON_union_raw"
union_final="${UNIONDIR}/PON_union_recurrence_ci2"
op_file="${UNIONDIR}/union_op.txt"

if [[ -f "${union_final}.kmc_pre" && -f "${union_final}.kmc_suf" ]]; then
    echo "[$(date +%T)] SKIP union — recurrence-filtered PON already exists"
else
    # ── Collect all binarized DBs present on disk ──────────────────────────
    bin_prefixes=()
    while IFS= read -r -d '' f; do
        bin_prefixes+=("${f%.kmc_pre}")
    done < <(find "$BINDIR" -maxdepth 1 -name '*_k31_bin.kmc_pre' -print0)

    n=${#bin_prefixes[@]}
    if (( n < 2 )); then
        echo "[$(date +%T)] ERROR: need >=2 binarized DBs for a union, found ${n}" >&2
        exit 1
    fi
    echo "[$(date +%T)] Building union of ${n} binarized DB(s)"

    # ── Write kmc_tools complex operation file ──────────────────────────────
    {
        echo "INPUT:"
        for i in "${!bin_prefixes[@]}"; do
            echo "s$((i+1)) = ${bin_prefixes[$i]}"
        done
        echo "OUTPUT:"
        expr=$(printf "+s%d" $(seq 1 "$n"))
        expr=${expr#+}   # drop leading '+'
        # rebuild with proper s1 + s2 + ... form
        expr=$(seq 1 "$n" | sed 's/^/s/' | paste -sd+ -)
        echo "union_out = ${expr}"
    } > "$op_file"

    echo "[$(date +%T)] Running kmc_tools complex"
    "$KMC_TOOLS" complex "$op_file"
    mv "${UNIONDIR}/union_out.kmc_pre" "${union_raw}.kmc_pre"
    mv "${UNIONDIR}/union_out.kmc_suf" "${union_raw}.kmc_suf"

    # ── Recurrence filter: keep k-mers present in >=2 samples ──────────────
    echo "[$(date +%T)] Applying recurrence filter (-ci2)"
    "$KMC_TOOLS" transform "$union_raw" -ci2 reduce "$union_final"

    rm -f "${union_raw}.kmc_pre" "${union_raw}.kmc_suf"
fi

echo "[$(date +%T)] Done. Final PON DB: ${union_final}"