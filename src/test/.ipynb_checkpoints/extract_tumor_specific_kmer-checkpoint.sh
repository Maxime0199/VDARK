#!/bin/bash
# =============================================================================
# extract_tumor_specific_kmers.sh
# Counts k-mers in tumour and normal samples, then extracts tumour-specific
# k-mers using a 2-tier germline filter (normal=0 auto-keep, normal 1-3
# via binomial test), and filters tumour reads accordingly.
# Working directory: /srv/home/mlef0011/VDARK/
# =============================================================================
set -euo pipefail

K=31
THREADS=12
NORMAL_GREY_MAX=3
P_THRESHOLD=0.01

WD="/srv/home/mlef0011/VDARK"
KMC="$WD/software/kmc/bin/kmc"
KMC_TOOLS="$WD/software/kmc/bin/kmc_tools"
RSCRIPT="Rscript"
FILTER_R="$WD/src/test/filter_tumour_specific_kmers.r"
SUFFIX="chr21_k${K}"

log() { echo "[$(date +%H:%M:%S)] $*"; }

# ── Count k-mers ──────────────────────────────────────────────────────────────
log "Counting tumour k-mers..."
"$KMC" -k${K} -t${THREADS} -ci1 -cs1000 -fq @"$WD/rawdata/reads/tumour_fastq_chr21.txt" \
    "$WD/rawdata/kmer/tumour_${SUFFIX}_kmc" "$WD/tmp"

log "Counting normal k-mers..."
"$KMC" -k${K} -t${THREADS} -ci1 -cs1000 -fq @"$WD/rawdata/reads/normal_fastq_chr21.txt" \
    "$WD/rawdata/kmer/normal_${SUFFIX}_kmc" "$WD/tmp"

# ── Histograms & dumps ────────────────────────────────────────────────────────
log "Generating histograms..."
"$KMC_TOOLS" transform "$WD/rawdata/kmer/tumour_${SUFFIX}_kmc" histogram \
    "$WD/rawdata/kmer/tumour_${SUFFIX}.histo" -cx1000
"$KMC_TOOLS" transform "$WD/rawdata/kmer/normal_${SUFFIX}_kmc" histogram \
    "$WD/rawdata/kmer/normal_${SUFFIX}.histo" -cx1000

"$KMC_TOOLS" transform "$WD/rawdata/kmer/tumour_${SUFFIX}_kmc" dump \
    "$WD/rawdata/kmer/tumour_${SUFFIX}_kmer.txt"
"$KMC_TOOLS" transform "$WD/rawdata/kmer/normal_${SUFFIX}_kmc" dump \
    "$WD/rawdata/kmer/normal_${SUFFIX}_kmer.txt"

# ── Tier 1 : normal=0, tumour>=2 -> garde direct ──────────────────────────────
log "Subtracting normal > 0 - tumour >= 2 (Tier 1)"
"$KMC_TOOLS" simple \
    "$WD/rawdata/kmer/tumour_${SUFFIX}_kmc" -ci3 \
    "$WD/rawdata/kmer/normal_${SUFFIX}_kmc" -ci1 \
    kmers_subtract \
    "$WD/rawdata/kmer/tumour_specific_filt1_${SUFFIX}"

"$KMC_TOOLS" transform "$WD/rawdata/kmer/tumour_specific_filt1_${SUFFIX}" dump \
    "$WD/rawdata/kmer/tumour_specific_filt1_${SUFFIX}.txt"

# ── Tier 2 : normal 1-3, tumour>=2 -> zone grise, comptes des deux côtés ─────
log "Extracting grey zone -- tumour counts..."
"$KMC_TOOLS" simple \
    "$WD/rawdata/kmer/tumour_${SUFFIX}_kmc" -ci3 \
    "$WD/rawdata/kmer/normal_${SUFFIX}_kmc" -ci1 -cx${NORMAL_GREY_MAX} \
    intersect \
    "$WD/rawdata/kmer/tumour_specific_grey_zone_tumour_count_${SUFFIX}" -ocleft

"$KMC_TOOLS" transform "$WD/rawdata/kmer/tumour_specific_grey_zone_tumour_count_${SUFFIX}" dump \
    "$WD/rawdata/kmer/tumour_specific_grey_zone_tumour_count_${SUFFIX}.txt"

log "Extracting grey zone -- normal counts..."
"$KMC_TOOLS" simple \
    "$WD/rawdata/kmer/tumour_${SUFFIX}_kmc" -ci3 \
    "$WD/rawdata/kmer/normal_${SUFFIX}_kmc" -ci1 -cx${NORMAL_GREY_MAX} \
    intersect \
    "$WD/rawdata/kmer/tumour_specific_grey_zone_normal_count_${SUFFIX}" -ocright

"$KMC_TOOLS" transform "$WD/rawdata/kmer/tumour_specific_grey_zone_normal_count_${SUFFIX}" dump \
    "$WD/rawdata/kmer/tumour_specific_grey_zone_normal_count_${SUFFIX}.txt"

# ── Filtre statistique (R) : combine Tier 1 + Tier 2 significatif ───────────
log "Running statistical filter in R..."
"$RSCRIPT" "$FILTER_R" \
    --histo_tumour "$WD/rawdata/kmer/tumour_${SUFFIX}.histo" \
    --histo_normal "$WD/rawdata/kmer/normal_${SUFFIX}.histo" \
    --tier1        "$WD/rawdata/kmer/tumour_specific_filt1_${SUFFIX}.txt" \
    --grey_tumour  "$WD/rawdata/kmer/tumour_specific_grey_zone_tumour_count_${SUFFIX}.txt" \
    --grey_normal  "$WD/rawdata/kmer/tumour_specific_grey_zone_normal_count_${SUFFIX}.txt" \
    --p_threshold  "$P_THRESHOLD" \
    --out          "$WD/rawdata/kmer/tumour_specific_${SUFFIX}_filtered.txt"

N_KEPT=$(wc -l < "$WD/rawdata/kmer/tumour_specific_${SUFFIX}_filtered.txt")
log "${N_KEPT} k-mer(s) retained after 2-tier germline filter"

# ── Reconstruit une base KMC à partir de la liste finale ─────────────────────
log "Rebuilding KMC database from filtered k-mer list..."
FILTERED_FA="$WD/rawdata/kmer/tumour_specific_${SUFFIX}_filtered.fa"
awk '{print ">k_"NR; print $1}' "$WD/rawdata/kmer/tumour_specific_${SUFFIX}_filtered.txt" > "$FILTERED_FA"

"$KMC" -k${K} -t${THREADS} -ci1 -fm "$FILTERED_FA" \
    "$WD/rawdata/kmer/tumour_specific_${SUFFIX}_filtered" "$WD/tmp"

# ── Filter tumour-specific reads ──────────────────────────────────────────────
log "Filtering tumour-specific reads..."
"$KMC_TOOLS" filter "$WD/rawdata/kmer/tumour_specific_${SUFFIX}_filtered" -ci1 -cx1000 \
    "$WD/rawdata/reads/tumour_chr21_R1.fq" \
    "$WD/rawdata/reads/tumour_R1_tumour_specific_${SUFFIX}_filtered.fq"
"$KMC_TOOLS" filter "$WD/rawdata/kmer/tumour_specific_${SUFFIX}_filtered" -ci1 -cx1000 \
    "$WD/rawdata/reads/tumour_chr21_R2.fq" \
    "$WD/rawdata/reads/tumour_R2_tumour_specific_${SUFFIX}_filtered.fq"

log "Done."