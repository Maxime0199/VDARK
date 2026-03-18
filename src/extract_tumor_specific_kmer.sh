#!/bin/bash
# =============================================================================
# extract_tumor_specific_kmers.sh
# Counts k-mers in tumour and normal samples, then extracts tumour-specific
# k-mers (present in tumour above coverage threshold, absent in normal).
# Filters tumour reads to retain only those carrying tumour-specific k-mers.
# Working directory: /srv/home/mlef0011/VDARK/
# =============================================================================

set -euo pipefail


K=31
THREADS=12


WD="/srv/home/mlef0011/VDARK"
KMC="$WD/software/kmc/bin/kmc"
KMC_TOOLS="$WD/software/kmc/bin/kmc_tools"
GET_MIN="$WD/src/get_local_minimum.sh"
SUFFIX="k${K}"


log() { echo "[$(date +%H:%M:%S)] $*"; }

# ── Count k-mers ──────────────────────────────────────────────────────────────

log "Counting tumour k-mers..."
"$KMC" -k${K} -t${THREADS} -ci1 -cs1000 -fq @"$WD/rawdata/reads/tumour_fastq.txt" \
    "$WD/rawdata/kmer/tumour_${SUFFIX}_kmc" "$WD/tmp"

"$KMC" -k${K} -t${THREADS} -ci1 -cs1000 -fq @"$WD/rawdata/reads/normal_fastq.txt" \
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

# ── Extract tumour-specific k-mers ────────────────────────────────────────────

log "Detecting coverage threshold..."
MIN_T=$($GET_MIN $WD/rawdata/kmer/tumour.histo)
log "Coverage threshold: $MIN_T"

log "Subtracting normal k-mers..."
"$KMC_TOOLS" simple \
    "$WD/rawdata/kmer/tumour_${SUFFIX}_kmc" -ci"$MIN_T" \
    "$WD/rawdata/kmer/normal_${SUFFIX}_kmc" -ci2 \
    kmers_subtract \
    "$WD/rawdata/kmer/tumour_specific_${SUFFIX}"

"$KMC_TOOLS" transform "$WD/rawdata/kmer/tumour_specific_${SUFFIX}" dump \
    "$WD/rawdata/kmer/tumour_specific_${SUFFIX}.txt"
    
# ── Filter tumour-specific reads ──────────────────────────────────────────────

log "Filtering tumour-specific reads..."
"$KMC_TOOLS" filter "$WD/rawdata/kmer/tumour_specific_${SUFFIX}" -ci1 -cx1000 \
    "$WD/rawdata/reads/tumour_R1.fq" \
    "$WD/rawdata/reads/tumour_R1_tumour_specific_${SUFFIX}.fq"

"$KMC_TOOLS" filter "$WD/rawdata/kmer/tumour_specific_${SUFFIX}" -ci1 -cx1000 \
    "$WD/rawdata/reads/tumour_R2.fq" \
    "$WD/rawdata/reads/tumour_R2_tumour_specific_${SUFFIX}.fq"

log "Done."