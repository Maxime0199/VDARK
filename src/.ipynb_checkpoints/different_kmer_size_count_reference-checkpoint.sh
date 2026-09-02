#!/bin/bash
# =============================================================================
# count_kmers_ref_genome.sh
# Compte les k-mers du génome de référence pour un vecteur de valeurs de K.
# Produit, pour chaque K : un log KMC individuel.
# Produit, à la fin : un TSV récapitulatif (K, n_unique_kmers, n_total_kmers).
#
# Usage:
#   ./count_kmers_ref_genome.sh                # utilise K_VALUES ci-dessous
#   ./count_kmers_ref_genome.sh 21 25 31 41 51  # ou passe le vecteur en argument
# =============================================================================

set -euo pipefail

# ── Vecteur de K à parcourir — modifie directement ici si besoin ─────────────
K_VALUES=(11 13 15 17 19 21 23 25 27 29 31 33 35 37 40 45 50 55 60 65 70 75 80 85 90)
if [ $# -gt 0 ]; then
    K_VALUES=("$@")
fi

THREADS=12

WD="/srv/home/mlef0011/VDARK"
KMC="$WD/software/kmc/bin/kmc"
KMC_TOOLS="$WD/software/kmc/bin/kmc_tools"
REF="/srv/home/mlef0011/rawdata/ref_genome/reference_genome_GRCh37.fa"

OUT_DIR="$WD/output/count_ref_genome"
mkdir -p "$WD/tmp" "$OUT_DIR"

log() { echo "[$(date +%H:%M:%S)] $*"; }

SUMMARY_TSV="$OUT_DIR/summary_kmer_counts.tsv"
echo -e "K\tn_unique_kmers\tn_total_kmers" > "$SUMMARY_TSV"

for K in "${K_VALUES[@]}"; do

    if ! [[ "$K" =~ ^[0-9]+$ ]]; then
        log "[ERROR] K doit être un entier positif, reçu : $K -- ignoré"
        continue
    fi

    log "Comptage des k-mers, K=${K}..."

    LOG_FILE="$OUT_DIR/ref_genome_k${K}_kmc.log"

    "$KMC" -k${K} -t${THREADS} -m32 -ci1 -cs1000 -fm \
        "$REF" \
        "$OUT_DIR/ref_genome_k${K}_kmc" \
        "$WD/tmp" 2>&1 | tee "$LOG_FILE"

    N_UNIQUE=$(grep "No. of unique k-mers" "$LOG_FILE" | awk -F: '{print $2}' | tr -d ' ')
    N_TOTAL=$(grep "Total no. of k-mers"   "$LOG_FILE" | awk -F: '{print $2}' | tr -d ' ')

    if [ -z "$N_UNIQUE" ] || [ -z "$N_TOTAL" ]; then
        log "[ERROR] Impossible d'extraire les stats pour K=${K} -- verifie $LOG_FILE"
        continue
    fi

    echo -e "${K}\t${N_UNIQUE}\t${N_TOTAL}" >> "$SUMMARY_TSV"
    log "K=${K} -> ${N_UNIQUE} k-mers uniques / ${N_TOTAL} k-mers au total"

done

log "Termine. Recapitulatif dans $SUMMARY_TSV"