#!/bin/bash
# =============================================================================
# get_local_minimum.sh
# Detects the first local minimum in a KMC k-mer frequency histogram.
# Used to set the coverage threshold for tumour-specific k-mer extraction.
# Usage: ./get_local_minimum.sh <histo_file>
# =============================================================================

if [ $# -ne 1 ]; then
    echo "Usage: $0 <histo_file>" >&2
    exit 1
fi

HISTO="$1"

if [ ! -f "$HISTO" ]; then
    echo "[ERROR] File not found: $HISTO" >&2
    exit 1
fi

LMIN=$(awk '
{
    count[NR] = $2
    freq[NR]  = $1
}
END {
    for (i = 2; i < NR; i++)
        smooth[i] = (count[i-1] + count[i] + count[i+1]) / 3

    for (i = 3; i < NR; i++) {
        if (smooth[i-1] < smooth[i-2] && smooth[i] > smooth[i-1]) {
            print freq[i-1]
            exit
        }
    }
}' "$HISTO")

if [ -z "$LMIN" ]; then
    echo "[ERROR] No local minimum detected in $HISTO" >&2
    exit 1
fi

echo "$LMIN"