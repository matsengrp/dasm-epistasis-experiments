#!/bin/bash
# Run partis cache-parameters on original and simulated MRCA FASTAs.
#
# Runs on pooled FASTAs (all subjects combined), one run per sequence type.
#
# Usage:
#   # Preview commands without executing:
#   bash run_partis_cache_parameters.sh --preview
#
#   # Check if previous runs completed:
#   bash run_partis_cache_parameters.sh --check
#
#   # Run both:
#   bash run_partis_cache_parameters.sh
#
#   # Run just one:
#   bash run_partis_cache_parameters.sh original_mrcas
#   bash run_partis_cache_parameters.sh simulated_mrcas

set -euo pipefail

# Fix terminfo for partis (colored_traceback needs it)
export TERMINFO=/usr/lib/terminfo
export TERM=${TERM:-xterm-256color}

INPUT_DIR="${INPUT_DIR:-_output/partis_input}"
OUTPUT_DIR="${OUTPUT_DIR:-_output/partis_output}"
N_PROCS="${N_PROCS:-8}"

ALL_SEQ_TYPES=("original_mrcas" "simulated_mrcas" "combined_mrcas")

# Parse flags
PREVIEW=false
CHECK=false
SEQ_TYPE_ARG=""

for arg in "$@"; do
    case "$arg" in
        --preview) PREVIEW=true ;;
        --check)   CHECK=true ;;
        *)         SEQ_TYPE_ARG="$arg" ;;
    esac
done

if [ -n "$SEQ_TYPE_ARG" ]; then
    SEQ_TYPES=("$SEQ_TYPE_ARG")
else
    SEQ_TYPES=("${ALL_SEQ_TYPES[@]}")
fi

MODE="RUN"
if $PREVIEW; then MODE="PREVIEW"; fi
if $CHECK; then MODE="CHECK"; fi

echo "[$MODE] partis cache-parameters"
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "N_PROCS: $N_PROCS"
echo ""

for SEQ_TYPE in "${SEQ_TYPES[@]}"; do
    INFNAME="$INPUT_DIR/${SEQ_TYPE}.fasta"
    PARAM_DIR="$OUTPUT_DIR/${SEQ_TYPE}_parameters"
    SW_CACHE="$OUTPUT_DIR/${SEQ_TYPE}_sw-cache.yaml"

    if [ ! -f "$INFNAME" ]; then
        echo "ERROR: $INFNAME not found"
        continue
    fi

    PARTIS_CMD=(partis cache-parameters \
        --infname "$INFNAME" \
        --parameter-dir "$PARAM_DIR" \
        --sw-cachefname "$SW_CACHE" \
        --locus igh \
        --n-procs "$N_PROCS" \
        --cache-hmm-annotations \
        --leave-default-germline \
        --random-seed 0)

    if $PREVIEW; then
        N_SEQS=$(grep -c "^>" "$INFNAME")
        echo "=== $SEQ_TYPE ==="
        echo "  ${PARTIS_CMD[*]}"
        echo "  ($N_SEQS sequences)"
        echo ""

    elif $CHECK; then
        HMM_CACHE="$PARAM_DIR/hmm-cache.yaml"
        echo -n "$SEQ_TYPE: "
        if [ -f "$HMM_CACHE" ] && [ -d "$PARAM_DIR/sw" ] && [ -d "$PARAM_DIR/hmm" ]; then
            echo "OK (hmm-cache.yaml exists)"
        elif [ -d "$PARAM_DIR" ]; then
            echo "INCOMPLETE (param dir exists but hmm-cache.yaml missing)"
        else
            echo "NOT RUN (no param dir)"
        fi

    else
        if [ -d "$PARAM_DIR" ]; then
            echo "SKIP: $PARAM_DIR already exists (use --check to verify completion)"
            continue
        fi

        echo "=== $SEQ_TYPE ==="
        echo "  Input: $INFNAME"
        echo "  Params: $PARAM_DIR"

        mkdir -p "$OUTPUT_DIR"
        "${PARTIS_CMD[@]}"

        echo "  Done: $SEQ_TYPE"
        echo ""
    fi
done

echo "[$MODE] All done."
