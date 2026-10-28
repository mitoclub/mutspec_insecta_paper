#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

RAW_DIR="${RAW_DIR:-${REPO_ROOT}/data/mammals_generation_length/raw}"
DATA_DIR="${DATA_DIR:-${REPO_ROOT}/data/mammals_generation_length}"
FIGURE_DIR="${FIGURE_DIR:-${REPO_ROOT}/figures}"

MIDORI="${MIDORI:-${RAW_DIR}/MIDORI2_LONGEST_NUC_GB270_Cytb_RAW.fasta.gz}"
GENERATION_LENGTH="${GENERATION_LENGTH:-${RAW_DIR}/Generation+Lenght+for+Mammals.xlsx}"
TREE="${TREE:-${RAW_DIR}/Mammalia_species.nwk}"

mkdir -p "${DATA_DIR}" "${FIGURE_DIR}"

python "${SCRIPT_DIR}/scripts/01_prepare_mammals.py" \
	--midori "${MIDORI}" \
	--generation-length "${GENERATION_LENGTH}" \
	--tree "${TREE}" \
	--output "${DATA_DIR}/mammals_pgls_dataset.csv" \
	--qc-report "${DATA_DIR}/preprocess_qc.txt"

python "${SCRIPT_DIR}/scripts/02_run_pgls.py" \
	--dataset "${DATA_DIR}/mammals_pgls_dataset.csv" \
	--tree "${TREE}" \
	--output "${DATA_DIR}/pgls_results.csv"

python "${SCRIPT_DIR}/scripts/03_plot_figure.py" \
	--dataset "${DATA_DIR}/mammals_pgls_dataset.csv" \
	--results "${DATA_DIR}/pgls_results.csv" \
	--pdf "${FIGURE_DIR}/Fig3A_mammals.pdf" \
	--png "${FIGURE_DIR}/Fig3A_mammals.png"

echo
echo "Done."
echo "Primary outputs:"
echo "  ${DATA_DIR}/mammals_pgls_dataset.csv"
echo "  ${DATA_DIR}/pgls_results.csv"
echo "  ${FIGURE_DIR}/Fig3A_mammals.pdf"
