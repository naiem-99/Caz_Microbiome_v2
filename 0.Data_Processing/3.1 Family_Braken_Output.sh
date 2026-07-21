#!/bin/bash

set -euo pipefail
shopt -s nullglob

# ===============================
# LOAD PYTHON
# ===============================
module load python-3.9.18

# Make "python" available for Bracken
mkdir -p ~/bin
ln -sf $(which python3) ~/bin/python
export PATH=~/bin:$PATH

# ===============================
# BRACKEN PATH
# ===============================
export PATH="/nfs/users/nfs_m/ma32/tools/Bracken:$PATH"

# ===============================
# PATHS
# ===============================
KRAKEN_DB="/data/pam/team216/ma32/scratch/metagenome/caz/kraken2_braken_db"

REPORT_DIR="/data/pam/team216/ma32/scratch/metagenome/caz/report_kraken2_bacteria"

OUT_DIR="/data/pam/team216/ma32/scratch/metagenome/caz/bracken_Family_Bacteria_outputF"

READ_LEN=150
LEVEL="F"

mkdir -p "${OUT_DIR}"

# ===============================
# RUN BRACKEN
# ===============================
for REPORT in "${REPORT_DIR}"/*.kreport; do

    SAMPLE=$(basename "${REPORT}" .kreport)

    echo "Running Bracken for ${SAMPLE} ..."

    bracken \
        -d "${KRAKEN_DB}" \
        -i "${REPORT}" \
        -o "${OUT_DIR}/${SAMPLE}.bracken" \
        -r "${READ_LEN}" \
        -l "${LEVEL}"

done

echo "Bracken Family abundance completed!"
