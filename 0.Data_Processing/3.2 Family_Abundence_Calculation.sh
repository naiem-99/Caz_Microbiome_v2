
#!/bin/bash

set -euo pipefail
shopt -s nullglob

# ===============================
# LOAD PYTHON
# ===============================
module load python-3.9.18

# ===============================
# PATHS
# ===============================
IN_DIR="/data/pam/team216/ma32/scratch/metagenome/caz/bracken_Family_Bacteria_outputF"

SCRIPT="/nfs/users/nfs_m/ma32/tools/Bracken/analysis_scripts/combine_bracken_outputs.py"

OUT_FILE="${IN_DIR}/family_Bacteria_matrix.tsv"

# ===============================
# CHECK INPUT
# ===============================
cd "$IN_DIR"

FILES=( *.bracken )

if [ ${#FILES[@]} -eq 0 ]; then
    echo "No .bracken files found!"
    exit 1
fi

echo "Found ${#FILES[@]} Bracken files"

echo "First few files:"
printf '%s\n' "${FILES[@]:0:5}"

# ===============================
# RUN COMBINE
# ===============================
echo "Combining Bracken outputs..."

python3 "$SCRIPT" \
    --files "${FILES[@]}" \
    -o "$OUT_FILE"

echo "Matrix generated:"
echo "$OUT_FILE"
