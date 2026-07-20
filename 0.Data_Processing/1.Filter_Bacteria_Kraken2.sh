#!/bin/bash

REPORT_DIR="/data/pam/team216/ma32/scratch/metagenome/caz/report_kraken2"

OUT_DIR="/data/pam/team216/ma32/scratch/metagenome/caz/report_kraken2_bacteria"

mkdir -p "$OUT_DIR"

for file in "$REPORT_DIR"/*.kreport; do

    base=$(basename "$file")

    echo "Processing $base"

    awk '
    BEGIN{bacteria=0}

    # Keep root
    $6=="R" {
        print
    }

    # Keep cellular organisms
    $6=="R1" {
        print
    }

    # Start bacterial subtree
    $6=="R2" && $7==2 {
        bacteria=1
    }

    # Stop at next major R2 subtree
    $6=="R2" && $7!=2 && bacteria==1 {
        exit
    }

    # Print bacterial subtree
    bacteria {
        print
    }

    ' "$file" > "$OUT_DIR/$base"

done

echo "Bacterial Bracken-compatible reports generated"
