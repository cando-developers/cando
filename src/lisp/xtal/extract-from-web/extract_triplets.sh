#!/bin/bash

INPUT_DIR="symmetry_triplets"
OUTPUT_DIR="parsed_triplets"
AWK_SCRIPT="extract_triplets.awk"

mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/sg*.txt; do
  sgnum=$(basename "$file" .txt | sed 's/sg//')
  outfile="$OUTPUT_DIR/sg${sgnum}_triplets.txt"
  echo "Parsing space group $sgnum..."

  awk -f "$AWK_SCRIPT" "$file" > "$outfile"
done

echo "Triplet extraction complete. See ./parsed_triplets/"
