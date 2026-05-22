#!/bin/bash

INPUT_DIR="parsed_triplets"
NAMES_FILE="spacegroup_names.txt"
OUTPUT_FILE="all_triplets.inp"

# Load SG number → name mapping
declare -A SGNAME
while read -r number && read -r name; do
  SGNAME[$number]="$name"
done < "$NAMES_FILE"

> "$OUTPUT_FILE"

# Sort files numerically by SG number
for file in $(ls "$INPUT_DIR"/sg*_triplets.txt | sed -E 's/.*sg([0-9]+)_triplets\.txt/\1 &/' | sort -n | cut -d' ' -f2-); do
  sgnum=$(basename "$file" | sed -E 's/^sg([0-9]+)_triplets\.txt$/\1/')
  sgname="${SGNAME[$sgnum]}"
  echo "# $sgnum ${sgname:-UNKNOWN}" >> "$OUTPUT_FILE"
  cat "$file" >> "$OUTPUT_FILE"
  echo "" >> "$OUTPUT_FILE"
  echo "" >> "$OUTPUT_FILE"
done

echo "Triplets combined with a blank line after each block in $OUTPUT_FILE"
