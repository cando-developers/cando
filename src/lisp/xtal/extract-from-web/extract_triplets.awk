BEGIN { in_block = 0 }

/^[[:space:]]*x[[:space:]]*,/ { in_block = 1 }

/^[[:space:]]*$/ { in_block = 0 }

in_block && NF >= 1 && $0 ~ /,/ {
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
  print $0
}
