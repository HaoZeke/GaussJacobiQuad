#!/bin/bash

# Get the Git root directory
GITROOT=$(git rev-parse --show-toplevel)

# Function to print the header
print_header() {
  cat << EOF
! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! -----------------------------------------------------------------------------
! END_HEADER
EOF
}

# Iterate over all Fortran files in the SRC directory
for file in "${GITROOT}"/src/*.f90; do
  # Skip if no files found
  [ -e "$file" ] || continue

  # Same directory as the source, so the replace stays on one filesystem.
  temp_file=$(mktemp "$(dirname "$file")/.hdr.XXXXXX")

  # Check for the existing header and replace or append
  if grep -q 'BEGIN_HEADER' "$file"; then
    if ! grep -q 'END_HEADER' "$file"; then
      echo "warning: $file has BEGIN_HEADER but no END_HEADER; left unchanged" >&2
      rm -f "$temp_file"
      continue
    fi
    # Prefix, one new header, then the suffix a single time.
    # Keep any source glued onto the END_HEADER line.
    sed -n '1,/BEGIN_HEADER/{ /BEGIN_HEADER/q; p; }' "$file" > "$temp_file"
    print_header >> "$temp_file"
    awk '
      BEGIN { marker = "END_HEADER"; n = length(marker); found = 0 }
      found { print; next }
      {
        i = index($0, marker)
        if (i == 0) next
        found = 1
        tail = substr($0, i + n)
        sub(/\r$/, "", tail)
        if (length(tail) > 0) print tail
      }
    ' "$file" >> "$temp_file"
  else
    print_header > "$temp_file"
    cat "$file" >> "$temp_file"
  fi

  if ! chmod --reference="$file" "$temp_file" 2>/dev/null; then
    chmod "$(stat -c %a "$file" 2>/dev/null || stat -f %OLp "$file")" "$temp_file"
  fi
  mv "$temp_file" "$file"
done

echo "Headers added or updated in all Fortran files in $GITROOT/src."
