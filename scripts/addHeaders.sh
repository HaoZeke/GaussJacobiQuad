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
    # Prefix, one new header, then the suffix a single time.
    sed -n '1,/BEGIN_HEADER/{ /BEGIN_HEADER/q; p; }' "$file" > "$temp_file"
    print_header >> "$temp_file"
    sed -n '/END_HEADER/,$p' "$file" | tail -n +2 >> "$temp_file"
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
