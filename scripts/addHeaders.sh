#!/bin/bash

# Get the Git root directory
GITROOT=$(git rev-parse --show-toplevel)

# Function to print the header
print_header() {
  cat << EOF
# BEGIN_HEADER
# -----------------------------------------------------------------------------
# Gauss-Jacobi Quadrature Implementation
# Authors: Rohit Goswami (rog32[at]hi[dot]is)
# Source: GaussJacobiQuad Library
# License: MIT
# GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
# -----------------------------------------------------------------------------
# This code is part of the GaussJacobiQuad library, providing an efficient
# implementation for Gauss-Jacobi quadrature nodes and weights computation.
# -----------------------------------------------------------------------------
# END_HEADER
#
EOF
}

# Iterate over all Fortran files in the SRC directory
for file in "${GITROOT}"/src/*.f90; do
  # Skip if no files found
  [ -e "$file" ] || continue

  # Create a temporary file
  temp_file=$(mktemp)

  # Check for the existing header and replace or append
  if grep -q 'BEGIN_HEADER' "$file"; then
    # Print everything before BEGIN_HEADER to the temporary file
    sed -n '/BEGIN_HEADER/!p' "$file" > "$temp_file"

    # Print the new header to the temporary file
    print_header >> "$temp_file"

    # Print everything after END_HEADER to the temporary file
    sed -n '/END_HEADER/,$p' "$file" | tail -n +2 >> "$temp_file"
  else
    # Print the new header to the temporary file
    print_header > "$temp_file"

    # Concatenate the original file
    cat "$file" >> "$temp_file"
  fi

  # Replace the original file with the new one
  mv "$temp_file" "$file"
done

echo "Headers added or updated in all Fortran files in $GITROOT/src."
