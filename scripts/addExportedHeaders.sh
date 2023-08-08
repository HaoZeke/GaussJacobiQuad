#!/bin/bash

# Get the Git root directory
GITROOT=$(git rev-parse --show-toplevel)

# Path to the header file
HEADER="$GITROOT/scripts/export_header.txt"

# Iterate over all Fortran files in the SRC directory
for file in "${GITROOT}"/src/*.f90; do
  # Skip if no files found
  [ -e "$file" ] || continue

  # Create a temporary file
  temp_file=$(mktemp)

  # Concatenate the header and the original file
  cat "$HEADER" "$file" > "$temp_file"

  # Replace the original file with the new one
  mv "$temp_file" "$file"
done

echo "Headers added to all Fortran files in $GITROOT/src."
