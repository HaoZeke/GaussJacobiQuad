#!/usr/bin/env python3

import argparse
import os
import re
from datetime import datetime

from jinja2 import Template

# Get the Git root directory
GITROOT = os.popen("git rev-parse --show-toplevel").read().strip()

# Get the current date
current_date = datetime.now().strftime("%Y-%m-%d")

# Get the latest commit hash
commit_hash = os.popen("git rev-parse HEAD").read().strip()[:7]

# Define the header template using Jinja2
header_template = Template(
    """! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: {{ author }}
! Source: {{ library_name }} Library
! License: MIT
! GitHub Repository: {{ repository_url }}
! Date: {{ date }}
! Commit: {{ commit }}
! -----------------------------------------------------------------------------
! This code is part of the {{ library_name }} library, providing an efficient
! implementation for Gauss-Jacobi quadrature nodes and weights computation.
! -----------------------------------------------------------------------------
! END_HEADER

"""
)

# Fill in the template variables
header_text = header_template.render(
    author="Rohit Goswami <rgoswami[at]ieee.org>",
    repository_url="https://github.com/HaoZeke/GaussJacobiQuad",
    library_name="GaussJacobiQuad",
    date=current_date,
    commit=commit_hash,
)


def add_headers(directories, file_types):
    for directory in directories:
        # Iterate over all files in the specified directory
        for root, dirs, files in os.walk(directory):
            for filename in files:
                if any(filename.endswith(suffix) for suffix in file_types):
                    filepath = os.path.join(root, filename)

                    with open(filepath, "r") as file:
                        content = file.read()

                    # Check for the existing header and replace or append
                    if "BEGIN_HEADER" in content and "END_HEADER" in content:
                        content = re.sub(
                            r"! BEGIN_HEADER.*?! END_HEADER\n",
                            header_text,
                            content,
                            flags=re.DOTALL,
                        )
                    else:
                        content = header_text + content

                    with open(filepath, "w") as file:
                        file.write(content)

        print(f"Headers added or updated in all files in {directory}.")


# Command-line arguments
parser = argparse.ArgumentParser(description="Add headers to source files.")
parser.add_argument(
    "--dirs", type=str, required=True, nargs="+", help="Directory(s) to process"
)
parser.add_argument(
    "--ftypes",
    type=str,
    required=True,
    help="Comma-separated list of file types to process",
)

args = parser.parse_args()

directories = args.dirs
file_types = args.ftypes.split(",")

add_headers(directories, file_types)
