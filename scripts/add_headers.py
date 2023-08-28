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
    """{{cchar}} BEGIN_HEADER
{{cchar}} -----------------------------------------------------------------------------
{{cchar}} Gauss-Jacobi Quadrature Implementation
{{cchar}} Authors: {{author}}
{{cchar}} Source: {{library_name}} Library
{{cchar}} License: MIT
{{cchar}} GitHub Repository: {{repository_url}}
{{cchar}} Date: {{date}}
{{cchar}} Commit: {{commit}}
{{cchar}} -----------------------------------------------------------------------------
{{cchar}} This code is part of the {{library_name}} library, providing an efficient
{{cchar}} implementation for Gauss-Jacobi quadrature nodes and weights computation.
{{cchar}} -----------------------------------------------------------------------------
{{cchar}} To cite this software:
{{cchar}} Rohit Goswami (2023). HaoZeke/GaussJacobiQuad: v0.1.0.
{{cchar}} Zenodo: https://doi.org/10.5281/ZENODO.8285112
{{cchar}} ---------------------------------------------------------------------
{{cchar}} END_HEADER

"""
)


def add_headers(
    directories,
    file_types,
    cchar="!",
    author="Rohit Goswami <rgoswami[at]ieee.org>",
    repository_url="https://github.com/HaoZeke/GaussJacobiQuad",
    library_name="GaussJacobiQuad",
    date=current_date,
    commit=commit_hash,
):
    # Fill in the template variables
    header_text = header_template.render(
        author=author,
        repository_url=repository_url,
        library_name=library_name,
        date=current_date,
        commit=commit_hash,
        cchar=cchar,
    )
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
                            rf"{cchar} BEGIN_HEADER.*?{cchar} END_HEADER\n",
                            header_text,
                            content,
                            flags=re.DOTALL,
                        )
                    else:
                        content = header_text + content

                    with open(filepath, "w") as file:
                        file.write(content)

        print("Headers added or updated in all files in" f" {directory} using {cchar}.")


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
parser.add_argument("--cchar", type=str, required=True, help="Comment character")

args = parser.parse_args()
file_types = args.ftypes.split(",")
add_headers(args.dirs, file_types, args.cchar)
