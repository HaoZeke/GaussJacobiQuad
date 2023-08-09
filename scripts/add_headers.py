#!/usr/bin/env python3

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

# Iterate over all Fortran files in the SRC directory
for root, dirs, files in os.walk(os.path.join(GITROOT, "src")):
    for filename in files:
        if filename.endswith(".f90"):
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

print("Headers added or updated in all Fortran files in {}/src.".format(GITROOT))
