#!/usr/bin/env python3

# BEGIN_HEADER
# -----------------------------------------------------------------------------
# Gauss-Jacobi Quadrature Implementation
# Authors: Rohit Goswami <rgoswami[at]ieee.org>
# Source: GaussJacobiQuad Library
# License: MIT
# GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
# -----------------------------------------------------------------------------
# END_HEADER
import argparse
import os
import re

from jinja2 import Template

# Get the Git root directory
GITROOT = os.popen("git rev-parse --show-toplevel").read().strip()

# Same text as scripts/license_header.txt (what pre-commit inserts).
header_template = Template(
    """{{cchar}} BEGIN_HEADER
{{cchar}} -----------------------------------------------------------------------------
{{cchar}} Gauss-Jacobi Quadrature Implementation
{{cchar}} Authors: {{author}}
{{cchar}} Source: {{library_name}} Library
{{cchar}} License: MIT
{{cchar}} GitHub Repository: {{repository_url}}
{{cchar}} -----------------------------------------------------------------------------
{{cchar}} END_HEADER
""",
    keep_trailing_newline=True,
)


def _place_banner(header_text, content):
    """Put a new banner after a leading shebang, otherwise at the top."""
    if not content.startswith("#!"):
        return header_text + content
    split_at = content.find("\n")
    if split_at == -1:
        ending = "\r\n" if "\r\n" in header_text else "\n"
        return content + ending + header_text
    return content[: split_at + 1] + header_text + content[split_at + 1 :]


def add_headers(
    directories,
    file_types,
    cchar="!",
    author="Rohit Goswami <rgoswami[at]ieee.org>",
    repository_url="https://github.com/HaoZeke/GaussJacobiQuad",
    library_name="GaussJacobiQuad",
):
    # Fill in the template variables
    header_text = header_template.render(
        author=author,
        repository_url=repository_url,
        library_name=library_name,
        cchar=cchar,
    )
    # The template keeps its trailing newline; force one anyway so a write
    # cannot glue the next source line onto END_HEADER.
    header_text = header_text.rstrip("\n") + "\n"
    for directory in directories:
        # Iterate over all files in the specified directory
        for root, dirs, files in os.walk(directory):
            for filename in files:
                if any(filename.endswith(suffix) for suffix in file_types):
                    filepath = os.path.join(root, filename)

                    # newline="" keeps CRLF. Rewrite only when the text changes.
                    with open(filepath, "r", newline="") as file:
                        content = file.read()

                    line_ending = "\r\n" if "\r\n" in content else "\n"
                    file_header = header_text.replace("\n", line_ending)
                    # Replace an existing header, or prepend one. The pattern
                    # consumes the newline after END_HEADER; file_header puts
                    # exactly one back, and the following source is unchanged.
                    header_pattern = (
                        rf"{re.escape(cchar)} BEGIN_HEADER.*?"
                        rf"{re.escape(cchar)} END_HEADER\r?\n?"
                    )
                    if "BEGIN_HEADER" in content and "END_HEADER" in content:
                        updated = re.sub(
                            header_pattern,
                            lambda _match, banner=file_header: banner,
                            content,
                            flags=re.DOTALL,
                        )
                    else:
                        updated = _place_banner(file_header, content)

                    if updated != content:
                        with open(filepath, "w", newline="") as file:
                            file.write(updated)

        print("Headers added or updated in all files in" f" {directory} using {cchar}.")


if __name__ == "__main__":
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
