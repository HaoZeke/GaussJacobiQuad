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
import sys

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


_MODE_LINE = re.compile(r"^![ \t]*-\*-[ \t]*f90[ \t]*-\*-[ \t]*$")


def _is_mode_line(line):
    return _MODE_LINE.match(line.rstrip("\r\n")) is not None


def _place_banner(header_text, content):
    """Put a new banner after a leading shebang or f2py mode line.

    Blank lines before the mode line stay above the banner, including when
    the mode line itself has no trailing newline.
    """
    if content == "":
        return header_text
    lines = content.splitlines(keepends=True)
    idx = 0
    if lines[0].startswith("#!"):
        if "\n" not in lines[0]:
            ending = "\r\n" if "\r\n" in header_text else "\n"
            lines[0] = lines[0] + ending
        idx = 1
    mode_at = idx
    while mode_at < len(lines) and lines[mode_at].strip() == "":
        mode_at += 1
    if mode_at < len(lines) and _is_mode_line(lines[mode_at]):
        if "\n" not in lines[mode_at]:
            ending = "\r\n" if "\r\n" in header_text else "\n"
            lines[mode_at] = lines[mode_at] + ending
        idx = mode_at + 1
    if idx == 0:
        return header_text + content
    return "".join(lines[:idx]) + header_text + "".join(lines[idx:])


def _comment_line(line, cchar):
    body = line.lstrip(" \t")
    if cchar == "#" and body.startswith("#!"):
        return False
    return body.startswith(cchar)


def _marker_line(line, cchar, marker):
    return line.strip() == "%s %s" % (cchar, marker)


def _apply_header(content, file_header, cchar, filepath):
    """Replace one banner in the leading comment block, or insert one.

    BEGIN_HEADER and END_HEADER must both be comment lines in that block.
    An unclosed leading banner is left unchanged and reported.
    """
    if content == "":
        return file_header
    lines = content.splitlines(keepends=True)
    idx = 0
    if lines[0].startswith("#!"):
        idx = 1
    block_at = idx
    while idx < len(lines) and (
        lines[idx].strip() == "" or _comment_line(lines[idx], cchar)
    ):
        idx += 1
    begin = None
    for i in range(block_at, idx):
        if _marker_line(lines[i], cchar, "BEGIN_HEADER"):
            begin = i
            break
    if begin is None:
        return _place_banner(file_header, content)
    end = None
    for j in range(begin + 1, idx):
        if _marker_line(lines[j], cchar, "END_HEADER"):
            end = j
            break
    if end is None:
        print(
            "warning: %s has BEGIN_HEADER but no END_HEADER "
            "in the leading comment block; left unchanged" % filepath,
            file=sys.stderr,
        )
        return content
    return "".join(lines[:begin]) + file_header + "".join(lines[end + 1 :])


def _line_ending(content):
    """Use the first line break, so a later CRLF does not flip the banner."""
    nl = content.find("\n")
    if nl > 0 and content[nl - 1] == "\r":
        return "\r\n"
    return "\n"


def _extension(filename):
    dot = filename.rfind(".")
    if dot < 0:
        return ""
    return filename[dot + 1 :]


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
            dirs[:] = [d for d in dirs if d != "__pycache__"]
            for filename in files:
                # `.pyc` ends with `c`; match the extension, and skip bytecode.
                if filename.endswith(".pyc"):
                    continue
                if not any(_extension(filename) == suffix for suffix in file_types):
                    continue
                filepath = os.path.join(root, filename)

                # newline="" keeps CRLF. Rewrite only when the text changes.
                with open(filepath, "r", newline="", encoding="utf-8") as file:
                    content = file.read()

                line_ending = _line_ending(content)
                file_header = header_text.replace("\n", line_ending)
                updated = _apply_header(content, file_header, cchar, filepath)

                if updated != content:
                    with open(filepath, "w", newline="", encoding="utf-8") as file:
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
