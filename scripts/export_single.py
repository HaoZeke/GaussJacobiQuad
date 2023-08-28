#!/usr/bin/env python3
# Cannot use ! in statements with this if strip-comments is true

import argparse
import os
import re
import subprocess
from pathlib import Path

import add_headers
import fortdepend


def get_git_root():
    try:
        git_root = (
            subprocess.check_output(["git", "rev-parse", "--show-toplevel"])
            .strip()
            .decode("utf-8")
        )
        return Path(git_root)
    except subprocess.CalledProcessError:
        print("This directory is not part of a Git repository.")
        return None


def process_files(fname, strip_comments=True):
    with open(fname, "r") as f:
        content = f.read()
    # Remove everything between ! BEGIN_HEADER and ! END_HEADER
    content = re.sub(r"! BEGIN_HEADER.*?! END_HEADER\n", "", content, flags=re.DOTALL)
    # Optionally, remove lines starting with !
    if strip_comments:
        content = re.sub(r"(!.*?$)|(\s*!.*?$)", "", content, flags=re.MULTILINE)
    return content


def export_single(modulename, strip_comments=True, outname=None):
    if outname is None:
        outname = f"{modulename}_single.f90"

    prj_srcs = get_git_root() / "src"
    os.chdir(prj_srcs)  # Needed for fortdepend

    gaussjacobiquad = fortdepend.FortranProject()
    concat_files = gaussjacobiquad.get_all_used_files(modulename)
    print(f"Concatenating in reverse order {concat_files}")

    exp_dat = []

    for idx, fname in enumerate(reversed(concat_files)):
        processed_content = process_files(prj_srcs / fname, strip_comments)
        exp_dat.append(processed_content)

    with Path(get_git_root() / "dist" / outname).open("w") as f_out:
        f_out.write("".join(exp_dat))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export single file algorithms.")
    parser.add_argument(
        "--modulename", type=str, required=True, help="Module to export"
    )
    parser.add_argument(
        "--outname", type=str, required=False, help="Defaults to modulename_single.f90"
    )
    parser.add_argument(
        "--strip-comments",
        dest="strip_comments",
        action="store_true",
        help="Remove comments",
    )
    parser.add_argument(
        "--keep-comments",
        dest="strip_comments",
        action="store_false",
        help="Keep comments",
    )
    parser.set_defaults(strip_comments=True)

    args = parser.parse_args()
    export_single(args.modulename, args.strip_comments, args.outname)
    add_headers.add_headers([get_git_root() / "dist"], ["f90"])
