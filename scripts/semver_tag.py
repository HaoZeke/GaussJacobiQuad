#!/usr/bin/env python3
"""SemVer 2.0.0 checks for the release workflow.

The expression is the numbered-group pattern published at semver.org.
Tags take a leading v. NEWS headings use the same pattern with no v.
python3 re compiles it (re.ASCII, so \\d is [0-9]). POSIX ERE has no \\d
and no non-capturing groups, so release.yml calls this instead of grep -E.
"""

import re
import sys
from pathlib import Path

# https://semver.org/#is-there-a-suggested-regular-expression-regex-to-check-a-semver-string
OFFICIAL_SEMVER = (
    r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)"
    r"(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)"
    r"(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?"
    r"(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"
)

TAG_RE = re.compile("^v" + OFFICIAL_SEMVER[1:], re.ASCII)
NEWS_BOUNDARY_RE = re.compile(
    r"^#? ?" + OFFICIAL_SEMVER[1:-1] + r" ",
    re.ASCII,
)


def is_release_tag(tag: str) -> bool:
    return TAG_RE.fullmatch(tag) is not None


def extract_notes(text: str, version: str) -> str:
    """Return the NEWS section for version.

    The section runs through the line before the next SemVer heading.
    """
    start = re.compile(r"^#? ?" + re.escape(version) + r" ")
    grabbed = []
    started = False
    for line in text.splitlines():
        if not started and start.match(line):
            started = True
            grabbed.append(line)
            continue
        if started and NEWS_BOUNDARY_RE.match(line):
            break
        if started:
            grabbed.append(line)
    if not grabbed:
        return ""
    return "\n".join(grabbed) + "\n"


def main(argv: list[str]) -> int:
    if len(argv) >= 2 and argv[1] == "--notes":
        if len(argv) != 4:
            print("usage: semver_tag.py --notes VERSION NEWS.md", file=sys.stderr)
            return 2
        text = Path(argv[3]).read_text(encoding="utf-8")
        sys.stdout.write(extract_notes(text, argv[2]))
        return 0
    if len(argv) != 2:
        print("usage: semver_tag.py TAG", file=sys.stderr)
        return 2
    return 0 if is_release_tag(argv[1]) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
