#!/usr/bin/env python3
"""Export docs/orgmode/*.org → docs/site/*.html with site-root-relative nav.

Invoked by docs/scripts/build_docs.sh. Avoids bash while-read + heredoc stdin
races that stalled GHA after the first org file.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[2]
ORG = ROOT / "docs" / "orgmode"
SITE = ROOT / "docs" / "site"
API_HREF = os.environ.get("GJP_DOCS_API_HREF", "../../html/index.html")

NAV_TARGETS = [
    ("Hub", "index.html"),
    ("Quickstart", "tutorials/quickstart.html"),
    ("Bindings", "howto/bindings.html"),
    ("CAF", "howto/caf.html"),
    ("Methods", "explanation/methods.html"),
    ("Fortran API", "reference/fortran-api.html"),
    ("Methods list", "reference/method-names.html"),
]

CSS = (
    ":root { --fg: #1a1a1a; --muted: #555; --link: #0b5fff;\n"
    "  --bg: #fafafa; --code-bg: #f0f0f0; }\n"
    "@media (prefers-color-scheme: dark) {\n"
    "  :root { --fg: #e8e8e8; --muted: #aaa; --link: #6db3ff;\n"
    "    --bg: #121212; --code-bg: #1e1e1e; }\n"
    "}\n"
    "body { font-family: system-ui, sans-serif; line-height: 1.55;\n"
    "  max-width: 52rem; margin: 1.5rem auto; padding: 0 1rem;\n"
    "  color: var(--fg); background: var(--bg); }\n"
    "a { color: var(--link); }\n"
    "nav.site { font-size: 0.95rem; margin-bottom: 1.5rem;\n"
    "  padding-bottom: 0.75rem; border-bottom: 1px solid var(--muted); }\n"
    "nav.site a { margin-right: 0.75rem; }\n"
    "pre, code { font-family: ui-monospace, monospace;\n"
    "  background: var(--code-bg); }\n"
    "pre { padding: 0.75rem 1rem; overflow-x: auto; border-radius: 4px; }\n"
    "table { border-collapse: collapse; width: 100%; margin: 1rem 0; }\n"
    "th, td { border: 1px solid var(--muted); padding: 0.35rem 0.6rem;\n"
    "  text-align: left; }\n"
    "h1, h2, h3 { line-height: 1.25; }\n"
    ".muted { color: var(--muted); font-size: 0.9rem; }\n"
)



def export_one(src: Path) -> None:
    rel = src.relative_to(ORG).as_posix()
    dest = SITE / PurePosixPath(rel).with_suffix(".html")
    dest.parent.mkdir(parents=True, exist_ok=True)
    title = PurePosixPath(rel).stem
    print(f"  pandoc {rel} ...", flush=True)
    with tempfile.NamedTemporaryFile(
        mode="w+", encoding="utf-8", suffix=".html", delete=False
    ) as tmp:
        tmp_path = Path(tmp.name)
    try:
        subprocess.run(
            [
                "pandoc",
                str(src),
                "-f",
                "org",
                "-t",
                "html5",
                "-s",
                "--css=style.css",
                f"--metadata=title={title}",
                "-o",
                str(tmp_path),
            ],
            check=True,
            stdin=subprocess.DEVNULL,
            timeout=60,
        )
        parts = PurePosixPath(rel).parts
        depth = max(0, len(parts) - 1)
        prefix = "../" * depth
        links = [f'<a href="{prefix}{t}">{lab}</a>' for lab, t in NAV_TARGETS]
        api_link = f"{prefix}{API_HREF}"
        nav = (
            '<nav class="site">'
            + "\n".join(links)
            + '\n<span class="muted">·</span>\n'
            + f'<a href="{api_link}">Doxygen API</a></nav>'
        )
        html = tmp_path.read_text(encoding="utf-8")
        css = f"{prefix}style.css"
        html = re.sub(r'href="style\.css"', f'href="{css}"', html)
        html = re.sub(r'href="([^"]+)\.org"', r'href="\1.html"', html)
        html = re.sub(r"href='([^']+)\.org'", r"href='\1.html'", html)
        if "<body>" in html:
            html = html.replace("<body>", "<body>\n" + nav + "\n", 1)
        else:
            html = nav + html
        dest.write_text(html, encoding="utf-8")
        print(f"  depth={depth} api={api_link} -> {dest}", flush=True)
    finally:
        tmp_path.unlink(missing_ok=True)


def main() -> int:
    if not shutil_which("pandoc"):
        print("error: pandoc is required for narrative docs", file=sys.stderr)
        return 1
    SITE.mkdir(parents=True, exist_ok=True)
    (SITE / "style.css").write_text(CSS, encoding="utf-8")
    org_files = sorted(ORG.rglob("*.org"))
    print(
        f"==> Exporting {len(org_files)} org files → {SITE} "
        f"(API href from site root: {API_HREF})",
        flush=True,
    )
    for src in org_files:
        export_one(src)
    print(f"==> Narrative site: {SITE}", flush=True)
    return 0


def shutil_which(cmd: str) -> str | None:
    from shutil import which

    return which(cmd)


if __name__ == "__main__":
    sys.exit(main())
