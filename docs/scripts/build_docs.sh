#!/usr/bin/env bash
# Build narrative docs (orgmode → HTML) and optionally Doxygen API HTML.
#
# Site-root-relative navigation is rewritten per page depth so nested HTML
# (howto/, tutorials/, …) resolves hub/API links correctly.
#
# GJP_DOCS_API_HREF — path from *site root* to Doxygen index.html
#   Local default (site under docs/site/, Doxygen under repo html/):
#     ../../html/index.html
#   CI GitHub Pages (site under html/narrative/, Doxygen under html/):
#     ../index.html
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
ORG="$ROOT/docs/orgmode"
SITE="$ROOT/docs/site"
API_HREF="${GJP_DOCS_API_HREF:-../../html/index.html}"
mkdir -p "$SITE"

if ! command -v pandoc >/dev/null 2>&1; then
  echo "error: pandoc is required for narrative docs" >&2
  exit 1
fi

echo "==> Exporting orgmode → HTML under docs/site/ (API href from site root: $API_HREF)"

cat > "$SITE/style.css" << 'CSS'
:root { --fg: #1a1a1a; --muted: #555; --link: #0b5fff; --bg: #fafafa; --code-bg: #f0f0f0; }
@media (prefers-color-scheme: dark) {
  :root { --fg: #e8e8e8; --muted: #aaa; --link: #6db3ff; --bg: #121212; --code-bg: #1e1e1e; }
}
body { font-family: system-ui, sans-serif; line-height: 1.55; max-width: 52rem; margin: 1.5rem auto; padding: 0 1rem; color: var(--fg); background: var(--bg); }
a { color: var(--link); }
nav.site { font-size: 0.95rem; margin-bottom: 1.5rem; padding-bottom: 0.75rem; border-bottom: 1px solid var(--muted); }
nav.site a { margin-right: 0.75rem; }
pre, code { font-family: ui-monospace, monospace; background: var(--code-bg); }
pre { padding: 0.75rem 1rem; overflow-x: auto; border-radius: 4px; }
table { border-collapse: collapse; width: 100%; margin: 1rem 0; }
th, td { border: 1px solid var(--muted); padding: 0.35rem 0.6rem; text-align: left; }
h1, h2, h3 { line-height: 1.25; }
.muted { color: var(--muted); font-size: 0.9rem; }
CSS

export_one() {
  local src="$1"
  local rel="${src#$ORG/}"
  local dest="$SITE/${rel%.org}.html"
  mkdir -p "$(dirname "$dest")"
  local tmp
  tmp="$(mktemp)"
  pandoc "$src" -f org -t html5 -s --css=style.css \
    --metadata title="$(basename "${rel%.org}")" -o "$tmp"

  python3 - "$tmp" "$dest" "$rel" "$API_HREF" <<'PY'
import re
import sys
from pathlib import PurePosixPath

src, dest, rel, api_href = sys.argv[1:5]
# depth = number of directory components in rel (e.g. howto/bindings.org → 1)
parts = PurePosixPath(rel).parts
depth = max(0, len(parts) - 1)
prefix = "../" * depth

# Site-root-relative destinations → page-relative
nav_targets = [
    ("Hub", "index.html"),
    ("Quickstart", "tutorials/quickstart.html"),
    ("Bindings", "howto/bindings.html"),
    ("CAF", "howto/caf.html"),
    ("Methods", "explanation/methods.html"),
    ("Fortran API", "reference/fortran-api.html"),
    ("Methods list", "reference/method-names.html"),
]
links = []
for label, target in nav_targets:
    links.append(f'<a href="{prefix}{target}">{label}</a>')
api_link = f'{prefix}{api_href}'
nav = (
    '<nav class="site">'
    + "\n".join(links)
    + f'\n<span class="muted">·</span>\n'
    + f'<a href="{api_link}">Doxygen API</a></nav>'
)

css = f"{prefix}style.css"
html = open(src, encoding="utf-8").read()
html = re.sub(r'href="style\.css"', f'href="{css}"', html)
html = re.sub(r'href="([^"]+)\.org"', r'href="\1.html"', html)
html = re.sub(r"href='([^']+)\.org'", r"href='\1.html'", html)
if "<body>" in html:
    html = html.replace("<body>", "<body>\n" + nav + "\n", 1)
else:
    html = nav + html
open(dest, "w", encoding="utf-8").write(html)
print(f"  depth={depth} api={api_link} -> {dest}")
PY
  rm -f "$tmp"
  echo "  $rel -> ${dest#$ROOT/}"
}

while IFS= read -r -d '' f; do
  export_one "$f"
done < <(find "$ORG" -name '*.org' -print0 | sort -z)

echo "==> Narrative site: $SITE"
if command -v doxygen >/dev/null 2>&1; then
  echo "==> Doxygen API"
  (cd "$ROOT" && doxygen apidocs/Doxygen-GaussJacobiQuad.cfg)
  echo "==> Doxygen → html/"
else
  echo "==> skip doxygen (not installed); narrative site only"
fi
echo "BUILD_DOCS_OK"
