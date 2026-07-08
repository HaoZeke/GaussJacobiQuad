#!/usr/bin/env bash
# Build narrative docs (orgmode → HTML) and optionally Doxygen API HTML.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
ORG="$ROOT/docs/orgmode"
SITE="$ROOT/docs/site"
mkdir -p "$SITE"

if ! command -v pandoc >/dev/null 2>&1; then
  echo "error: pandoc is required for narrative docs" >&2
  exit 1
fi

echo "==> Exporting orgmode → HTML under docs/site/"
# CSS
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

nav='<nav class="site"><a href="index.html">Hub</a>
<a href="tutorials/quickstart.html">Quickstart</a>
<a href="howto/bindings.html">Bindings</a>
<a href="howto/caf.html">CAF</a>
<a href="explanation/methods.html">Methods</a>
<a href="reference/fortran-api.html">Fortran API</a>
<a href="reference/method-names.html">Methods list</a>
<span class="muted">·</span>
<a href="../html/index.html">Doxygen API</a></nav>'

export_one() {
  local src="$1"
  local rel="${src#$ORG/}"
  local dest="$SITE/${rel%.org}.html"
  mkdir -p "$(dirname "$dest")"
  # rewrite file: links to .html for same-tree org links
  local tmp
  tmp="$(mktemp)"
  # pandoc with standalone; inject nav via header-include is awkward — postprocess
  pandoc "$src" -f org -t html5 -s --css=style.css --metadata title="$(basename "${rel%.org}")" -o "$tmp"
  # Fix relative css for nested pages
  local depth="${rel//[^\/]/}"
  local css="style.css"
  if [[ -n "$depth" ]]; then
    css="$(printf '../%.0s' $(seq 1 ${#depth}))style.css"
  fi
  sed -i "s|href=\"style.css\"|href=\"$css\"|g" "$tmp"
  # prepend nav after <body>
  python3 - "$tmp" "$dest" "$nav" <<'PY'
import sys
src, dest, nav = sys.argv[1], sys.argv[2], sys.argv[3]
html = open(src, encoding="utf-8").read()
# rewrite .org file links to .html
import re
html = re.sub(r'href="([^"]+)\.org"', r'href="\1.html"', html)
html = re.sub(r"href='([^']+)\.org'", r"href='\1.html'", html)
if "<body>" in html:
    html = html.replace("<body>", "<body>\n" + nav + "\n", 1)
else:
    html = nav + html
open(dest, "w", encoding="utf-8").write(html)
PY
  rm -f "$tmp"
  echo "  $rel -> ${dest#$ROOT/}"
}

while IFS= read -r -d '' f; do
  export_one "$f"
done < <(find "$ORG" -name '*.org' -print0 | sort -z)

# copy css already written
# fix nested css paths: copy style to root only; sed already adjusted

echo "==> Narrative site: $SITE"
if command -v doxygen >/dev/null 2>&1; then
  echo "==> Doxygen API"
  (cd "$ROOT" && doxygen apidocs/Doxygen-GaussJacobiQuad.cfg)
  echo "==> Doxygen → html/"
else
  echo "==> skip doxygen (not installed); narrative site only"
fi
echo "BUILD_DOCS_OK"
