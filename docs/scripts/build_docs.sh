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
# GJP_DOCS_SKIP_DOXYGEN=1 — skip the optional Doxygen pass
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export GJP_DOCS_API_HREF="${GJP_DOCS_API_HREF:-../../html/index.html}"

python3 "$ROOT/docs/scripts/export_narrative.py"

if [[ "${GJP_DOCS_SKIP_DOXYGEN:-0}" == "1" ]]; then
  echo "==> skip doxygen (GJP_DOCS_SKIP_DOXYGEN=1)"
elif command -v doxygen >/dev/null 2>&1; then
  echo "==> Doxygen API"
  (cd "$ROOT" && doxygen apidocs/Doxygen-GaussJacobiQuad.cfg)
  echo "==> Doxygen → html/"
else
  echo "==> skip doxygen (not installed); narrative site only"
fi
echo "BUILD_DOCS_OK"
