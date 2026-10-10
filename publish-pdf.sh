#!/usr/bin/env bash
# Build the book PDF locally and put it on the site.
#
# The GitHub deploy (.github/workflows/publish.yml) renders HTML only and
# copies affine-toda.pdf from the pdf-latest release into the site. This
# script renders the PDF, replaces it in that release (creating the release
# the first time) and redeploys the site. The PDF is not rendered on GitHub
# because Quarto hangs there drawing the mermaid figure of chapter 1.
#
# Needs quarto with a LaTeX installation, and gh logged in to GitHub.
set -euo pipefail
cd "$(dirname "$0")"

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "Warning: uncommitted changes will be in the PDF." >&2
fi
git fetch -q origin main
if ! git merge-base --is-ancestor HEAD origin/main; then
  echo "Warning: HEAD is not on origin/main, so the PDF will not match the site." >&2
fi

quarto render --to pdf
pdf=_site/affine-toda.pdf

notes="Built locally from $(git rev-parse --short HEAD) on $(date -u +%Y-%m-%d)."
if gh release view pdf-latest >/dev/null 2>&1; then
  gh release upload pdf-latest "$pdf" --clobber
  gh release edit pdf-latest --notes "$notes"
else
  gh release create pdf-latest "$pdf" --target main --title "Latest PDF" --notes "$notes"
fi
gh workflow run publish.yml --ref main
echo "Uploaded $pdf; the site redeploys in about two minutes."
