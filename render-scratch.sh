#!/usr/bin/env bash
# Render the book in a scratch copy of this checkout and report problems.
#
# Do not run `quarto render` in a checkout where `quarto preview` may be
# running: the two race (missing *_files directories, a failed rename of
# index.html). `--output-dir` does not help, because Quarto still writes each
# page next to its source and only moves it at the end. This script copies the
# working tree, uncommitted changes included, to a scratch directory, renders
# there, and lists unresolved cross-references and render warnings.
#
# Usage: ./render-scratch.sh [--pdf] [--project DIR] [scratch-dir]
#   default        HTML only; scratch dir ${TMPDIR:-/tmp}/affine-toda-render-<checkout name>
#   --pdf          HTML and PDF (about 7 minutes; needs LaTeX)
#   --project DIR  render the Quarto project nested in DIR, such as
#                  lecture-notes/defects, instead of the book
# The site lands in <scratch-dir>/book/[DIR/]_site and the log in <scratch-dir>/render.log.
# To look at it: python3 -m http.server -d <scratch-dir>/book/[DIR/]_site 8000
# Exit status: non-zero if the render fails or a cross-reference is unresolved.
set -euo pipefail
cd "$(dirname "$0")"

formats=(--to html)
project=.
while [ $# -gt 0 ]; do
  case "$1" in
    --pdf) formats=(); shift ;;
    --project) project=${2:?--project needs a directory}; shift 2 ;;
    *) break ;;
  esac
done
scratch=${1:-${TMPDIR:-/tmp}/affine-toda-render-$(basename "$PWD")}
book=$scratch/book
log=$scratch/render.log
mkdir -p "$book"

rsync -a --delete \
  --exclude .git --exclude .claude --exclude _site --exclude .quarto \
  --exclude referee_scratch \
  ./ "$book/"

echo "Rendering in $book (log: $log) ..."
start=$SECONDS
status=0
(cd "$book/$project" && quarto render "${formats[@]}") >"$log" 2>&1 || status=$?
echo "quarto render exited with $status after $((SECONDS - start)) s."

# Quarto 1.9 writes an unresolved reference as <span class="quarto-unresolved-ref">?sec-x</span>;
# older versions wrote ?@sec-x.
site=$book/$project/_site
unresolved=$(grep -rhoE 'quarto-unresolved-ref">[^<]*|\?@[a-z]+-[A-Za-z0-9_-]+' "$site" --include='*.html' 2>/dev/null \
  | sed 's/^quarto-unresolved-ref">//' | sort | uniq -c || true)
if [ -n "$unresolved" ]; then
  echo "Unresolved cross-references in the HTML (count, target):"
  echo "$unresolved" | sed 's/^/  /'
  if [ "$status" -eq 0 ]; then status=1; fi
else
  echo "No unresolved cross-references in the HTML."
fi

warnings=$(sed 's/\x1b\[[0-9;]*m//g' "$log" | grep -iE 'warn|unable to resolve|error' | sort | uniq -c | sort -rn || true)
if [ -n "$warnings" ]; then
  echo "Warnings and errors in the log (count, line):"
  echo "$warnings" | head -20 | sed 's/^/  /'
fi

if [ ${#formats[@]} -eq 0 ]; then
  for pdf in "$site"/*.pdf; do
    [ -f "$pdf" ] || continue
    pages=$(pdfinfo "$pdf" 2>/dev/null | awk '/^Pages:/ {print $2}')
    marks=$(pdftotext "$pdf" - 2>/dev/null | grep -c '??' || true)
    echo "PDF ($(basename "$pdf")): ${pages:-?} pages; $marks lines contain '??' (unresolved LaTeX references, or genuine text)."
  done
fi

exit "$status"
