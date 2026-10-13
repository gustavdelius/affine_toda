#!/usr/bin/env bash
# Run a referee brief with a non-Claude model through the Antigravity CLI (agy),
# in a throwaway copy of the repository, so the referee cannot modify the worktree.
#
# Usage:  briefs/run_referee_agy.sh <brief-name> [model]
#   <brief-name>  file name in briefs/ without .md, e.g. referee_defects_36_37
#   [model]       an `agy models` id; default gemini-3.8-flash-high
#
# The report (the referee's final answer, on stdout) is written to
# briefs/<brief-name>_report.md in the real repository; the CLI's stderr goes to
# <copy>/agy.log. Paper texts are copied from $PAPERS (default ~/referee-copy/papers)
# to <copy>/papers, where the briefs expect them.
set -euo pipefail

name=${1:?usage: run_referee_agy.sh <brief-name> [model]}
model=${2:-gemini-3.8-flash-high}
papers=${PAPERS:-$HOME/referee-copy/papers}

repo=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
brief="$repo/briefs/$name.md"
[ -f "$brief" ] || { echo "no brief $brief" >&2; exit 1; }

copy="$HOME/referee-copy/$name"
case "$copy" in "$HOME"/referee-copy/?*) ;; *) echo "refusing copy dir $copy" >&2; exit 1;; esac
rm -rf "$copy"
mkdir -p "$copy"
rsync -a --exclude .git --exclude .claude --exclude _site --exclude .quarto "$repo"/ "$copy"/
if [ -d "$papers" ]; then cp -r "$papers" "$copy"/papers; else echo "warning: no papers at $papers" >&2; fi

cd "$copy"
echo "running $name with $model in $copy" >&2
agy -p "$(cat "$brief")" --model "$model" --dangerously-skip-permissions \
    > "$repo/briefs/${name}_report.md" 2> "$copy/agy.log"
echo "report: $repo/briefs/${name}_report.md" >&2
