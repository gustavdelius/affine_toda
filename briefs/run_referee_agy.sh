#!/usr/bin/env bash
# Run a referee brief with Codex (`codex exec`) or Gemini (`agy -p`) in a
# throwaway copy of the repository, so the referee cannot modify the worktree.
#
# Usage:  briefs/run_referee_agy.sh <brief-name> [agent-or-model]
#   <brief-name>      file name in briefs/ without .md, e.g. referee_classical_09_11
#                     (or referee_classical_09_11_items for phase 2, which reuses
#                     the phase-1 copy and resumes a Codex session if present)
#   [agent-or-model]  codex (default) | gemini | both | <codex-model> | <agy-model>
#                     - codex, gpt-*, o3*, o4*: runs `codex exec` (reasoning effort high)
#                     - gemini, agy, gemini-*, claude-*: runs `agy -p`
#                       (default model when `gemini` or `agy` is given: gemini-3.8-flash-high)
#                     - both: runs codex and gemini in parallel (for blind phase-1 passes)
#
# Paper texts are copied from $PAPERS (default ~/referee-copy/papers) to
# <copy>/papers, where the briefs expect them. When the referee finishes, the
# report file it wrote in <copy>/briefs/ (or its final stdout message if it did
# not write the file directly) is copied back to briefs/ in the real repository,
# and any scratch scripts under <copy>/referee_scratch/ are copied back to
# referee_scratch/.
set -euo pipefail

name=${1:?usage: run_referee_agy.sh <brief-name> [codex|gemini|both|<model>]}
target=${2:-codex}
papers=${PAPERS:-$HOME/referee-copy/papers}

if [ "$target" = "both" ]; then
    "$0" "$name" codex &
    pid_codex=$!
    "$0" "$name" gemini &
    pid_gemini=$!
    rc=0
    wait "$pid_codex" || rc=$?
    wait "$pid_gemini" || rc=$?
    exit "$rc"
fi

case "$target" in
    codex)
        backend=codex
        model=""
        ;;
    gpt-*|o3*|o4*)
        backend=codex
        model="$target"
        ;;
    gemini|agy)
        backend=agy
        model="gemini-3.8-flash-high"
        ;;
    *)
        backend=agy
        model="$target"
        ;;
esac

repo=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
brief="$repo/briefs/$name.md"
[ -f "$brief" ] || { echo "no brief $brief" >&2; exit 1; }

is_phase2=0
base="$name"
if [[ "$name" == *_items ]]; then
    is_phase2=1
    base="${name%_items}"
fi

copy="$HOME/referee-copy/${base}-${backend}"
case "$copy" in "$HOME"/referee-copy/?*) ;; *) echo "refusing copy dir $copy" >&2; exit 1;; esac

if [ "$is_phase2" -eq 1 ] && [ -d "$copy" ]; then
    # Reuse the phase-1 copy and refresh briefs/ so the committed phase-1 report
    # and the phase-2 *_items.md brief are available inside the copy.
    mkdir -p "$copy/briefs"
    cp -a "$repo/briefs/${base}"*.md "$copy/briefs/" 2>/dev/null || true
    cp -a "$brief" "$copy/briefs/"
else
    rm -rf "$copy"
    mkdir -p "$copy"
    rsync -a --exclude .git --exclude .claude --exclude _site --exclude .quarto "$repo"/ "$copy"/
    if [ -d "$papers" ]; then cp -r "$papers" "$copy"/papers; else echo "warning: no papers at $papers" >&2; fi
fi

if [ "$is_phase2" -eq 1 ]; then
    expected_report="${base}_report.md"
elif grep -q "${base}_phase1_report\.md" "$brief"; then
    expected_report="${base}_phase1_report.md"
else
    expected_report="${base}_report.md"
fi

cd "$copy"
touch "$copy/.run_start"
rm -f "$copy/last_message.md"

if [ "$backend" = "codex" ]; then
    if command -v codex >/dev/null 2>&1; then
        codex_bin=codex
    elif [ -x /usr/lib/chatgpt/resources/codex ]; then
        codex_bin=/usr/lib/chatgpt/resources/codex
    else
        echo "codex CLI not found in PATH or /usr/lib/chatgpt/resources/codex" >&2
        exit 1
    fi
    echo "running $name with codex ${model:+($model) }in $copy" >&2
    if [ "$is_phase2" -eq 1 ] && [ -f "$copy/.codex_phase1_done" ]; then
        "$codex_bin" exec resume --last \
            --skip-git-repo-check --dangerously-bypass-approvals-and-sandbox \
            ${model:+-m "$model"} -c 'model_reasoning_effort="high"' \
            -o "$copy/last_message.md" - < "$brief" \
            > "$copy/codex.out" 2> "$copy/codex.log"
    else
        "$codex_bin" exec --cd "$copy" \
            --skip-git-repo-check --dangerously-bypass-approvals-and-sandbox \
            ${model:+-m "$model"} -c 'model_reasoning_effort="high"' \
            -o "$copy/last_message.md" - < "$brief" \
            > "$copy/codex.out" 2> "$copy/codex.log"
        touch "$copy/.codex_phase1_done"
    fi
else
    echo "running $name with agy ($model) in $copy" >&2
    agy -p "$(cat "$brief")" --model "$model" --dangerously-skip-permissions \
        > "$copy/last_message.md" 2> "$copy/agy.log"
fi

if [ -s "$copy/briefs/$expected_report" ] && [ "$copy/briefs/$expected_report" -nt "$copy/.run_start" ]; then
    src_report="$copy/briefs/$expected_report"
elif [ -s "$copy/last_message.md" ]; then
    src_report="$copy/last_message.md"
else
    echo "error: no report produced in $copy/briefs/$expected_report or $copy/last_message.md" >&2
    exit 1
fi

dst_report="$repo/briefs/$expected_report"
if [ -s "$dst_report" ]; then
    # Avoid overwriting a non-empty report already written by another model.
    if [[ "$expected_report" == *_phase1_report.md ]]; then
        dst_report="$repo/briefs/${base}_${backend}_phase1_report.md"
    else
        dst_report="$repo/briefs/${base}_${backend}_report.md"
    fi
fi

cp "$src_report" "$dst_report"
if [ -d "$copy/referee_scratch" ]; then
    mkdir -p "$repo/referee_scratch"
    cp -rn "$copy/referee_scratch/"* "$repo/referee_scratch/" 2>/dev/null || true
fi
echo "report: $dst_report" >&2

