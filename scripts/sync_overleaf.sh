#!/usr/bin/env bash
# Bidirectional helper to synchronize manuscript/ (Overleaf) with GitHub (veg/unalignable-arms-race).
#
# Usage:
#   ./scripts/sync_overleaf.sh          # Pull Overleaf changes, stage tracked manuscript files in parent repo, show status
#   ./scripts/sync_overleaf.sh --push   # Also commit local manuscript edits (if any), push to Overleaf, and push synced commits to GitHub

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MS="$ROOT/manuscript"
DO_PUSH=0

for arg in "$@"; do
  case "$arg" in
    --push) DO_PUSH=1 ;;
    -h|--help)
      sed -n '2,6p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown option: $arg" >&2
      exit 1
      ;;
  esac
done

if [ ! -d "$MS/.git" ]; then
  echo "Error: $MS/.git does not exist. Wire the Overleaf clone first." >&2
  exit 1
fi

echo "==> Pulling latest commits from Overleaf (origin/main) into manuscript/..."
if [ "$DO_PUSH" -eq 1 ] && [ -n "$(git -C "$MS" status --porcelain)" ]; then
  git -C "$MS" add -A
  git -C "$MS" commit -m "Sync local manuscript edits"
fi
git -C "$MS" pull --rebase --autostash origin main

if [ "$DO_PUSH" -eq 1 ]; then
  if [ -n "$(git -C "$MS" rev-list @{u}..HEAD 2>/dev/null || true)" ]; then
    echo "==> Pushing manuscript/ commits to Overleaf..."
    git -C "$MS" push origin main
  fi
fi

# Stage individual tracked manuscript files in the parent GitHub repo instead of
# adding the directory directly, which avoids creating a 160000 submodule gitlink
# when manuscript/.git is present.
echo "==> Mirroring Overleaf-tracked manuscript/ files into GitHub repo index..."
git -C "$MS" ls-files -z | xargs -0 -I{} git -C "$ROOT" add -- "manuscript/{}"
git -C "$ROOT" add -u -- manuscript

if [ "$DO_PUSH" -eq 1 ]; then
  if ! git -C "$ROOT" diff --cached --quiet -- manuscript; then
    git -C "$ROOT" commit -m "Sync manuscript from Overleaf" -- manuscript
  fi
  if [ -n "$(git -C "$ROOT" rev-list @{u}..HEAD 2>/dev/null || true)" ]; then
    echo "==> Pushing synced commits to GitHub (origin/main)..."
    git -C "$ROOT" push origin main
  fi
fi

echo
echo "==> manuscript/ (Overleaf) status:"
git -C "$MS" status --short --branch
echo
echo "==> Root repo (GitHub) status:"
git -C "$ROOT" status --short --branch
