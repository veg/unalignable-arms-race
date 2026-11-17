#!/usr/bin/env bash
# Bidirectional helper to synchronize manuscript/ (Overleaf + veg/unalignable-arms-race-paper)
# with the full study repository (veg/unalignable-arms-race).
#
# Usage:
#   ./scripts/sync_overleaf.sh          # Pull Overleaf & paper-repo changes, stage tracked manuscript files in root repo, show status
#   ./scripts/sync_overleaf.sh --push   # Also commit local manuscript edits (if any), push to Overleaf + paper repo, and push root repo

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MS="$ROOT/manuscript"
DO_PUSH=0

for arg in "$@"; do
  case "$arg" in
    --push) DO_PUSH=1 ;;
    -h|--help)
      sed -n '2,7p' "$0"
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

if [ "$DO_PUSH" -eq 1 ] && [ -n "$(git -C "$MS" status --porcelain)" ]; then
  echo "==> Committing local changes inside manuscript/..."
  git -C "$MS" add -A
  git -C "$MS" commit -m "Sync local manuscript edits"
fi

echo "==> Pulling latest commits from Overleaf (overleaf/main) and GitHub (origin/main) into manuscript/..."
git -C "$MS" pull --rebase --autostash overleaf main
git -C "$MS" pull --rebase --autostash origin main

if [ "$DO_PUSH" -eq 1 ]; then
  if [ -n "$(git -C "$MS" rev-list overleaf/main..HEAD 2>/dev/null || true)" ]; then
    echo "==> Pushing manuscript/ commits to Overleaf (overleaf/main)..."
    git -C "$MS" push overleaf main
  fi
  if [ -n "$(git -C "$MS" rev-list origin/main..HEAD 2>/dev/null || true)" ]; then
    echo "==> Pushing manuscript/ commits to GitHub (veg/unalignable-arms-race-paper)..."
    git -C "$MS" push origin main
  fi
fi

# Stage individual tracked manuscript files in the parent GitHub repo instead of
# adding the directory directly, which avoids creating a 160000 submodule gitlink
# when manuscript/.git is present.
echo "==> Mirroring Overleaf-tracked manuscript/ files into root GitHub repo index..."
git -C "$MS" ls-files -z | xargs -0 -I{} git -C "$ROOT" add -- "manuscript/{}"
git -C "$ROOT" add -u -- manuscript

if [ "$DO_PUSH" -eq 1 ]; then
  if ! git -C "$ROOT" diff --cached --quiet -- manuscript; then
    git -C "$ROOT" commit -m "Sync manuscript from Overleaf" -- manuscript
  fi
  if [ -n "$(git -C "$ROOT" rev-list @{u}..HEAD 2>/dev/null || true)" ]; then
    echo "==> Pushing synced commits to root GitHub (veg/unalignable-arms-race)..."
    git -C "$ROOT" push origin main
  fi
fi

echo
echo "==> manuscript/ (origin=veg/unalignable-arms-race-paper, overleaf=Overleaf) status:"
git -C "$MS" status --short --branch
git -C "$MS" remote -v
echo
echo "==> Root repo (veg/unalignable-arms-race) status:"
git -C "$ROOT" status --short --branch
