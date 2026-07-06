#!/usr/bin/env bash
# Push this repo to https://github.com/xulinpan/scabm  (run from this folder)
set -euo pipefail
REMOTE="https://github.com/xulinpan/scabm.git"
git rev-parse --is-inside-work-tree >/dev/null 2>&1 || git init
git config user.name  >/dev/null 2>&1 || git config user.name  "Xulin Pan"
git config user.email >/dev/null 2>&1 || git config user.email "xulinpanias@gmail.com"
git add -A
git commit -m "SC-ABM submission: manuscript, code, data, NEON application" || echo "nothing to commit"
git remote get-url origin >/dev/null 2>&1 && git remote set-url origin "$REMOTE" || git remote add origin "$REMOTE"
git branch -M main
git push -u origin main
