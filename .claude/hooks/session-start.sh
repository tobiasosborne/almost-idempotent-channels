#!/bin/bash
# SessionStart hook — provision the C build dependencies on ephemeral Claude Code
# on the web containers, so `make test` / `make lib` work out of the box.
#
# The fresh web container ships WITHOUT the build deps (verified: no
# flint/acb.h, no lapacke.h, no lib{flint,lapack,blas}). This installs FLINT 3.x
# (bundles arb/acb — the certified arb path) plus LAPACK/BLAS/LAPACKE (the fast
# double path). See HANDOFF.md "How to build & test" and the Makefile link line
# (-lflint -llapacke -llapack -lblas -lm).
#
# Properties (all deliberate):
#   - WEB-ONLY: skips local machines (guards on $CLAUDE_CODE_REMOTE) so it never
#     touches a developer's own system.
#   - IDEMPOTENT + FAST: no-op once the headers are present, so warm/cached
#     containers pay nothing (the web container state is cached after the hook).
#   - NON-FATAL: never blocks a session if apt or the network is unavailable —
#     `make` then simply reports the missing dep, exactly as before this hook.
#
# Out of scope (not needed for `make test`): `bd` (beads — tracking only; the
# project reads .beads/issues.jsonl directly and tolerates its absence) and
# Julia (only `make fixtures` regenerates the committed MOSEK golden master).
set -uo pipefail

# Local machine? Leave the developer's system alone — only provision web boxes.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

# Already provisioned (warm container)? Fast no-op.
if [ -f /usr/include/flint/acb_mat.h ] && [ -f /usr/include/lapacke.h ]; then
  exit 0
fi

command -v apt-get >/dev/null 2>&1 || {
  echo "session-start: no apt-get; skipping C dep install" >&2; exit 0; }

SUDO=""
if [ "$(id -u)" -ne 0 ]; then
  if command -v sudo >/dev/null 2>&1; then
    SUDO="sudo"
  else
    echo "session-start: not root and no sudo; skipping C dep install" >&2
    exit 0
  fi
fi

echo "session-start: installing FLINT + LAPACK/BLAS/LAPACKE ..." >&2
export DEBIAN_FRONTEND=noninteractive
PKGS="libflint-dev liblapacke-dev liblapack-dev libblas-dev"

# Plain install first (package lists are usually warm); refresh the index once
# and retry if that fails. All apt chatter goes to stderr to keep the hook's
# stdout (added to session context) clean.
if ! $SUDO apt-get install -y --no-install-recommends $PKGS >&2; then
  $SUDO apt-get update -qq >&2 || true
  $SUDO apt-get install -y --no-install-recommends $PKGS >&2 || {
    echo "session-start: C dep install failed (offline?); 'make' will need them" >&2
    exit 0
  }
fi

echo "session-start: build deps ready" >&2
exit 0
