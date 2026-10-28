#!/usr/bin/env bash
# Lancia la batteria delle prove del motore Julia.
# Uso: scripts/prove.sh
set -euo pipefail
cd "$(dirname "$0")/.."
if ! command -v julia >/dev/null 2>&1; then
    echo "julia non trovato: vedi la sezione 'Installare Julia' in CLAUDE.md" >&2
    exit 1
fi
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
