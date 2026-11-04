#!/usr/bin/env bash
# Lancia la batteria delle prove del motore Julia.
# Uso: scripts/prove.sh            tutte le prove (circa 10 minuti)
#      scripts/prove.sh veloci     salta i gruppi marcati @lenta (test/lente.jl)
set -euo pipefail
cd "$(dirname "$0")/.."
if ! command -v julia >/dev/null 2>&1; then
    echo "julia non trovato: vedi la sezione 'Installare Julia' in CLAUDE.md" >&2
    exit 1
fi
case "${1:-tutte}" in
    tutte)  export MORFOGENETICA_PROVE=tutte ;;
    veloci) export MORFOGENETICA_PROVE=veloci ;;
    *) echo "uso: scripts/prove.sh [veloci]" >&2; exit 2 ;;
esac
# Manifest.toml non è versionato: lo si rigenera dalle dipendenze di Project.toml
# --threads=auto: il solutore sparso calcola la matrice in parallelo (15_sparso.jl)
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.test(; julia_args = `--threads=auto`)'
