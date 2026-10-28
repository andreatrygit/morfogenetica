"""
    Morfogenetica

Motore morfogenetico: cellule con una macchina a stati che fanno crescere forme da sole regole
locali. La specifica è in `docs/specifica-motore-julia.md`.

Il codice vive nelle cartelle descritte nel README, non in `src/`: questo file le raccoglie.
- `motore/`     geometria, energia, topologia, chimica, passata
- `primitive/`  una regola per file, pura, con il suo predicato di fine
"""
module Morfogenetica

const RADICE = normpath(joinpath(@__DIR__, ".."))

# Include in ordine alfabetico ogni file .jl di una cartella del repository.
function _includi_cartella(cartella)
    percorso = joinpath(RADICE, cartella)
    isdir(percorso) || return
    for f in sort(readdir(percorso))
        endswith(f, ".jl") && include(joinpath(percorso, f))
    end
end

_includi_cartella("motore")
_includi_cartella("primitive")

end # module
