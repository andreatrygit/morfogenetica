"""
    Morfogenetica

Motore morfogenetico: cellule con una macchina a stati che fanno crescere forme da sole regole
locali. La specifica è in `docs/specifica-motore-julia.md`, il piano in `docs/piano.md`.

Il codice vive nelle cartelle descritte nel README, non in `src/`: questo file le raccoglie.
Ogni file nuovo va aggiunto qui sotto, in ordine.
"""
module Morfogenetica

const RADICE = normpath(joinpath(@__DIR__, ".."))

# motore/: geometria, energia, dinamica, misure, contatto, integratore, guscio 3D, scambi di vicine, estrusione, scutoidi
include(joinpath(RADICE, "motore", "01_geometria.jl"))
include(joinpath(RADICE, "motore", "02_energia.jl"))
include(joinpath(RADICE, "motore", "03_rilassamento.jl"))
include(joinpath(RADICE, "motore", "04_misure.jl"))
include(joinpath(RADICE, "motore", "05_contatto.jl"))
include(joinpath(RADICE, "motore", "06_integratore.jl"))
include(joinpath(RADICE, "motore", "07_guscio.jl"))
include(joinpath(RADICE, "motore", "08_scambi.jl"))
include(joinpath(RADICE, "motore", "09_estrusione.jl"))
include(joinpath(RADICE, "motore", "10_scutoidi.jl"))
include(joinpath(RADICE, "motore", "11_operazioni.jl"))

end # module
