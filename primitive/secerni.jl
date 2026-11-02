# secerni: emettere un segnale (repertorio, "emettere e leggere un segnale"). La cellula marcata
# secerne la sostanza a tasso costante alla propria posizione; le altre no.
#
# Regola, sul canale "dentro di sé": se stato[marca] == 1 secerno `tasso` (quantità/min).
# Predicato di fine: nessuno.
#
# Meccanismo: le cellule sorgente esprimono e secernono una proteina segnale, che diffonde e si
#   degrada nel tessuto (Dpp nel disco dell'ala di Drosophila). Chi è marcata lo decide una fase
#   precedente, o la condizione iniziale.
# Fonte: Kicheva et al., "Kinetics of morphogen gradient formation", Science 315, 521–525 (2007),
#   doi:10.1126/science.1135774.

export secerni

function _regola_secerni(p::Percezione, par)
    q = get(p.stato, par.marca, 0.0) == 1 ? par.tasso : 0.0
    return Uscita(; secrezioni = Dict{Symbol,Float64}(par.sostanza => q))
end

_fine_secerni(p::Percezione, par) = false

"""
    secerni(; sostanza, tasso, marca = :sorgente)

La fase in cui la cellula con stato[marca] == 1 secerne `sostanza` a `tasso` (quantità/min).
"""
secerni(; sostanza, tasso, marca = :sorgente) =
    Fase(:secerni, _regola_secerni, (sostanza = sostanza, tasso = Float64(tasso), marca = marca),
         _fine_secerni; scrive = (:secrezioni,))
