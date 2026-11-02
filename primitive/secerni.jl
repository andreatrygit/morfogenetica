# secerni: emettere un segnale (repertorio, "emettere e leggere un segnale"). La cellula marcata
# secerne la sostanza a tasso costante alla propria posizione; le altre no.
#
# Regola, sul canale "dentro di sé": se stato[marca] == 1 secerno `tasso` (quantità/min).
# Predicato di fine: nessuno.
#
# Cosa si chiederebbe a una cellula vera: esprimere e secernere una proteina segnale (come le
# cellule di un organizzatore). Chi è marcata lo decide una fase precedente, o la condizione
# iniziale.

export secerni

function _regola_secerni(p::Percezione, par)
    q = get(p.stato, par.marca, 0.0) == 1 ? par.tasso : 0.0
    return Uscita(copy(p.stato), p.intenzioni, Dict{Symbol,Float64}(par.sostanza => q),
                  Dict{Symbol,Float64}())
end

_fine_secerni(p::Percezione, par) = false

"""
    secerni(; sostanza, tasso, marca = :sorgente)

La fase in cui la cellula con stato[marca] == 1 secerne `sostanza` a `tasso` (quantità/min).
"""
secerni(; sostanza, tasso, marca = :sorgente) =
    Fase(:secerni, _regola_secerni, (sostanza = sostanza, tasso = Float64(tasso), marca = marca),
         _fine_secerni; scrive = (:secrezioni,))
