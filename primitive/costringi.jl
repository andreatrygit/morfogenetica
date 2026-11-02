# costringi: la costrizione apicale (specifica, sezione 3). Non è una forza inventata: è alzare
# la tensione dell'anello apicale della cellula. Minimizzare l'energia a volume costante con
# l'anello teso è diventare un cuneo.
#
# Regola, sui canali "dentro di sé" e "misura assoluta": se la concentrazione della sostanza
# alla mia posizione supera la soglia, la mia tensione di linea apicale è Λ; altrimenti torna
# quella a riposo. La condizione si rivaluta a ogni passata: una marca è un'ipotesi sul
# presente, non un ricordo (specifica, sezione 8). Predicato di fine: nessuno (la fase dura
# finché un'altra fase non la segue per l'onda).
#
# Cosa si chiederebbe a una cellula vera: contrarre l'anello di actomiosina apicale dove legge
# il segnale sopra una soglia (come nel solco ventrale di Drosophila).

export costringi

function _regola_costringi(p::Percezione, par)
    attiva = get(p.concentrazioni, par.sostanza, 0.0) > par.soglia
    r = p.intenzioni
    intenzioni = Intenzioni3D(r.σa, r.σb, r.σl, attiva ? par.Λ : r.Λa, r.V0)
    return Uscita(copy(p.stato), intenzioni, Dict{Symbol,Float64}(), Dict{Symbol,Float64}())
end

_fine_costringi(p::Percezione, par) = false

"""
    costringi(; Λ, sostanza, soglia)

La fase che alza la tensione dell'anello apicale a `Λ` (nN) dove la concentrazione di
`sostanza` supera `soglia`.
"""
costringi(; Λ, sostanza, soglia) =
    Fase(:costringi, _regola_costringi, (Λ = Float64(Λ), sostanza = sostanza, soglia = Float64(soglia)),
         _fine_costringi; scrive = (:Λa,))
