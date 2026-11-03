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
# Meccanismo: la rete di actomiosina apicale (miosina II) si contrae a impulsi e riduce l'area
#   apicale della cellula, come nel solco ventrale di Drosophila; l'attivazione dipende da un
#   segnale letto dalla cellula.
# Fonte: Martin, Kaschube e Wieschaus, "Pulsed contractions of an actin–myosin network drive
#   apical constriction", Nature 457, 495–499 (2009), doi:10.1038/nature07522.

export costringi

function _regola_costringi(p::Percezione, par)
    attiva = get(p.concentrazioni, par.sostanza, 0.0) > par.soglia
    return Uscita(; intenzioni = (Λa = attiva ? par.Λ : p.intenzioni.Λa,))
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
