# L'invaginazione di una placca (docs/roadmap.md, R5): la prima forma.
#
#   fase 1  invagina   secerni ‖ costringi: l'organizzatore secerne il segnale; ogni cellula
#                      che lo legge sopra soglia alla propria posizione tende l'anello apicale.
#                      Senza `rilascio` non finisce mai: è `sempre mantieni(costringi)`.
#   fase 2  rilascia   (solo con `rilascio` = T) la costrizione si spegne; la fase 1 finisce
#                      dopo T passate sull'orologio di ogni cellula.
#
# La placca non è disegnata: sono le cellule a cui arriva abbastanza segnale. L'organizzatore
# (stato :sorgente = 1) è per ora una condizione iniziale; l'elezione che lo sceglie senza seme
# è della catena di gastrulazione.

export invaginazione

"""
    invaginazione(; Λ, soglia, tasso, sostanza = :segnale, rilascio = nothing)

La tabella di fasi dell'invaginazione: tensione dell'anello apicale `Λ` (nN) dove `sostanza`
supera `soglia`, secreta dall'organizzatore a `tasso` (quantità/min); con `rilascio` = T la
costrizione si spegne dopo T passate.
"""
function invaginazione(; Λ, soglia, tasso, sostanza = :segnale, rilascio = nothing)
    attiva = insieme(secerni(; sostanza, tasso), costringi(; Λ, sostanza, soglia); nome = :invagina)
    rilascio === nothing && return [attiva]
    spenta = insieme(secerni(; sostanza, tasso), costringi(; Λ = 0.0, sostanza, soglia);
                     nome = :rilascia)
    return [con_fine(attiva, (p, par) -> p.orologio >= rilascio - 1), spenta]
end
