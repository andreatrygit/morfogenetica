# La catena dell'organizzatore e dell'invaginazione a scala vera (docs/roadmap.md, R5): il primo
# pezzo della gastrulazione compilata (specifica, sezione 5), con meccanismi plausibili al posto
# dell'elezione a biglietti e del campo di distanza.
#
#   fase 1  organizza   attiva_inibisci ‖ induci: Nodal e Lefty fanno nascere senza seme una
#                       macchia di Nodal (giorni, circa 1 mm); ogni cellula accumula il prodotto di
#                       un gene bersaglio. Finisce, per la cellula, quando il prodotto supera
#                       `impegno` (il destino di organizzatore); l'onda porta poi tutte le cellule
#                       nella fase 2.
#   fase 2  invagina    lo stesso, più costringi sul prodotto del gene: le cellule impegnate
#                       tendono l'anello apicale, e la regione si invagina.
#
# La soglia è sul prodotto del gene, non sulla concentrazione: la cellula integra il segnale
# (Dubrulle et al. 2015), così i picchi di rumore e le cellule al bordo della macchia non
# accendono e spengono la costrizione.

export gastrula

"""
    gastrula(; Da, ka, Dh, kh, κ, rumore = 0.2, c0 = 1.0, β, γ, K, n, impegno, Λ, dt)

La tabella di fasi organizza ; invagina. `ka`, `kh`, `κ`, `rumore`, `c0` sono quelli di
`attiva_inibisci` (i campi di Nodal e Lefty, con `Da` e `Dh`, sono dell'embrione); `β`, `γ`,
`K`, `n` quelli di `induci`; `impegno` la soglia del prodotto del gene; `Λ` (nN) la tensione
dell'anello apicale delle cellule impegnate; `dt` la passata (min).
"""
function gastrula(; ka, kh, κ, rumore = 0.2, c0 = 1.0, β, γ, K, n, impegno, Λ, dt, Da = nothing,
                  Dh = nothing)
    ai = attiva_inibisci(; attivatore = :nodal, inibitore = :lefty, ka, kh, κ, rumore, c0)
    ge = induci(; sostanza = :nodal, voce = :gene, β, γ, K, n, dt)
    organizza = con_fine(insieme(ai, ge; nome = :organizza),
                         (p, par) -> get(p.stato, :gene, 0.0) > impegno)
    invagina = insieme(insieme(ai, ge; nome = :segnale), costringi(; Λ, soglia = impegno, voce = :gene);
                       nome = :invagina)
    return [organizza, invagina]
end
