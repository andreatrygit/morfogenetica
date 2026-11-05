# attiva_inibisci: attivazione locale e inibizione a lungo raggio (Gierer e Meinhardt). È il modo
# plausibile di far nascere un organizzatore senza seme: da un rumore piccolo, una sola macchia.
#
# Regola, sui canali "dentro di sé" (volume, rumore di espressione) e "misura assoluta" (le due
# concentrazioni alla propria posizione): con â = a/c0 e ĥ = h/c0, la cellula secerne
#
#   attivatore   ka c0 (1 + κ) â² / ((1 + κ â²) ĥ) · V · r
#   inibitore    kh c0 â² · V · r
#
# dove V è il suo volume e r = 1 + rumore·(2u − 1), con u il suo numero a caso della passata (il
# rumore di espressione). Le costanti sono scelte perché lo stato uniforme sia a = h = c0 quando
# le cellule riempiono il tessuto. κ satura l'autoattivazione: senza, con l'inibitore che si
# elimina più lentamente dell'attivatore (come Lefty rispetto a Nodal), lo stato uniforme
# oscilla invece di formare una macchia. Predicato di fine: nessuno.
#
# Attenzione al passo della passata: le secrezioni restano fisse per tutta la passata, e per
# questo circuito è un passo esplicito, stabile solo se la passata è più corta di 2δ/ω² (δ e ω
# smorzamento e frequenza dello stato uniforme linearizzato). Con Cyclops/Lefty2 e κ = 0.35 il
# limite è 12.9 minuti, con κ = 0.29 1.9 minuti (risultato negativo #17).
#
# Meccanismo: Nodal induce la propria espressione e quella di Lefty; Lefty, più diffusivo,
#   inibisce la via di Nodal a distanza. Nodal e Lefty sono una coppia attivatore–inibitore dal
#   riccio di mare al topo; la loro diffusività diversa è la condizione di Turing.
# Fonte: Müller et al., "Differential diffusivity of Nodal and Lefty underlies a
#   reaction-diffusion patterning system", Science 336, 721–724 (2012),
#   doi:10.1126/science.1221920; Gierer e Meinhardt, "A theory of biological pattern
#   formation", Kybernetik 12, 30–39 (1972).

export attiva_inibisci

function _regola_attiva_inibisci(p::Percezione, par)
    a = max(get(p.concentrazioni, par.attivatore, 0.0), 0.0) / par.c0
    h = max(get(p.concentrazioni, par.inibitore, 0.0), 1e-9 * par.c0) / par.c0
    r = 1 + par.rumore * (2 * p.caso - 1)
    V = p.misure.volume
    qa = par.ka * par.c0 * (1 + par.κ) * a^2 / ((1 + par.κ * a^2) * h) * V * r
    qh = par.kh * par.c0 * a^2 * V * r
    return Uscita(; secrezioni = Dict{Symbol,Float64}(par.attivatore => qa, par.inibitore => qh))
end

_fine_attiva_inibisci(p::Percezione, par) = false

"""
    attiva_inibisci(; attivatore, inibitore, ka, kh, κ, rumore, c0)

La fase in cui ogni cellula produce l'`attivatore` in modo autocatalitico e saturante (κ),
frenato dall'`inibitore`, e l'inibitore in proporzione al quadrato dell'attivatore. `ka` e
`kh` (1/min) sono le costanti di eliminazione delle due sostanze (le stesse dei loro campi),
`c0` la concentrazione dello stato uniforme, `rumore` l'ampiezza relativa del rumore di
espressione.
"""
attiva_inibisci(; attivatore, inibitore, ka, kh, κ, rumore, c0) =
    Fase(:attiva_inibisci, _regola_attiva_inibisci,
         (attivatore = attivatore, inibitore = inibitore, ka = Float64(ka), kh = Float64(kh),
          κ = Float64(κ), rumore = Float64(rumore), c0 = Float64(c0)),
         _fine_attiva_inibisci; scrive = (:secrezioni,))
