# R5 (docs/roadmap.md): la prima forma, l'invaginazione di una placca.
#
#   secerni(; sostanza, tasso, marca = :sorgente)   primitiva (primitive/secerni.jl): secerne
#                                                   `tasso` se stato[marca] == 1
#   insieme(a, b; nome)                             ‖ dell'algebra (programmi/composizione.jl):
#                                                   due regole attive insieme; se toccano la
#                                                   stessa intenzione è un conflitto, che non ha
#                                                   ancora una legge, e si rifiuta
#   con_fine(f, fine)                               la stessa fase con un altro predicato di fine
#   invaginazione(; Λ, soglia, tasso, rilascio = nothing)
#                                                   il programma (programmi/invaginazione.jl):
#                                                   l'organizzatore secerne, chi legge il segnale
#                                                   sopra soglia costringe; con `rilascio` = T la
#                                                   costrizione si spegne dopo T passate
#                                                   sull'orologio di ogni cellula
#
# Le prove sono comportamenti (prove/attesi.md, gastrulazione; quaderno, invaginazione), non
# valori: la placca scende dal lato basale, più forte la costrizione più profonda, rilasciata si
# riapre (la piega è elastica), mantenuta resta. L'organizzatore è una condizione iniziale (lo
# stato :sorgente di una cellula): l'elezione che lo sceglie senza seme è della catena di
# gastrulazione, dopo.
#
# Ogni prova è una funzione che restituisce un Bool. Scritte prima del codice (commit rosso, con
# @test_broken), passano dal commit verde. Corrette nel verde tre scelte della prova: il foglio
# parte all'equilibrio, il passo meccanico è di 10 minuti (la piega si muove con un tempo di
# circa 35 minuti), e "mantenuta" si misura contro l'equilibrio statico.

const P_FORMA = ParametriMeccanici(K = 0.1, η = 1.0, lunghezza_minima = 0.05,
                                   k_contatto = 0.0, d_contatto = 0.5)
const ORGANIZZATORE = 28            # una cellula interna del foglio periodico 8 × 8

# Il foglio parte all'equilibrio: i prismi costruiti hanno 974 µm³ contro un bersaglio di
# 1000, e senza questo la prima passata è un assestamento che si confonde con la forma.
const _FOGLIO = Ref{Any}(nothing)
function foglio_forma()
    if _FOGLIO[] === nothing
        t = foglio_periodico(fill(PRISMA, 64), 8, 8; lato = 5.0, altezza = 15.0)
        minimizza!(t, P_FORMA; tolleranza = 1e-9, iterazioni_massime = 300)
        _FOGLIO[] = t
    end
    t = deepcopy(_FOGLIO[])
    Lx, Ly = t.x[end-1], t.x[end]
    campo = Campo(Sostanza(:segnale, 6.0, 0.1); origine = (0.0, 0.0, -6.0),
                  h = (Lx / 24, Ly / 24, 3.0), n = (24, 24, 9),
                  bordi = (bordo_periodico, bordo_periodico, bordo_chiuso))
    return t, campo
end

function embrione_forma(programma; meccanica = true)
    t, campo = foglio_forma()
    stati = [Dict{Symbol,Float64}(:sorgente => (c == ORGANIZZATORE ? 1.0 : 0.0)) for c in eachindex(t.cellule)]
    return Embrione(t, programma; parametri = P_FORMA, dt = 10.0, campi = Dict(:segnale => campo),
                    stati, interruttori = Interruttori(; meccanica, topologia = false))
end

# La soglia fra il primo e il secondo anello attorno all'organizzatore, letta dall'osservatore a
# regime (prepara la prova: la placca sono l'organizzatore e le sue sei vicine).
const _SOGLIA = Ref(NaN)
function soglia_placca()
    isnan(_SOGLIA[]) || return _SOGLIA[]
    e = embrione_forma(invaginazione(; Λ = 0.0, soglia = Inf, tasso = 20.0); meccanica = false)
    for _ in 1:12
        passata!(e)
    end
    t = e.tessuto
    c = [concentrazione(e.campi[:segnale], centro_cellula(t, k)) for k in eachindex(t.cellule)]
    primo = vicine_di(t, ORGANIZZATORE)
    secondo = setdiff(union((vicine_di(t, v) for v in primo)...), primo, [ORGANIZZATORE])
    _SOGLIA[] = (minimum(c[primo]) + maximum(c[secondo])) / 2
    return _SOGLIA[]
end

placca(t) = vcat([ORGANIZZATORE], vicine_di(t, ORGANIZZATORE))

# Profondità della placca: quota media della faccia apicale delle cellule della placca meno
# quella delle cellule lontane (negativa: verso il lato basale).
function _quota_apicale(t, c)
    cel = t.cellule[c]
    k = findfirst(fi -> t.facce[fi].etichetta == apicale, cel.facce)
    return Morfogenetica._baricentro(Morfogenetica._vertici(t, t.x, t.facce[cel.facce[k]],
                                                            cel.traslazioni[k]))[3]
end
function profondita(t)
    P = placca(t)
    vicini = union(P, (vicine_di(t, v) for v in P)...)
    lontane = setdiff(eachindex(t.cellule), vicini)
    return sum(_quota_apicale(t, c) for c in P) / length(P) -
           sum(_quota_apicale(t, c) for c in lontane) / length(lontane)
end

function corsa_forma(Λ; passate = 12, rilascio = nothing)
    e = embrione_forma(invaginazione(; Λ, soglia = soglia_placca(), tasso = 20.0, rilascio))
    quote_ = Float64[]
    for _ in 1:passate
        passata!(e)
        push!(quote_, profondita(e.tessuto))
    end
    return e, quote_
end

const _CORSE = Dict{Any,Any}()
corsa(Λ; kw...) = get!(() -> corsa_forma(Λ; kw...), _CORSE, (Λ, sort(collect(kw))))

function prova_la_placca_si_invagina()
    e, q = corsa(3.0; passate = 24)
    attive = [c for c in eachindex(e.cellule) if e.cellule[c].intenzioni.Λa > 0]
    return sort(attive) == sort(placca(e.tessuto)) && q[12] < -2.0 &&
           isempty(problemi_topologici(e.tessuto)) &&
           all(volume_cellula(e.tessuto, c) > 0 for c in eachindex(e.cellule))
end

function prova_senza_costrizione_resta_piatta()
    _, q = corsa(0.0)
    return maximum(abs, q) < 0.05
end

function prova_piu_forte_piu_profonda()
    _, q1 = corsa(1.5); _, q3 = corsa(3.0; passate = 24)
    return q3[12] < q1[12] < 0
end

function prova_mantenuta_resta()
    # `sempre mantieni(costringi)`: la piega non torna mai indietro e va all'equilibrio statico
    # con la costrizione accesa (la piega si muove con un tempo di circa 35 minuti)
    e, q = corsa(3.0; passate = 24)
    all(q[k+1] <= q[k] + 1e-9 for k in 1:length(q)-1) || return false
    t = deepcopy(e.tessuto)
    minimizza!(t, P_FORMA; tolleranza = 1e-9, iterazioni_massime = 300)
    return abs(q[end] - profondita(t)) < 0.05 * abs(profondita(t))
end

function prova_rilasciata_si_riapre()
    # la piega è elastica: spenta la costrizione dopo 12 passate, la placca torna piatta
    e, q = corsa(3.0; passate = 24, rilascio = 12)
    return q[12] < -2.0 && abs(q[end]) < 0.1 * abs(q[12]) &&
           all(c.intenzioni.Λa == 0 for c in e.cellule)
end

function prova_insieme_rifiuta_i_conflitti()
    # due costringi scrivono entrambe Λa: rifiutate quando si compongono (test/composizione.jl)
    a = costringi(; Λ = 1.0, sostanza = :s, soglia = 0.0)
    b = costringi(; Λ = 2.0, sostanza = :s, soglia = 0.0)
    conflitto = try
        insieme(a, b; nome = :due); false
    catch err
        err isa ArgumentError
    end
    # senza conflitto: secerne e costringe insieme
    g = insieme(secerni(; sostanza = :s, tasso = 3.0), a; nome = :ok)
    u = g.regola(Percezione(1, 0, 0, 0.5, Dict(:sorgente => 1.0), PRISMA,
                            (volume = 1.0, area_apicale = 1.0, area_basale = 1.0), Vicina[],
                            Dict{Symbol,Float64}(), Dict(:s => 1.0)), g.parametri)
    return conflitto && u.secrezioni[:s] == 3.0 && u.intenzioni[:Λa] == 1.0
end

@testset "R5 invaginazione" begin
    @test prova_insieme_rifiuta_i_conflitti()
    @test prova_la_placca_si_invagina()
    @test prova_senza_costrizione_resta_piatta()
    @test prova_piu_forte_piu_profonda()
    @test prova_mantenuta_resta()
    @test prova_rilasciata_si_riapre()
end
