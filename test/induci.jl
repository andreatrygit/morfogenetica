# induci: un gene bersaglio che si accumula con il segnale (R5, la catena a scala vera).
#
#   induci(; sostanza, voce, β, γ, K, n)   primitive/induci.jl: il prodotto g (stato[voce]) segue
#                                          dg/dt = β·H(s) − γ·g, H(s) = sⁿ/(Kⁿ + sⁿ), con s la
#                                          concentrazione alla propria posizione, integrato in
#                                          modo esatto nella passata (s costante nella passata)
#   costringi(; Λ, soglia, voce)           la stessa costrizione, comandata da stato[voce] invece
#                                          che da una concentrazione
#
# Ogni prova è una funzione che restituisce un Bool. Scritte prima del codice (commit rosso, con
# @test_broken), passano dal commit verde.

const GENE = (β = 0.02, γ = log(2) / 60, K = 1.3, n = 4)     # 1/min; t½ = 60 min

function embrione_segnale(s; dt = 5.0, fase = induci(; sostanza = :segnale, voce = :gene, GENE...))
    t = foglio_programma(1)
    campo = Campo(Sostanza(:segnale, 0.0, 0.0); origine = (-20.0, -20.0, -5.0), h = 5.0, n = (8, 8, 5))
    campo.c .= s
    return Embrione(t, [fase]; parametri = P_PROGRAMMA, dt, campi = Dict(:segnale => campo),
                    interruttori = Interruttori(; SPENTO...))
end

function prova_accumula_e_decade_esattamente()
    s = 3.0
    e = embrione_segnale(s)
    for _ in 1:24                                        # 2 ore
        passata!(e)
    end
    H = s^GENE.n / (GENE.K^GENE.n + s^GENE.n)
    atteso = GENE.β / GENE.γ * H * (1 - exp(-GENE.γ * 120))
    all(isapprox(c.stato[:gene], atteso; rtol = 1e-12) for c in e.cellule) || return false
    # spento il segnale, decade con il tempo di dimezzamento
    e.campi[:segnale].c .= 0.0
    for _ in 1:12                                        # 1 ora
        passata!(e)
    end
    return all(isapprox(c.stato[:gene], atteso * exp(-GENE.γ * 60); rtol = 1e-12) for c in e.cellule)
end

function prova_un_impulso_breve_non_basta()
    # una passata di segnale alto (un picco di rumore) porta il prodotto molto sotto quello di
    # regime: la cellula filtra le fluttuazioni brevi
    e = embrione_segnale(5.0)
    passata!(e)
    regime = GENE.β / GENE.γ * 5.0^GENE.n / (GENE.K^GENE.n + 5.0^GENE.n)
    return all(c.stato[:gene] < 0.1 * regime for c in e.cellule)
end

function prova_costringi_dal_gene()
    f = costringi(; Λ = 2.0, soglia = 0.5, voce = :gene)
    p(g) = Percezione(1, 0, 0, 0.5, Dict(:gene => g), PRISMA,
                      (volume = 1.0, area_apicale = 1.0, area_basale = 1.0), Vicina[],
                      Dict{Symbol,Float64}(), Dict{Symbol,Float64}())
    return f.regola(p(0.7), f.parametri).intenzioni[:Λa] == 2.0 &&
           f.regola(p(0.3), f.parametri).intenzioni[:Λa] == PRISMA.Λa
end

@testset "R5 induci: il gene bersaglio" begin
    @test prova_accumula_e_decade_esattamente()
    @test prova_un_impulso_breve_non_basta()
    @test prova_costringi_dal_gene()
end
