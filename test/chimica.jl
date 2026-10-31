# R3 (docs/roadmap.md): la chimica. Una sostanza che diffonde su una griglia cartesiana 3D, con
# degradazione, secrezione e assorbimento delle cellule alla propria posizione.
#
#   Sostanza(nome, D, k)                     D in µm²/min, k in 1/min
#   Campo(sostanza; origine, h, n, bordi)    h per asse (µm) o uno solo; bordi bordo_chiuso o
#                                            bordo_periodico per asse; c[i, j, l] concentrazioni
#   centro_voxel(campo, i, j, l)             centro di un voxel (µm)
#   massa(campo)                             Σ c · volume del voxel
#   concentrazione(campo, p)                 lettura alla posizione p (interpolazione trilineare)
#   Scambio(p, secrezione, assorbimento)     quantità/min; µm³/min (assorbimento del primo ordine)
#   avanza!(campo, Δt, scambi; passo)        → (secreto, assorbito, passi)
#   passo_massimo(campo)                     il passo esplicito più lungo che resta stabile
#   centro_cellula(t, c)                     dove una cellula del guscio secerne e legge
#
# Le prove sono analitiche: conservazione, varianza 6Dt, decadimento esponenziale, lunghezza
# √(D/k), bilancio stazionario, sorgente puntiforme, convergenza con il passo. Parametri di
# riferimento (Dpp nel disco dell'ala di Drosophila, Kicheva et al., Science 315, 521, 2007):
# D = 0.1 µm²/s = 6 µm²/min e λ ≈ 20 µm, quindi k = D/λ² = 0.015/min.
#
# Ogni prova è una funzione che restituisce un Bool. Commit rosso: marcate @test_broken finché
# il codice (motore/12_chimica.jl) non c'è.

const DPP = (D = 6.0, k = 0.015)          # µm²/min, 1/min: λ = 20 µm

function _momenti(campo)
    M = 0.0; g = (0.0, 0.0, 0.0)
    for I in CartesianIndices(campo.c)
        p = centro_voxel(campo, Tuple(I)...)
        M += campo.c[I]; g = g .+ campo.c[I] .* p
    end
    g = g ./ M
    v = 0.0
    for I in CartesianIndices(campo.c)
        p = centro_voxel(campo, Tuple(I)...)
        v += campo.c[I] * sum(abs2, p .- g)
    end
    return M, g, v / M
end

function prova_massa_conservata()
    rng = MersenneTwister(7)
    for bordi in ((bordo_chiuso, bordo_chiuso, bordo_chiuso),
                  (bordo_periodico, bordo_periodico, bordo_chiuso))
        campo = Campo(Sostanza(:s, DPP.D, 0.0); origine = (0.0, 0.0, 0.0), h = 2.0, n = (10, 8, 6), bordi)
        campo.c .= rand(rng, size(campo.c)...)
        M = massa(campo)
        avanza!(campo, 30.0)
        isapprox(massa(campo), M; rtol = 1e-12) || return false
    end
    return true
end

function prova_varianza_cresce_come_6Dt()
    campo = Campo(Sostanza(:s, DPP.D, 0.0); origine = (0.0, 0.0, 0.0), h = 1.0, n = (41, 41, 41))
    campo.c[21, 21, 21] = 1.0
    _, _, v0 = _momenti(campo)
    t = 0.5
    avanza!(campo, t)
    _, g, v = _momenti(campo)
    return v0 == 0 && isapprox(v, 6 * DPP.D * t; rtol = 1e-9) &&
           all(isapprox.(g, centro_voxel(campo, 21, 21, 21); atol = 1e-9))
end

function prova_degradazione_esatta()
    campo = Campo(Sostanza(:s, DPP.D, DPP.k); origine = (0.0, 0.0, 0.0), h = 2.0, n = (6, 6, 6))
    campo.c .= 3.0
    avanza!(campo, 50.0)
    return all(isapprox.(campo.c, 3.0 * exp(-DPP.k * 50.0); rtol = 1e-12))
end

# Profilo stazionario in 1D da una sorgente all'estremo: c(x) ∝ exp(-x/λ), λ = √(D/k).
function _lambda_misurata(h)
    n = round(Int, 200 / h)
    campo = Campo(Sostanza(:s, DPP.D, DPP.k); origine = (0.0, 0.0, 0.0), h = h, n = (n, 1, 1))
    sorgente = [Scambio(centro_voxel(campo, 1, 1, 1), 1.0, 0.0)]
    avanza!(campo, 15 / DPP.k, sorgente)
    i = round(Int, 40 / h)
    return h / log(campo.c[i, 1, 1] / campo.c[i + 1, 1, 1])
end

function prova_lunghezza_di_decadimento()
    λ = sqrt(DPP.D / DPP.k)
    e1 = abs(_lambda_misurata(1.0) - λ) / λ
    e2 = abs(_lambda_misurata(2.0) - λ) / λ
    return e1 < 1e-3 && 3 < e2 / e1 < 5            # secondo ordine nella griglia
end

function prova_bilancio_stazionario()
    # allo stato stazionario, con bordi chiusi o periodici, ciò che entra esce per degradazione
    k = 0.1; q = 5.0
    campo = Campo(Sostanza(:s, DPP.D, k); origine = (0.0, 0.0, 0.0), h = 2.0, n = (12, 12, 4),
                  bordi = (bordo_periodico, bordo_periodico, bordo_chiuso))
    avanza!(campo, 20 / k, [Scambio((7.3, 15.1, 3.2), q, 0.0)])
    return isapprox(k * massa(campo), q; rtol = 1e-6)
end

function prova_sorgente_puntiforme()
    # c(r) = q / (4π D r) · exp(-r/λ); con h = 1 µm e λ = 5 µm la tolleranza del 5% è quella
    # della discretizzazione vicino alla sorgente (regressione, non un valore della fisica)
    D = DPP.D; λ = 5.0; k = D / λ^2; q = 10.0
    campo = Campo(Sostanza(:s, D, k); origine = (0.0, 0.0, 0.0), h = 1.0, n = (41, 41, 41))
    p0 = centro_voxel(campo, 21, 21, 21)
    avanza!(campo, 14 / k, [Scambio(p0, q, 0.0)])
    for r in (5, 8)
        atteso = q / (4π * D * r) * exp(-r / λ)
        isapprox(campo.c[21 + r, 21, 21], atteso; rtol = 0.05) || return false
    end
    return true
end

function prova_secrezione_e_lettura_alla_posizione()
    # senza diffusione: la secrezione va nei voxel vicini con il baricentro nella posizione
    campo = Campo(Sostanza(:s, 0.0, 0.0); origine = (0.0, 0.0, 0.0), h = 2.0, n = (8, 8, 8))
    p = (5.3, 7.9, 6.4)
    r = avanza!(campo, 2.0, [Scambio(p, 1.5, 0.0)])
    M, g, _ = _momenti(campo)
    (isapprox(massa(campo), 3.0; rtol = 1e-12) && isapprox(r.secreto, 3.0; rtol = 1e-12)) || return false
    all(isapprox.(g, p; atol = 1e-12)) || return false
    # lettura: esatta su un campo lineare
    for I in CartesianIndices(campo.c)
        x, y, z = centro_voxel(campo, Tuple(I)...)
        campo.c[I] = 1 + 0.1x + 0.2y - 0.05z
    end
    return isapprox(concentrazione(campo, p), 1 + 0.1p[1] + 0.2p[2] - 0.05p[3]; rtol = 1e-12)
end

function prova_assorbimento_positivo_e_registrato()
    campo = Campo(Sostanza(:s, DPP.D, 0.0); origine = (0.0, 0.0, 0.0), h = 2.0, n = (8, 8, 8))
    campo.c .= 1.0
    M = massa(campo)
    r = avanza!(campo, 5.0, [Scambio((7.0, 8.2, 9.1), 0.0, 1e6)])     # un assorbitore enorme
    return minimum(campo.c) >= 0 && r.assorbito > 0 &&
           isapprox(M - massa(campo), r.assorbito; rtol = 1e-12)
end

function prova_convergenza_con_il_passo()
    function corsa(passo)
        campo = Campo(Sostanza(:s, DPP.D, 0.2); origine = (0.0, 0.0, 0.0), h = 1.0, n = (30, 1, 1))
        avanza!(campo, 3.0, [Scambio(centro_voxel(campo, 10, 1, 1), 1.0, 0.0)]; passo)
        return copy(campo.c)
    end
    dmax = passo_massimo(Campo(Sostanza(:s, DPP.D, 0.2); origine = (0.0, 0.0, 0.0), h = 1.0, n = (30, 1, 1)))
    rif = corsa(dmax / 64)
    e2 = maximum(abs, corsa(dmax / 2) .- rif)
    e4 = maximum(abs, corsa(dmax / 4) .- rif)
    instabile = try
        corsa(2dmax); false
    catch e
        e isa ArgumentError
    end
    return instabile && 1.7 < e2 / e4 < 2.3            # primo ordine nel tempo
end

function prova_cellule_secernono()
    # ogni cellula del guscio secerne alla propria posizione; la griglia copre il foglio
    q = 0.5; T = 10.0
    t = foglio_esagonale(fill(PRISMA, 19), 2; lato = 5.0, altezza = 15.0)
    campo = Campo(Sostanza(:s, DPP.D, 0.0); origine = (-30.0, -30.0, -5.0), h = 2.0, n = (30, 30, 13))
    scambi = [Scambio(centro_cellula(t, c), q, 0.0) for c in eachindex(t.cellule)]
    avanza!(campo, T, scambi)
    isapprox(massa(campo), 19 * q * T; rtol = 1e-12) || return false
    # nel foglio periodico la griglia ha la scatola come periodo, e i centri vi sono riportati
    s = foglio_periodico(fill(PRISMA, 16), 4, 4; lato = 5.0, altezza = 15.0)
    Lx, Ly = s.x[end-1], s.x[end]
    campo = Campo(Sostanza(:s, DPP.D, 0.0); origine = (0.0, 0.0, -5.0), h = (Lx / 16, Ly / 16, 2.0),
                  n = (16, 16, 13), bordi = (bordo_periodico, bordo_periodico, bordo_chiuso))
    centri = [centro_cellula(s, c) for c in eachindex(s.cellule)]
    all(0 <= p[1] < Lx && 0 <= p[2] < Ly for p in centri) || return false
    avanza!(campo, T, [Scambio(p, q, 0.0) for p in centri])
    isapprox(massa(campo), 16 * q * T; rtol = 1e-12) || return false
    p = (1.3, 2.7, 7.5)
    return isapprox(concentrazione(campo, p), concentrazione(campo, (p[1] + Lx, p[2] - Ly, p[3]));
                    rtol = 1e-12)
end

@testset "R3 chimica" begin
    @test_broken prova_massa_conservata()
    @test_broken prova_varianza_cresce_come_6Dt()
    @test_broken prova_degradazione_esatta()
    @test_broken prova_lunghezza_di_decadimento()
    @test_broken prova_bilancio_stazionario()
    @test_broken prova_sorgente_puntiforme()
    @test_broken prova_secrezione_e_lettura_alla_posizione()
    @test_broken prova_assorbimento_positivo_e_registrato()
    @test_broken prova_convergenza_con_il_passo()
    @test_broken prova_cellule_secernono()
end
