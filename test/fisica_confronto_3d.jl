# F7 — confronto esterno 3D con tyssue: un pezzo di foglio di 37 prismi con il bordo bloccato,
# in cui la cellula centrale (o il centro e il primo anello) stringe l'anello apicale.
# Geometria e tensioni esportate da confronto/esporta_guscio.jl; i numeri di tyssue sono in
# confronto/risultati/guscio_<scenario>_tyssue.csv (confronto/tyssue_guscio.py).

include(joinpath(Morfogenetica.RADICE, "confronto", "scenari_guscio.jl"))

function leggi_tyssue_guscio(nome)
    righe = readlines(joinpath(Morfogenetica.RADICE, "confronto", "risultati",
                               "guscio_$(nome)_tyssue.csv"))
    campi = [parse.(Float64, split(r, ",")) for r in righe[2:end]]
    return (area_apicale = [c[2] for c in campi], altezza = [c[3] for c in campi],
            volume = [c[4] for c in campi], z_apicale = [c[5] for c in campi],
            energia = campi[1][6])
end

@testset "F7 confronto 3D con tyssue" begin
    for nome in SCENARI_GUSCIO
        @testset "$nome" begin
            t = scenario_guscio(nome)
            fermi = copy(t.x[repeat(t.fissi, inner = 3)])
            @test minimizza!(t, P_GUSCIO_CONFRONTO; tolleranza = 1e-9,
                             iterazioni_massime = 200).convergito
            # i vertici del bordo non si sono mossi
            @test t.x[repeat(t.fissi, inner = 3)] == fermi
            loro = leggi_tyssue_guscio(nome)
            n = length(t.cellule)
            za(c) = Morfogenetica._baricentro(Morfogenetica._vertici(t, t.x,
                        Morfogenetica._faccia_di(t, c, apicale)))[3]
            @test maximum(abs(area_apicale(t, c) - loro.area_apicale[c]) for c in 1:n) < 1e-5
            @test maximum(abs(altezza_prisma(t, c) - loro.altezza[c]) for c in 1:n) < 1e-5
            @test maximum(abs(volume_cellula(t, c) - loro.volume[c]) for c in 1:n) < 1e-4
            @test maximum(abs(za(c) - loro.z_apicale[c]) for c in 1:n) < 1e-5
            @test isapprox(energia(t, P_GUSCIO_CONFRONTO), loro.energia; rtol = 1e-10)
        end
    end
end
