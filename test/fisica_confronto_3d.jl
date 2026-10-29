# F7 — confronto esterno 3D con tyssue: un pezzo di foglio di 37 prismi con il bordo bloccato,
# in cui la cellula centrale (o il centro e il primo anello) stringe l'anello apicale.
# Geometria e tensioni esportate da confronto/esporta_guscio.jl; lo stato finale di tyssue
# (posizioni, gradiente, energia) è in confronto/risultati/guscio_<scenario>_tyssue_vertici.csv.
#
# Il minimizzatore di tyssue qui non converge, e il motivo è misurato: il suo gradiente non è
# la derivata esatta della sua energia (scarto dell'1% nel termine di volume, dello 0.1% in
# quello di superficie, esatto per l'anello apicale), probabilmente perché non deriva il
# baricentro delle facce, che conta solo quando le facce si deformano. Quindi il confronto non
# è fra due minimi: è fra le due energie sulle stesse posizioni deformate, più la verifica che
# le nostre forze siano la derivata esatta di quell'energia comune, che il nostro minimo sia più
# basso del punto in cui tyssue si ferma e che ci si arrivi anche partendo da lì.

include(joinpath(Morfogenetica.RADICE, "confronto", "scenari_guscio.jl"))

function stato_tyssue(nome)
    righe = readlines(joinpath(Morfogenetica.RADICE, "confronto", "risultati",
                               "guscio_$(nome)_tyssue_vertici.csv"))
    campi = [parse.(Float64, split(r, ",")) for r in righe[2:end]]
    x = reduce(vcat, [c[1:3] for c in campi])
    g = reduce(vcat, [c[4:6] for c in campi])
    return (x = x, g = g, energia = campi[1][7])
end

forma(t) = [(area_apicale(t, c), altezza_prisma(t, c), volume_cellula(t, c))
            for c in eachindex(t.cellule)]

@testset "F7 confronto 3D con tyssue" begin
    for nome in SCENARI_GUSCIO
        @testset "$nome" begin
            loro = stato_tyssue(nome)
            # 1. stesse posizioni deformate: stessa energia e stesse forze
            t = scenario_guscio(nome)
            # il bordo bloccato è rimasto fermo anche in tyssue (a meno dell'arrotondamento del CSV)
            @test maximum(abs.(loro.x[repeat(t.fissi, inner = 3)] .- t.x[repeat(t.fissi, inner = 3)])) < 1e-12
            t.x .= loro.x
            @test isapprox(energia(t, P_GUSCIO_CONFRONTO), loro.energia; rtol = 1e-12)
            # le nostre forze sono la derivata esatta di quell'energia comune …
            @test errore_gradiente(t, P_GUSCIO_CONFRONTO) < 1e-7
            # … quelle di tyssue no (risultato registrato: se un tyssue corretto lo smentisce,
            # il confronto va ripensato come confronto fra minimi)
            g = gradiente!(similar(t.x), t, P_GUSCIO_CONFRONTO)
            @test maximum(abs.(g .- loro.g)) > 1e-2

            # 2. il nostro minimo, dalla stessa partenza di tyssue
            a = scenario_guscio(nome)
            fermi = copy(a.x[repeat(a.fissi, inner = 3)])
            @test minimizza!(a, P_GUSCIO_CONFRONTO; tolleranza = 1e-9,
                             iterazioni_massime = 200).convergito
            @test a.x[repeat(a.fissi, inner = 3)] == fermi
            Ea = energia(a, P_GUSCIO_CONFRONTO)
            @test Ea < loro.energia                       # più basso di dove tyssue si ferma
            @test (loro.energia - Ea) / Ea < 1e-5

            # 3. partendo da dove si ferma tyssue si arriva allo stesso minimo
            @test minimizza!(t, P_GUSCIO_CONFRONTO; tolleranza = 1e-9,
                             iterazioni_massime = 200).convergito
            @test isapprox(energia(t, P_GUSCIO_CONFRONTO), Ea; rtol = 1e-12)
            @test maximum(maximum(abs.(collect(u) .- collect(v))) for (u, v) in zip(forma(t), forma(a))) < 1e-6
        end
    end
end
