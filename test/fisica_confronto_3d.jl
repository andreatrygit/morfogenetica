# F7 — confronto esterno 3D con tyssue: un pezzo di foglio di 37 prismi con il bordo bloccato,
# in cui la cellula centrale (o il centro e il primo anello) stringe l'anello apicale.
# Geometria e tensioni esportate da confronto/esporta_guscio.jl; lo stato finale di tyssue
# (posizioni, gradiente, energia) è in confronto/risultati/guscio_<scenario>_tyssue_vertici.csv.
#
# Il minimizzatore di tyssue qui non converge, e il motivo è misurato (confronto/LEGGIMI.md): il
# suo gradiente non è la derivata esatta della sua energia. Dal 29 ottobre le due energie
# coincidono solo dove le facce sono piane: il nostro centro di triangolazione è il baricentro
# del contorno, il suo quello dei vertici (vedi il test "suddividere uno spigolo non cambia la
# geometria"). Il confronto quindi verifica: stessa energia sulla geometria di partenza; scarto
# piccolo, dovuto solo alla convenzione, sulle forme deformate; un nostro minimo più basso del
# punto in cui tyssue si ferma, raggiunto anche partendo da lì.

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

# Energia di tyssue sulla geometria di partenza, stampata da confronto/tyssue_guscio.py
# (confronto/LEGGIMI.md). Lì le facce sono piane e regolari.
const ENERGIA_INIZIALE_TYSSUE = Dict("anello_centrale" => 14748.680433476375,
                                     "placca" => 14876.180433476375)

@testset "F7 confronto 3D con tyssue" begin
    for nome in SCENARI_GUSCIO
        @testset "$nome" begin
            loro = stato_tyssue(nome)
            # 1. geometria di partenza, facce piane: stessa energia
            t = scenario_guscio(nome)
            @test isapprox(energia(t, P_GUSCIO_CONFRONTO), ENERGIA_INIZIALE_TYSSUE[nome]; rtol = 1e-12)
            # il bordo bloccato è rimasto fermo anche in tyssue (a meno dell'arrotondamento del CSV)
            @test maximum(abs.(loro.x[repeat(t.fissi, inner = 3)] .- t.x[repeat(t.fissi, inner = 3)])) < 1e-12

            # 2. posizioni deformate dove tyssue si ferma: le facce non sono più piane, e le due
            # energie differiscono per la convenzione del centro da cui si triangola una faccia
            # (tyssue: baricentro dei vertici; noi: baricentro del contorno, che non cambia
            # suddividendo uno spigolo). Lo scarto resta piccolo.
            t.x .= loro.x
            Et = energia(t, P_GUSCIO_CONFRONTO)
            @test abs(Et - loro.energia) / loro.energia < 1e-4
            @test errore_gradiente(t, P_GUSCIO_CONFRONTO) < 1e-7

            # 3. il nostro minimo, dalla stessa partenza di tyssue, è più basso del punto in cui
            # tyssue si ferma, misurato con la nostra energia
            a = scenario_guscio(nome)
            fermi = copy(a.x[repeat(a.fissi, inner = 3)])
            @test minimizza!(a, P_GUSCIO_CONFRONTO; tolleranza = 1e-9,
                             iterazioni_massime = 200).convergito
            @test a.x[repeat(a.fissi, inner = 3)] == fermi
            Ea = energia(a, P_GUSCIO_CONFRONTO)
            @test Ea < Et
            @test (Et - Ea) / Ea < 1e-5

            # 4. partendo da dove si ferma tyssue si arriva allo stesso minimo
            @test minimizza!(t, P_GUSCIO_CONFRONTO; tolleranza = 1e-9,
                             iterazioni_massime = 200).convergito
            @test isapprox(energia(t, P_GUSCIO_CONFRONTO), Ea; rtol = 1e-12)
            @test maximum(maximum(abs.(collect(u) .- collect(v))) for (u, v) in zip(forma(t), forma(a))) < 1e-6
        end
    end
end
