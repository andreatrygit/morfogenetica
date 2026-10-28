# F5 — confronto esterno con tyssue (docs/piano.md, sezione 1).
# Stessi tessuti, stessa energia, due motori indipendenti. I numeri di tyssue sono in
# confronto/risultati/*.csv, prodotti da confronto/tyssue_laterale.py; qui non serve Python.

function leggi_tyssue(nome)
    righe = readlines(joinpath(Morfogenetica.RADICE, "confronto", "risultati", nome * ".csv"))
    campi = [parse.(Float64, split(r, ",")) for r in righe[2:end]]
    return (semiangolo = [c[2] for c in campi], altezza = [c[3] for c in campi],
            area = [c[4] for c in campi], energia = campi[1][5])
end

# Parametri condivisi con confronto/tyssue_laterale.py
const N_CONFRONTO = 20
const P_CONFRONTO = ParametriMeccanici(K = 0.01, η = 1.0, lunghezza_minima = 0.05,
    k_contatto = 0.0, d_contatto = 0.5)
intenzioni_confronto(costrette) =
    [Intenzioni(c in costrette ? 1.3 : 1.0, 1.0, 1.0, 100.0) for c in 1:N_CONFRONTO]

function misure(t)
    n = length(t.cellule)
    return (semiangolo = [semiangolo(t, c) for c in 1:n], altezza = [altezza(t, c) for c in 1:n],
            area = [area(t.x, cel) for cel in t.cellule], energia = energia(t, P_CONFRONTO))
end

function confronta(nostro, loro; ordine = 1:N_CONFRONTO)
    @test maximum(abs.(nostro.semiangolo .- loro.semiangolo[ordine])) < 1e-5
    @test maximum(abs.(nostro.altezza .- loro.altezza[ordine])) < 1e-5
    @test maximum(abs.(nostro.area .- loro.area[ordine])) < 1e-5
    @test isapprox(nostro.energia, loro.energia; rtol = 1e-10)
end

@testset "F5 confronto con tyssue" begin
    rilassa_confronto!(t) = rilassa!(t, P_CONFRONTO; dt = 0.2, tolleranza = 1e-9,
        passi_massimi = 5_000_000)

    @testset "striscia piatta, una cellula costretta" begin
        t = striscia(intenzioni_confronto([10]); larghezza = 7.0, altezza = 14.0)
        @test rilassa_confronto!(t).convergito
        confronta(misure(t), leggi_tyssue("striscia_k1"))
    end

    @testset "striscia piatta, quattro cellule costrette" begin
        t = striscia(intenzioni_confronto(9:12); larghezza = 7.0, altezza = 14.0)
        @test rilassa_confronto!(t).convergito
        confronta(misure(t), leggi_tyssue("striscia_k4"))
    end

    @testset "anello, quattro cellule costrette" begin
        t = anello(intenzioni_confronto(1:4); r_apicale = 25.0, r_basale = 39.0)
        @test rilassa_confronto!(t).convergito
        # tyssue numera le facce in senso orario, noi in senso antiorario: la nostra cellula
        # 1 + m corrisponde alla sua faccia 3 − m (da 0), cioè all'indice mod1(4 − m, N)
        ordine = [mod1(4 - (j - 1), N_CONFRONTO) for j in 1:N_CONFRONTO]
        confronta(misure(t), leggi_tyssue("anello_k4"); ordine)
    end
end
