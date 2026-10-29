# Estrusione (motore/09_estrusione.jl): una cellula apoptotica esce dal foglio.

const P_ESTRUSIONE = ParametriMeccanici(K = 0.1, η = 1.0, lunghezza_minima = 0.05,
    k_contatto = 0.0, d_contatto = 0.5)

# Foglio periodico 6 × 6 in cui la cellula 15 entra in apoptosi: il suo volume bersaglio
# scende da 1000 a 1 µm³ ("morire: apoptosi", repertorio della specifica).
function foglio_apoptosi()
    ints = fill(PRISMA, 36)
    ints[15] = Intenzioni3D(1.0, 1.0, 1.0, 0.0, 1.0)
    return foglio_periodico(ints, 6, 6; lato = 5.0, altezza = 15.0)
end

topologia_estrusione(t) = rilassa_con_topologia!(t, P_ESTRUSIONE; soglia = 0.3,
    lunghezza_nuova = 1.0, area_soglia = 0.5, tolleranza = 1e-9, cicli_massimi = 30)

@testset "estrusione" begin
    @testset "interruttore: senza topologia la cellula apoptotica si blocca" begin
        t = foglio_apoptosi()
        @test !minimizza!(t, P_ESTRUSIONE; tolleranza = 1e-9, iterazioni_massime = 300).convergito
        @test length(poligoni(t)[15]) == 6
    end

    @testset "la cellula apoptotica diventa un triangolo ed esce dal foglio" begin
        t = foglio_apoptosi()
        nv = (length(t.x) - 2) ÷ 6
        r = topologia_estrusione(t)
        @test r.convergito
        @test length(r.scambi) == 3                      # da esagono a triangolo
        @test all(15 in s.perdono for s in r.scambi)
        @test length(r.estrusioni) == 1 && r.estrusioni[1].cellula == 15
        @test length(t.cellule) == 35
        @test (length(t.x) - 2) ÷ 6 == nv - 2            # tre vertici diventano uno
        W = poligoni(t)
        @test all(length(cellule_attorno(W, j)) == 3 for j in 1:(length(t.x) - 2) ÷ 6)
        @test all(volume_cellula(t, c) > 0 for c in eachindex(t.cellule))
        @test all(t.intenzioni[c].V0 == 1000.0 for c in eachindex(t.cellule))   # è uscita lei
        @test errore_gradiente(t, P_ESTRUSIONE) < 1e-7
    end

    @testset "l'estrusione libera energia, poco, e il tessuto si richiude" begin
        # Fino al triangolo con i soli scambi, poi l'estrusione a mano. Il triangolo ha apice e
        # base minuscoli ma è ancora alto ~15 µm: con lui spariscono pareti per ~26 nN·µm, e la
        # fusione dei tre vertici sposta anche le tre pareti fra le vicine. È un evento
        # discreto: il salto non è zero, ma è una caduta, piccola rispetto al totale, e il
        # rilassamento successivo abbassa ancora l'energia.
        t = foglio_apoptosi()
        rilassa_con_scambi!(t, P_ESTRUSIONE; soglia = 0.3, lunghezza_nuova = 1.0,
                            tolleranza = 1e-9, cicli_massimi = 30)
        @test length(poligoni(t)[15]) == 3
        @test area_apicale(t, 15) < 0.5 && Morfogenetica._area_basale(t, 15) < 0.5
        E = energia(t, P_ESTRUSIONE)
        @test estrudi!(t, 15) isa Estrusione
        E1 = energia(t, P_ESTRUSIONE)
        @test E1 < E
        @test (E - E1) / E < 0.005
        @test minimizza!(t, P_ESTRUSIONE; tolleranza = 1e-9, iterazioni_massime = 200).convergito
        @test energia(t, P_ESTRUSIONE) < E1
    end

    @testset "si estrude solo un triangolo" begin
        t = foglio_apoptosi()
        @test estrudi!(t, 15) === nothing
        @test length(t.cellule) == 36
    end

    @testset "risultato negativo #15: una cellula che chiude solo l'apice non esce" begin
        # Anello apicale forte: l'apice si chiude (≈ 0.015 µm²) ma la base si allarga
        # (≈ 175 µm²), la cellula a bottiglia. Nel guscio a prismi non la si può estrudere:
        # sparirebbe di colpo una base grande.
        ints = fill(PRISMA, 36)
        ints[15] = Intenzioni3D(1.0, 1.0, 1.0, 8.0, 1000.0)
        t = foglio_periodico(ints, 6, 6; lato = 5.0, altezza = 15.0)
        minimizza!(t, P_ESTRUSIONE; tolleranza = 1e-9, iterazioni_massime = 100)
        @test area_apicale(t, 15) < 0.5
        @test Morfogenetica._area_basale(t, 15) > 100
        @test isempty(estrusioni!(t; area_soglia = 0.5))
        @test length(t.cellule) == 36
    end
end
