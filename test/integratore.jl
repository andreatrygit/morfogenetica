# Integratore (motore/06_integratore.jl): Newton smorzato per l'equilibrio, Eulero implicito
# per la dinamica. Deve ritrovare gli stessi equilibri e, per dt → 0, la traiettoria di Eulero.

minimizza_prova!(t, p = P_PROVA) = minimizza!(t, p; tolleranza = 1e-10, iterazioni_massime = 3000)

@testset "integratore" begin
    @testset "Newton: stessi equilibri di Eulero, molte meno iterazioni" begin
        lieve = Intenzioni(1.05, 1.0, 1.0, 100.0)
        a = striscia(fill(lieve, 10); larghezza = 7.0, altezza = 14.0)
        b = striscia(fill(lieve, 10); larghezza = 7.0, altezza = 14.0)
        # Eulero fino a 1e-12 nN: con 1e-9 lo scarto residuo nella piegatura, che è morbida,
        # vale da solo ~1e-8 di forma
        ra = rilassa!(a, P_PROVA; dt = 0.2, tolleranza = 1e-12, passi_massimi = 5_000_000)
        rb = minimizza_prova!(b)
        @test ra.convergito && rb.convergito
        # una striscia libera può traslare e ruotare senza costo: si confronta la forma, non
        # le posizioni
        for c in 1:10
            @test abs(semiangolo(a, c) - semiangolo(b, c)) < 1e-8
            @test abs(altezza(a, c) - altezza(b, c)) < 1e-8
            @test abs(area(a.x, a.cellule[c]) - area(b.x, b.cellule[c])) < 1e-7
        end
        @test isapprox(energia(a, P_PROVA), energia(b, P_PROVA); rtol = 1e-12)
        @test rb.iterazioni < ra.passi / 100
    end

    @testset "Newton: F2 anello esatto" begin
        t = anello(fill(CUNEO, 12); r_apicale = 15.0, r_basale = 25.0)
        @test minimizza_prova!(t).convergito
        α = π / 12
        for c in 1:12
            A = area(t.x, t.cellule[c])
            @test isapprox(altezza(t, c), cuneo_altezza(CUNEO.Γa, CUNEO.Γb, CUNEO.Γl, A, α); rtol = 1e-9)
        end
    end

    @testset "Newton: F3 con strisce che Eulero non riesce a rilassare" begin
        # con Eulero, 40 cellule non convergono in 5 milioni di passi (docs/piano.md, sezione 6)
        lieve = Intenzioni(1.05, 1.0, 1.0, 100.0)
        α_teor = cuneo_semiangolo(lieve.Γa, lieve.Γb, lieve.Γl)
        t = striscia(fill(lieve, 40); larghezza = 7.0, altezza = 14.0)
        @test minimizza_prova!(t).convergito
        @test autointersezioni(t) == 0
        @test abs(semiangolo(t, 20) - α_teor) < 1e-9 * α_teor
    end

    @testset "Newton: F5 contro tyssue" begin
        t = striscia(intenzioni_confronto(9:12); larghezza = 7.0, altezza = 14.0)
        @test minimizza!(t, P_CONFRONTO; tolleranza = 1e-10, iterazioni_massime = 3000).convergito
        confronta(misure(t), leggi_tyssue("striscia_k4"))
    end

    @testset "Newton con il contatto: non si attraversa" begin
        t = arrotolata()
        @test minimizza_prova!(t, P_CONTATTO).convergito
        @test autointersezioni(t) == 0
        @test energia_contatto(t, P_CONTATTO) > 0
    end

    @testset "passo implicito: l'energia non sale mai, anche con passi lunghi" begin
        # Eulero esplicito con questi parametri è stabile solo sotto dt ≈ 0.4 min
        t = striscia(fill(Intenzioni(1.05, 1.0, 1.0, 100.0), 20); larghezza = 7.0, altezza = 14.0)
        E = energia(t, P_PROVA)
        for _ in 1:20
            passo_implicito!(t, P_PROVA, 50.0)
            En = energia(t, P_PROVA)
            @test En <= E + 1e-12 * abs(E)
            E = En
        end
    end

    @testset "passo implicito: per dt → 0 ritrova Eulero, al primo ordine" begin
        T = 20.0
        riferimento = tessuto_disturbato(6)
        for _ in 1:round(Int, T / 0.0125)
            passo!(riferimento, P_PROVA, 0.0125)
        end
        errore(dt) = begin
            t = tessuto_disturbato(6)
            evolvi!(t, P_PROVA; durata = T, dt)
            maximum(abs.(t.x .- riferimento.x))
        end
        e2, e1, e05 = errore(2.0), errore(1.0), errore(0.5)
        @test e1 < e2 && e05 < e1
        @test 1.6 < e2 / e1 < 2.5
        @test 1.6 < e1 / e05 < 2.5
    end
end
