# Prove fisiche F1–F4 del modello laterale contro le previsioni esatte (docs/piano.md).
# Si confronta usando l'area misurata: con l'incomprimibilità come penalità l'area all'equilibrio
# non è esattamente A0, e la previsione vale per l'area che c'è.

rilassa_prova!(t, p = P_PROVA) =
    rilassa!(t, p; dt = 0.2, tolleranza = 1e-9, passi_massimi = 3_000_000)

@testset "fisica laterale" begin
    Γa, Γb, Γl = CUNEO.Γa, CUNEO.Γb, CUNEO.Γl

    @testset "F1 foglio piatto: h/w = Σ/Γl" begin
        errori = Float64[]
        for n in (10, 20, 40)
            t = striscia(fill(Intenzioni(1.0, 1.0, 1.0, 100.0), n); larghezza = 7.0, altezza = 14.0)
            @test rilassa_prova!(t).convergito
            c = n ÷ 2
            @test abs(semiangolo(t, c)) < 1e-9
            @test autointersezioni(t) == 0
            push!(errori, abs(altezza(t, c) / larghezza_media(t, c) - piatto_rapporto(1.0, 1.0, 1.0)))
        end
        @test errori[3] < errori[2] < errori[1]
        @test errori[3] < 1e-5
    end

    @testset "F2 anello chiuso: altezza esatta" begin
        for n in (10, 12, 16)
            t = anello(fill(CUNEO, n); r_apicale = 15.0, r_basale = 25.0)
            @test rilassa_prova!(t).convergito
            α = π / n
            for c in 1:n
                A = area(t.x, t.cellule[c])
                @test isapprox(semiangolo(t, c), α; rtol = 1e-8)
                @test isapprox(altezza(t, c), cuneo_altezza(Γa, Γb, Γl, A, α); rtol = 1e-8)
            end
        end
    end

    @testset "F3 striscia costretta: sin α = Δ/Γl al centro" begin
        # Δ piccolo, così la striscia resta sotto un giro e non si sovrappone a sé stessa.
        lieve = Intenzioni(1.05, 1.0, 1.0, 100.0)
        α_teor = cuneo_semiangolo(lieve.Γa, lieve.Γb, lieve.Γl)
        errori = Float64[]
        for n in (10, 20)
            t = striscia(fill(lieve, n); larghezza = 7.0, altezza = 14.0)
            @test rilassa_prova!(t).convergito
            @test autointersezioni(t) == 0
            push!(errori, abs(semiangolo(t, n ÷ 2) - α_teor))
        end
        @test errori[2] < errori[1]
        @test errori[2] < 1e-5 * α_teor
    end

    @testset "F4 invarianza di scala e indipendenza da K" begin
        n = 20
        misura(A0, K) = begin
            p = ParametriMeccanici(K = K, η = 1.0, lunghezza_minima = 0.05)
            s = sqrt(A0 / 100)
            t = striscia(fill(Intenzioni(1.05, 1.0, 1.0, A0), n); larghezza = 7.0s, altezza = 14.0s)
            @test rilassa_prova!(t, p).convergito
            c = n ÷ 2
            (α = semiangolo(t, c), h_su_radiceA = altezza(t, c) / sqrt(area(t.x, t.cellule[c])))
        end
        base = misura(100.0, 0.01)
        # lunghezze ×2: l'energia di linea scala ×2, quella d'area come K·L⁴, quindi K/8
        # rende le due strisce esattamente simili
        doppia = misura(400.0, 0.01 / 8)
        più_rigida = misura(100.0, 0.02)
        @test isapprox(doppia.α, base.α; rtol = 1e-6)
        @test isapprox(doppia.h_su_radiceA, base.h_su_radiceA; rtol = 1e-6)
        # cambiare K cambia solo l'effetto dei bordi, piccolo al centro
        @test isapprox(più_rigida.α, base.α; rtol = 1e-4)
    end
end
