# Contatto fra parti lontane del tessuto (motore/05_contatto.jl).
# Regressione: senza contatto una striscia che si arrotola oltre un giro si attraversa.

const P_CONTATTO = ParametriMeccanici(K = 0.01, η = 1.0, lunghezza_minima = 0.05,
    k_contatto = 5.0, d_contatto = 0.5)

arrotolata() = striscia(fill(CUNEO, 20); larghezza = 7.0, altezza = 14.0)

@testset "contatto" begin
    @testset "interruttore: spento si attraversa, acceso no" begin
        spento = arrotolata()
        @test rilassa_prova!(spento).convergito          # P_PROVA ha k_contatto = 0
        @test autointersezioni(spento) > 0
        acceso = arrotolata()
        @test rilassa_prova!(acceso, P_CONTATTO).convergito
        @test autointersezioni(acceso) == 0
        @test energia_contatto(acceso, P_CONTATTO) > 0   # si tocca davvero
        # una configurazione già attraversata è rifiutata appena si accende il contatto
        @test_throws ConfigurazioneNonAmmissibile aggiorna_contatti!(spento, P_CONTATTO)
    end

    t = arrotolata()
    rilassa_prova!(t, P_CONTATTO)
    t.x .+= 0.02 .* randn(Xoshiro(7), length(t.x))
    t.contatti.valida = false

    @testset "le forze di contatto sono il gradiente dell'energia" begin
        g = gradiente!(similar(t.x), t, P_CONTATTO)
        ε = 1e-6
        for k in eachindex(t.x)
            xp = copy(t.x); xp[k] += ε
            xm = copy(t.x); xm[k] -= ε
            gd = (energia(t, P_CONTATTO, xp) - energia(t, P_CONTATTO, xm)) / 2ε
            @test isapprox(g[k], gd; atol = 1e-5, rtol = 1e-5)
        end
    end

    @testset "le forze di contatto sono appaiate" begin
        f = forze(t, P_CONTATTO)
        n = numero_vertici(t)
        scala = maximum(abs, f)
        @test abs(sum(f[2k-1] for k in 1:n)) < 1e-10 * scala * n
        @test abs(sum(f[2k] for k in 1:n)) < 1e-10 * scala * n
        @test abs(sum(t.x[2k-1] * f[2k] - t.x[2k] * f[2k-1] for k in 1:n)) < 1e-10 * scala * n * 100
    end

    @testset "dove nulla si tocca il contatto non cambia niente" begin
        lieve = Intenzioni(1.05, 1.0, 1.0, 100.0)
        a = striscia(fill(lieve, 10); larghezza = 7.0, altezza = 14.0)
        b = striscia(fill(lieve, 10); larghezza = 7.0, altezza = 14.0)
        rilassa_prova!(a)
        rilassa_prova!(b, P_CONTATTO)
        @test energia_contatto(b, P_CONTATTO) == 0
        @test a.x == b.x
    end
end
