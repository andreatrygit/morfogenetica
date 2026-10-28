# Prove di meccanica pura: non dipendono dal modello biologico.

const P_PROVA = ParametriMeccanici(K = 0.01, η = 1.0, lunghezza_minima = 0.05)
const CUNEO = Intenzioni(1.3, 1.0, 1.0, 100.0)

function tessuto_disturbato(seme)
    rng = Xoshiro(seme)
    t = anello(fill(CUNEO, 12); r_apicale = 15.0, r_basale = 25.0)
    t.x .+= 0.5 .* randn(rng, length(t.x))
    return t
end

@testset "meccanica" begin
    @testset "le forze sono il gradiente dell'energia" begin
        for seme in 1:3
            t = tessuto_disturbato(seme)
            g = gradiente!(similar(t.x), t, P_PROVA)
            ε = 1e-6
            for k in eachindex(t.x)
                xp = copy(t.x); xp[k] += ε
                xm = copy(t.x); xm[k] -= ε
                gd = (energia(t, P_PROVA, xp) - energia(t, P_PROVA, xm)) / 2ε
                @test isapprox(g[k], gd; atol = 1e-6, rtol = 1e-6)
            end
        end
    end

    @testset "le forze sono appaiate: risultante e momento nulli" begin
        t = tessuto_disturbato(4)
        f = forze(t, P_PROVA)
        n = numero_vertici(t)
        Fx = sum(f[2k-1] for k in 1:n)
        Fy = sum(f[2k] for k in 1:n)
        M = sum(t.x[2k-1] * f[2k] - t.x[2k] * f[2k-1] for k in 1:n)
        scala = maximum(abs, f)
        @test abs(Fx) < 1e-10 * scala * n
        @test abs(Fy) < 1e-10 * scala * n
        @test abs(M) < 1e-10 * scala * n * 30
    end

    @testset "l'energia scende a ogni passo" begin
        t = tessuto_disturbato(5)
        E = energia(t, P_PROVA)
        sale = 0
        for _ in 1:2000
            passo!(t, P_PROVA, 0.2)
            En = energia(t, P_PROVA)
            En > E + 1e-12 * abs(E) && (sale += 1)
            E = En
        end
        @test sale == 0
    end

    @testset "dimezzare il passo dimezza l'errore (Eulero, primo ordine)" begin
        # posizione al tempo T con passi dt, dt/2, dt/4
        function traiettoria(dt)
            t = tessuto_disturbato(6)
            for _ in 1:round(Int, 20 / dt)
                passo!(t, P_PROVA, dt)
            end
            return t.x
        end
        x1, x2, x4 = traiettoria(0.2), traiettoria(0.1), traiettoria(0.05)
        e1 = maximum(abs, x1 .- x4)
        e2 = maximum(abs, x2 .- x4)
        @test 2.5 < e1 / e2 < 3.5   # con riferimento dt/4 il rapporto atteso è 3
    end

    @testset "una configurazione non ammissibile è rifiutata" begin
        # Anello con α = π/N imposto: il lato apicale si annulla se sin(π/N) > Γl/(2Γa)
        # (docs/piano.md). Con Γa = 1.3 e Γl = 1 la soglia cade fra N = 7 e N = 9.
        for N in (6, 7)
            @test sin(π / N) > CUNEO.Γl / (2CUNEO.Γa)
            t = anello(fill(CUNEO, N); r_apicale = 15.0, r_basale = 25.0)
            @test_throws ConfigurazioneNonAmmissibile rilassa!(t, P_PROVA; dt = 0.2,
                tolleranza = 1e-9, passi_massimi = 2_000_000)
        end
        t = anello(fill(CUNEO, 9); r_apicale = 15.0, r_basale = 25.0)
        @test sin(π / 9) < CUNEO.Γl / (2CUNEO.Γa)
        @test rilassa!(t, P_PROVA; dt = 0.2, tolleranza = 1e-9, passi_massimi = 2_000_000).convergito
    end
end
