# Il guscio 3D (motore/07_guscio.jl): geometria, forze, F6.

const P_3D = ParametriMeccanici(K = 1e-3, η = 1.0, lunghezza_minima = 0.05,
    k_contatto = 0.0, d_contatto = 0.5)
const PRISMA = Intenzioni3D(1.0, 1.0, 1.0, 0.0, 1000.0)

lato_esagono(t, c) = sqrt(area_apicale(t, c) / ESAGONO)
minimizza_3d!(t, p = P_3D) = minimizza!(t, p; tolleranza = 1e-9, iterazioni_massime = 500)

# Passo delle differenze finite 1e-4: con energie di migliaia di nN·µm un passo più piccolo è
# dominato dall'arrotondamento (l'errore cresce come 1/ε), uno più grande dal troncamento.
function errore_gradiente(t, p)
    g = gradiente!(similar(t.x), t, p)
    ε = 1e-4
    e = 0.0
    for k in eachindex(t.x)
        xp = copy(t.x); xp[k] += ε
        xm = copy(t.x); xm[k] -= ε
        e = max(e, abs(g[k] - (energia(t, p, xp) - energia(t, p, xm)) / 2ε) / max(1, abs(g[k])))
    end
    return e
end

@testset "guscio 3D" begin
    @testset "geometria: prismi regolari esatti" begin
        for t in (foglio_esagonale(fill(PRISMA, 19), 2; lato = 5.0, altezza = 15.0),
                  foglio_periodico(fill(PRISMA, 16), 4, 4; lato = 5.0, altezza = 15.0))
            for c in eachindex(t.cellule)
                @test isapprox(volume_cellula(t, c), ESAGONO * 25 * 15; rtol = 1e-12)
                @test isapprox(area_apicale(t, c), ESAGONO * 25; rtol = 1e-12)
                @test isapprox(altezza_prisma(t, c), 15.0; rtol = 1e-12)
            end
        end
    end

    @testset "le forze sono il gradiente dell'energia" begin
        varia = Intenzioni3D(1.2, 0.9, 1.1, 0.4, 1000.0)
        a = foglio_esagonale(fill(varia, 7), 1; lato = 5.0, altezza = 15.0)
        a.x .+= 0.3 .* randn(Xoshiro(11), length(a.x))
        @test errore_gradiente(a, P_3D) < 1e-7
        # periodico: anche le derivate rispetto alle due dimensioni del periodo
        b = foglio_periodico(fill(varia, 16), 4, 4; lato = 5.0, altezza = 15.0)
        b.x[1:end-2] .+= 0.3 .* randn(Xoshiro(12), length(b.x) - 2)
        b.x[end-1] *= 1.03
        @test errore_gradiente(b, P_3D) < 1e-7
    end

    @testset "le forze sono appaiate: risultante e momento nulli" begin
        t = foglio_esagonale(fill(Intenzioni3D(1.2, 0.9, 1.1, 0.4, 1000.0), 7), 1;
                             lato = 5.0, altezza = 15.0)
        t.x .+= 0.3 .* randn(Xoshiro(13), length(t.x))
        f = forze(t, P_3D)
        n = length(t.x) ÷ 3
        scala = maximum(abs, f) * n
        for d in 1:3
            @test abs(sum(f[3k-3+d] for k in 1:n)) < 1e-10 * scala
        end
        r(k) = (t.x[3k-2], t.x[3k-1], t.x[3k])
        F(k) = (f[3k-2], f[3k-1], f[3k])
        M = reduce((a, b) -> a .+ b, (Morfogenetica._croce(r(k), F(k)) for k in 1:n))
        @test maximum(abs, M) < 1e-10 * scala * 30
    end

    @testset "bordo bloccato: i vertici fissi non sono incognite" begin
        t = foglio_esagonale(fill(Intenzioni3D(1.2, 0.9, 1.1, 0.4, 1000.0), 19), 2;
                             lato = 5.0, altezza = 15.0, bordo_fisso = true)
        @test count(t.fissi) > 0
        t.x .+= 0.3 .* randn(Xoshiro(14), length(t.x))
        t.ancore .= t.x
        @test errore_gradiente(t, P_3D) < 1e-7
        g = gradiente!(similar(t.x), t, P_3D)
        @test all(g[repeat(t.fissi, inner = 3)] .== 0)
        # cambiare la voce di x di un vertice fisso non cambia l'energia
        j = findfirst(t.fissi)
        y = copy(t.x); y[3j] += 1.0
        @test energia(t, P_3D, y) == energia(t, P_3D)
    end

    @testset "F6a prisma isolato: lato e rapporto esatti" begin
        # Solo con σa = σb e σl ≥ σa. Un prisma isolato non è costretto ad avere apice e base
        # uguali: con σa ≠ σb diventa un tronco di piramide (nel foglio è il reticolo a imporli
        # uguali). E con σl < σa il prisma regolare non è un minimo (prova sotto).
        for (σa, σb, σl) in ((1.0, 1.0, 1.0), (1.0, 1.0, 1.1), (1.2, 1.2, 1.2))
            t = foglio_esagonale([Intenzioni3D(σa, σb, σl, 0.0, 1000.0)], 0; lato = 5.0, altezza = 15.0)
            @test minimizza_3d!(t).convergito
            V = volume_cellula(t, 1)
            @test isapprox(lato_esagono(t, 1), prisma_isolato_lato(σa, σb, σl, V); rtol = 1e-9)
            @test isapprox(altezza_prisma(t, 1) / lato_esagono(t, 1),
                           prisma_foglio_rapporto(σa, σb, σl) / 2; rtol = 1e-9)
        end
    end

    @testset "risultato: con σl < σa il prisma isolato si arrotonda" begin
        # Se i fianchi costano meno di apice e base, la cellula isolata lascia la forma di
        # prisma: tre spigoli apicali su sei collassano (fino alla lunghezza minima, dove
        # servirebbe uno scambio topologico) e l'energia scende sotto quella del miglior prisma
        # regolare, 458.95 nN·µm contro 456.03 per σa = σb = 1, σl = 0.9.
        t = foglio_esagonale([Intenzioni3D(1.0, 1.0, 0.9, 0.0, 1000.0)], 0; lato = 5.0, altezza = 15.0)
        minimizza_3d!(t)
        P = Morfogenetica._vertici(t, t.x, t.facce[1])
        lati = sort([Morfogenetica._norma(Morfogenetica._meno(P[mod1(k + 1, 6)], P[k])) for k in 1:6])
        @test all(lati[1:3] .< 2 * P_3D.lunghezza_minima)
        @test energia(t, P_3D) < 458.9
    end

    @testset "F6b foglio periodico: altezza/lato = √3·(σa+σb)/σl" begin
        for (σa, σb, σl) in ((1.0, 1.0, 1.0), (1.0, 0.5, 1.5), (1.4, 0.8, 1.2))
            t = foglio_periodico(fill(Intenzioni3D(σa, σb, σl, 0.0, 1000.0), 16), 4, 4;
                                 lato = 5.0, altezza = 15.0)
            @test minimizza_3d!(t).convergito
            for c in eachindex(t.cellule)
                @test isapprox(altezza_prisma(t, c) / lato_esagono(t, c),
                               prisma_foglio_rapporto(σa, σb, σl); rtol = 1e-9)
            end
        end
    end

    @testset "F6c invarianza di scala e indipendenza da K" begin
        misura(s, K) = begin
            p = ParametriMeccanici(K = K, η = 1.0, lunghezza_minima = 0.05,
                                   k_contatto = 0.0, d_contatto = 0.5)
            t = foglio_periodico(fill(Intenzioni3D(1.2, 0.8, 1.1, 0.0, 1000.0 * s^3), 16), 4, 4;
                                 lato = 5.0s, altezza = 15.0s)
            @test minimizza_3d!(t, p).convergito
            altezza_prisma(t, 1) / lato_esagono(t, 1)
        end
        base = misura(1.0, 1e-3)
        # lunghezze ×2: le tensioni scalano ×4, il volume come K·L⁶, quindi K/16
        @test isapprox(misura(2.0, 1e-3 / 16), base; rtol = 1e-9)
        @test isapprox(misura(1.0, 5e-3), base; rtol = 1e-9)
    end

    @testset "una cellula che stringe l'anello apicale si restringe all'apice" begin
        ints = fill(PRISMA, 36)
        centro = 3 * 6 + 4
        ints[centro] = Intenzioni3D(1.0, 1.0, 1.0, 2.0, 1000.0)
        t = foglio_periodico(ints, 6, 6; lato = 5.0, altezza = 15.0)
        @test minimizza_3d!(t).convergito
        @test area_apicale(t, centro) < 0.4 * area_apicale(t, 1)
    end

    @testset "risultato negativo: un pezzo di foglio con bordo libero non ha equilibrio" begin
        # Al bordo le facce laterali sono di una cellula sola: la loro tensione accorcia gli
        # spigoli del bordo fino alla lunghezza minima, e senza scambi di vicine non c'è un
        # equilibrio ammissibile. Per questo F6b usa il foglio periodico.
        t = foglio_esagonale(fill(PRISMA, 7), 1; lato = 5.0, altezza = 15.0)
        r = minimizza_3d!(t)
        @test !r.convergito
        corto = minimum(f -> begin
                P = Morfogenetica._vertici(t, t.x, f)
                minimum(Morfogenetica._norma(Morfogenetica._meno(P[mod1(k + 1, length(P))], P[k]))
                        for k in eachindex(P))
            end, t.facce)
        @test corto < 2 * P_3D.lunghezza_minima
    end
end
