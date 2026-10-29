# Scutoidi (motore/10_scutoidi.jl): lo scambio apicale con scutoide.

foglio_scutoidi() = foglio_esagonale(fill(PRISMA, 37), 3; lato = 5.0, altezza = 15.0,
                                     bordo_fisso = true)
const P_SCUTOIDI = ParametriMeccanici(K = 0.1, η = 1.0, lunghezza_minima = 0.05,
    k_contatto = 0.0, d_contatto = 0.5)

function spigolo_apicale_12(t)
    ap(c) = t.facce[Morfogenetica._faccia_apicale_indice(t, c)].vertici
    return Tuple(intersect(ap(1), ap(2)))
end

# Avvicina i due vertici apicali a distanza ε, sul loro asse; la base non si muove.
function accorcia_apicale!(t, a1, a2, ε)
    pa, pb = Morfogenetica._pos(t, t.x, a1, (0, 0)), Morfogenetica._pos(t, t.x, a2, (0, 0))
    m = (pa .+ pb) ./ 2
    e = Morfogenetica._meno(pb, pa) ./ Morfogenetica._norma(Morfogenetica._meno(pb, pa))
    t.x[3a1-2:3a1] .= m .- e .* (ε / 2)
    t.x[3a2-2:3a2] .= m .+ e .* (ε / 2)
    return t
end

@testset "scutoidi" begin
    @testset "lo scambio apicale crea uno scutoide" begin
        t = foglio_scutoidi()
        a1, a2 = spigolo_apicale_12(t)
        nv = length(t.x) ÷ 3
        r = scambio_apicale!(t, a1, a2; lunghezza_nuova = 1.0)
        @test r isa ScambioApicale
        @test r.perdono == (1, 2) || r.perdono == (2, 1)
        @test length(t.x) ÷ 3 == nv + 1                      # il vertice interno
        # il vertice interno sta nel foglio, sotto l'apice
        @test 0 < t.x[3r.v] < 15
        # ogni cellula resta un poliedro chiuso, con volume positivo
        @test all(caratteristica_di_eulero(t, c) == 2 for c in eachindex(t.cellule))
        @test all(volume_cellula(t, c) > 0 for c in eachindex(t.cellule))
        # in alto A e B non si toccano più, C e D sì; la parete fra A e B c'è ancora
        A, B = r.perdono; C, D = r.guadagnano
        ap(c) = t.facce[Morfogenetica._faccia_apicale_indice(t, c)].vertici
        @test isempty(intersect(ap(A), ap(B)))
        @test length(intersect(ap(C), ap(D))) == 2
        @test any(f -> f.etichetta == laterale && A in f.cellule && B in f.cellule, t.facce)
        @test errore_gradiente(t, P_SCUTOIDI) < 1e-7
    end

    @testset "l'energia è continua e il volume si conserva" begin
        salto(ε) = begin
            t = accorcia_apicale!(foglio_scutoidi(), spigolo_apicale_12(foglio_scutoidi())..., ε)
            E = energia(t, P_SCUTOIDI)
            V = sum(volume_cellula(t, c) for c in eachindex(t.cellule))
            scambio_apicale!(t, spigolo_apicale_12(foglio_scutoidi())...;
                             lunghezza_nuova = ε, profondita = ε)
            @test isapprox(sum(volume_cellula(t, c) for c in eachindex(t.cellule)), V; rtol = 1e-12)
            abs(energia(t, P_SCUTOIDI) - E)
        end
        s2, s3 = salto(1e-2), salto(1e-3)
        @test 9 < s2 / s3 < 11
        @test s3 < 1.0
    end

    @testset "al bordo bloccato e nei tessuti periodici non è ammesso" begin
        t = foglio_scutoidi()
        ap = t.facce[Morfogenetica._faccia_apicale_indice(t, 37)].vertici
        k = findfirst(i -> t.fissi[ap[i]], eachindex(ap))
        @test scambio_apicale!(t, ap[k], ap[mod1(k + 1, length(ap))]; lunghezza_nuova = 1.0) === nothing
        q = foglio_periodico(fill(PRISMA, 16), 4, 4; lato = 5.0, altezza = 15.0)
        w = q.facce[Morfogenetica._faccia_apicale_indice(q, 6)].vertici
        @test scambio_apicale!(q, w[1], w[2]; lunghezza_nuova = 1.0) === nothing
    end
end
