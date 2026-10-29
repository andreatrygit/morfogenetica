# Scambi di vicine (motore/08_scambi.jl).

const P_SCAMBI = P_3D
foglio_scambi() = foglio_periodico(fill(PRISMA, 16), 4, 4; lato = 5.0, altezza = 15.0)
# uno spigolo interno qualunque: il primo lato della cellula 6
spigolo_di_prova(t) = (P = poligoni(t); (P[6][1][1], P[6][2][1]))

# Avvicina i due vertici dello spigolo p–q (apicali e basali) a distanza ε, sul loro asse.
function accorcia!(t, p, q, ε)
    P = poligoni(t)
    c = only(intersect(cellule_attorno(P, p), cellule_attorno(P, q))[1:1])
    ip = P[c][findfirst(e -> e[1] == p, P[c])][2]
    iq = P[c][findfirst(e -> e[1] == q, P[c])][2]
    Lx, Ly = Morfogenetica._periodo(t, t.x)
    for (va, vb) in ((2p - 1, 2q - 1), (2p, 2q))
        a = Morfogenetica._pos(t, t.x, va, ip); b = Morfogenetica._pos(t, t.x, vb, iq)
        m = (a .+ b) ./ 2
        e = Morfogenetica._meno(b, a) ./ Morfogenetica._norma(Morfogenetica._meno(b, a))
        na = m .- e .* (ε / 2); nb = m .+ e .* (ε / 2)
        t.x[3va-2:3va] .= (na[1] - ip[1] * Lx, na[2] - ip[2] * Ly, na[3])
        t.x[3vb-2:3vb] .= (nb[1] - iq[1] * Lx, nb[2] - iq[2] * Ly, nb[3])
    end
    return t
end

insiemi(t) = sort([sort(first.(w)) for w in poligoni(t)])
volume_totale(t) = sum(volume_cellula(t, c) for c in eachindex(t.cellule))

@testset "scambi di vicine" begin
    @testset "topologia: due cellule perdono un lato, due lo guadagnano" begin
        t = foglio_scambi()
        p, q = spigolo_di_prova(t)
        lati = length.(poligoni(t))
        r = scambia!(t, p, q; lunghezza_nuova = 5.0)
        @test r isa ScambioDiVicine
        dopo = length.(poligoni(t))
        A, B = r.perdono; C, D = r.guadagnano
        @test dopo[A] == lati[A] - 1 && dopo[B] == lati[B] - 1
        @test dopo[C] == lati[C] + 1 && dopo[D] == lati[D] + 1
        @test sum(dopo) == sum(lati)
        @test length(t.facce) == length(foglio_scambi().facce)
        P = poligoni(t)
        @test all(length(cellule_attorno(P, j)) == 3 for j in 1:(length(t.x) - 2) ÷ 6)
        @test all(volume_cellula(t, c) > 0 for c in eachindex(t.cellule))
        @test errore_gradiente(t, P_SCAMBI) < 1e-7
    end

    @testset "lo scambio è reversibile" begin
        t = foglio_scambi()
        prima = insiemi(t)
        p, q = spigolo_di_prova(t)
        scambia!(t, p, q; lunghezza_nuova = 5.0)
        @test insiemi(t) != prima
        scambia!(t, p, q; lunghezza_nuova = 5.0)
        # stessa rete di cellule, con i due vertici dello spigolo che si sono scambiati il nome
        rinomina(j) = j == p ? q : j == q ? p : j
        @test insiemi(t) == prima || sort([sort(rinomina.(w)) for w in insiemi(t)]) == prima
    end

    @testset "l'energia è continua attraverso lo scambio" begin
        salto(ε) = begin
            t = accorcia!(foglio_scambi(), spigolo_di_prova(foglio_scambi())..., ε)
            E = energia(t, P_SCAMBI)
            V = volume_totale(t)
            scambia!(t, spigolo_di_prova(foglio_scambi())...; lunghezza_nuova = ε)
            @test isapprox(volume_totale(t), V; rtol = 1e-12)
            abs(energia(t, P_SCAMBI) - E)
        end
        s2, s3 = salto(1e-2), salto(1e-3)
        @test 9 < s2 / s3 < 11          # il salto va a zero come ε
        @test s3 < 0.1
    end

    @testset "uno spigolo sotto soglia viene scambiato, gli altri no" begin
        t = foglio_scambi()
        p, q = spigolo_di_prova(t)
        accorcia!(t, p, q, 0.1)
        r = scambi_di_vicine!(t, P_SCAMBI; soglia = 0.3, lunghezza_nuova = 1.0)
        @test length(r) == 1
        @test minmax(r[1].p, r[1].q) == minmax(p, q)
    end

    @testset "al bordo bloccato lo scambio non è ammesso" begin
        t = foglio_esagonale(fill(PRISMA, 19), 2; lato = 5.0, altezza = 15.0, bordo_fisso = true)
        P = poligoni(t)
        bordo = findfirst(w -> any(e -> t.fissi[2e[1] - 1], w), P)
        w = P[bordo]
        k = findfirst(i -> t.fissi[2w[i][1] - 1], eachindex(w))
        @test scambia!(t, w[k][1], w[mod1(k + 1, length(w))][1]; lunghezza_nuova = 1.0) === nothing
    end

    @testset "un difetto da scambio forzato è metastabile" begin
        # Dopo uno scambio forzato nel reticolo regolare, le due cellule a 5 lati e le due a 7
        # restano: il tessuto si ferma in un minimo locale più alto, senza tornare indietro.
        # Tornare indietro con uno scambio esplicito riporta al reticolo regolare.
        regolare = foglio_scambi()
        @test minimizza_3d!(regolare).convergito
        Ereg = energia(regolare, P_SCAMBI)
        t = foglio_scambi()
        p, q = spigolo_di_prova(t)
        scambia!(t, p, q; lunghezza_nuova = 1.0)
        r = rilassa_con_scambi!(t, P_SCAMBI; soglia = 0.3, lunghezza_nuova = 1.0, tolleranza = 1e-9)
        @test r.convergito
        @test isempty(r.registro)
        @test sort(unique(length.(poligoni(t)))) == [5, 6, 7]
        @test energia(t, P_SCAMBI) > Ereg + 10
        scambia!(t, p, q; lunghezza_nuova = 1.0)
        @test minimizza_3d!(t).convergito
        @test isapprox(energia(t, P_SCAMBI), Ereg; rtol = 1e-9)
    end
end
