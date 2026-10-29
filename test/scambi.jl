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

    @testset "una parete sotto soglia viene scambiata, le altre no" begin
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
    @testset "risultato negativo #14: una giunzione solo apicale non si scambia" begin
        # Due cellule vicine stringono l'anello apicale: la giunzione apicale fra loro si chiude,
        # quella basale no. Lo scambio non deve partire, perché cambierebbe anche la base.
        P = ParametriMeccanici(K = 0.1, η = 1.0, lunghezza_minima = 0.05,
                               k_contatto = 0.0, d_contatto = 0.5)
        ints = fill(PRISMA, 36)
        ints[15] = ints[16] = Intenzioni3D(1.0, 1.0, 1.0, 4.0, 1000.0)
        t = foglio_periodico(ints, 6, 6; lato = 5.0, altezza = 15.0)
        minimizza!(t, P; tolleranza = 1e-9, iterazioni_massime = 100)
        W = poligoni(t)
        p, q = intersect(first.(W[15]), first.(W[16]))
        w = W[15]
        ip, iq = w[findfirst(e -> e[1] == p, w)][2], w[findfirst(e -> e[1] == q, w)][2]
        L(a, b) = Morfogenetica._norma(Morfogenetica._meno(Morfogenetica._pos(t, t.x, b, iq),
                                                           Morfogenetica._pos(t, t.x, a, ip)))
        @test L(2p - 1, 2q - 1) < 0.3          # apicale chiusa
        @test L(2p, 2q) > 5.0                  # basale lunga
        @test isempty(scambi_di_vicine!(t, P; soglia = 0.3, lunghezza_nuova = 1.0))
    end
    @testset "scambi spontanei: un tessuto stirato oltre la soglia scorre" begin
        # Pezzo di 37 prismi, bordo bloccato, tessuto "fluido" (σl = 0.4). Estensione convergente
        # a area costante, x / s e y · s. Fino a s = 2.1 la risposta è elastica, nessuno scambio;
        # a s = 2.15 parte una valanga di scambi e l'energia scende sotto lo stato elastico: il
        # tessuto scorre invece di strapparsi. La valanga arriva al bordo, dove lo scambio non è
        # ancora ammesso (docs/piano.md, sezione 9): lì restano pareti bloccate, tutte di bordo.
        P = ParametriMeccanici(K = 0.1, η = 1.0, lunghezza_minima = 0.05,
                               k_contatto = 0.0, d_contatto = 0.5)
        stirato(s) = begin
            t = foglio_esagonale(fill(Intenzioni3D(1.0, 1.0, 0.4, 0.0, 1000.0), 37), 3;
                                 lato = 5.0, altezza = 15.0, bordo_fisso = true)
            for k in 1:length(t.x)÷3
                t.x[3k-2] /= s; t.x[3k-1] *= s
            end
            t.ancore .= t.x
            t
        end
        # sotto la soglia: elastico
        t = stirato(2.1)
        r = rilassa_con_scambi!(t, P; soglia = 0.3, lunghezza_nuova = 1.0, tolleranza = 1e-9)
        @test r.convergito && isempty(r.registro)
        # sopra la soglia: scorre
        elastico = stirato(2.15)
        @test minimizza!(elastico, P; tolleranza = 1e-9, iterazioni_massime = 300).convergito
        t = stirato(2.15)
        r = rilassa_con_scambi!(t, P; soglia = 0.3, lunghezza_nuova = 1.0, tolleranza = 1e-9,
                                cicli_massimi = 5)
        @test length(r.registro) >= 10
        @test energia(t, P) < energia(elastico, P) - 100
        W = poligoni(t)
        interni = [j for j in 1:length(t.x)÷6 if !t.fissi[2j - 1]]
        @test all(length(cellule_attorno(W, j)) == 3 for j in interni)
        @test all(volume_cellula(t, c) > 0 for c in eachindex(t.cellule))
        # quello che resta bloccato tocca il bordo
        for w in W, k in eachindex(w)
            (j1, i1), (j2, i2) = w[k], w[mod1(k + 1, length(w))]
            L = Morfogenetica._norma(Morfogenetica._meno(Morfogenetica._pos(t, t.x, 2j2 - 1, i2),
                                                         Morfogenetica._pos(t, t.x, 2j1 - 1, i1)))
            L < 0.3 && @test t.fissi[2j1 - 1] || t.fissi[2j2 - 1]
        end
    end
end
