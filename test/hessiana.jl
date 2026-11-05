# La matrice delle derivate seconde esatta (R6, 5 novembre): `Morfogenetica._hessiana_esatta(t, p,
# x; ancora, peso)` in motore/16_hessiana.jl, per differenziazione automatica in avanti
# (ForwardDiff) dei termini locali dell'energia: faccia per faccia (tensione per area e anello
# apicale), cellula per cellula (volume: K·(∇V ∇Vᵀ + (V − V0) ∇²V)). Usata da
# `minimizza_sparso!(...; hessiana = :esatta)` e da `Embrione(...; solutore = :esatto)`.
#
# Ogni prova è una funzione che restituisce un Bool. Commit rosso: marcate @test_broken finché
# il codice non c'è.

function prova_esatta_come_differenze()
    P = P_PROGRAMMA
    for t in tessuti_per_impronta(), ancorato in (false, true)
        kw = ancorato ? (ancora = t.x .+ 0.1, peso = 0.7) : NamedTuple()
        E = Matrix(Morfogenetica._hessiana_esatta(t, P, t.x; kw...))
        D = hessiana_piena(t, P, t.x)
        ancorato && for k in axes(D, 1); D[k, k] += 0.7; end
        maximum(abs, E .- D) <= 1e-5 * max(1, maximum(abs, D)) || return false
        E == E' || return false
    end
    return true
end

function prova_stesso_minimo_esatto()
    for periodico in (false, true), ancorato in (false, true)
        a = tessuto_disturbato(periodico); b = deepcopy(a)
        kw = ancorato ? (ancora = copy(a.x), peso = 1.0) : NamedTuple()
        ra = minimizza!(a, P_PROGRAMMA; tolleranza = 1e-9, iterazioni_massime = 300, kw...)
        rb = minimizza_sparso!(b, P_PROGRAMMA; tolleranza = 1e-9, iterazioni_massime = 300,
                               hessiana = :esatta, kw...)
        (ra.convergito && rb.convergito) || return false
        nv = Morfogenetica._nvert(a)
        d = [a.x[3k-3+c] - b.x[3k-3+c] for k in 1:nv, c in 1:3]
        periodico && !ancorato && (d .-= sum(d; dims = 1) ./ nv)
        maximum(abs, d) < 1e-6 || return false
        isapprox(energia(a, P_PROGRAMMA), energia(b, P_PROGRAMMA); rtol = 1e-10) || return false
    end
    return true
end

function prova_passata_esatta_come_sparsa()
    # costruita con la parola chiave: un solutore che non esiste è rifiutato dal costruttore
    con(solutore) = (e = caso_passata();
                     Embrione(e.tessuto, invaginazione(; Λ = 3.0, soglia = 0.001, tasso = 20.0);
                              parametri = P_PROGRAMMA, dt = 1.0, campi = e.campi,
                              stati = [c.stato for c in e.cellule], solutore,
                              interruttori = Interruttori(; topologia = false)))
    a = con(:sparso); b = con(:esatto)
    for _ in 1:3
        passata!(a); passata!(b)
    end
    return maximum(abs, a.tessuto.x .- b.tessuto.x) < 1e-6
end

function prova_esatta_piu_veloce()
    # sul foglio di 1672 cellule, nello stesso processo: almeno due volte più veloce della matrice
    # per differenze del gradiente
    t = foglio_grande_rilassato()
    s = Morfogenetica._struttura(t)
    Morfogenetica._hessiana_sparsa(t, P_PROGRAMMA, t.x; struttura = s)
    Morfogenetica._hessiana_esatta(t, P_PROGRAMMA, t.x)
    td = minimum(@elapsed(Morfogenetica._hessiana_sparsa(t, P_PROGRAMMA, t.x; struttura = s)) for _ in 1:3)
    te = minimum(@elapsed(Morfogenetica._hessiana_esatta(t, P_PROGRAMMA, t.x)) for _ in 1:3)
    return te < td / 2
end

@testset "R6 la matrice delle derivate seconde esatta" begin
    @test_broken prova_esatta_come_differenze()
    @test_broken prova_stesso_minimo_esatto()
    @test_broken prova_passata_esatta_come_sparsa()
    @lenta @test_broken prova_esatta_piu_veloce()
end
