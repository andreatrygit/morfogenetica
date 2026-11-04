# R6 (docs/roadmap.md): il solutore sparso, per portare all'equilibrio tessuti di migliaia di
# cellule (la catena di gastrulazione a scala vera, circa 1 mm).
#
#   minimizza_sparso!(t, p; tolleranza, iterazioni_massime, spostamento_massimo, ancora, peso)
#       lo stesso metodo di `minimizza!` (06_integratore.jl: Newton smorzato, passo accettato solo
#       se l'energia scende), con la matrice delle derivate seconde sparsa: due vertici sono
#       accoppiati solo se stanno in una stessa cellula, e la matrice si calcola per differenze del
#       gradiente esatto perturbando insieme i vertici che non hanno vicini in comune (colorazione
#       a distanza due). Le due voci della scatola dei fogli periodici sono righe e colonne piene.
#   Morfogenetica._hessiana_sparsa(t, p, x)      la matrice, per le prove
#   Embrione(...; solutore = :pieno o :sparso)  il solutore della passata (predefinito :pieno)
#
# Ogni prova è una funzione che restituisce un Bool. Scritte prima del codice
# (motore/15_sparso.jl; commit rosso, con @test_broken), passano dal commit verde.

# La matrice delle derivate seconde piena, per differenze del gradiente esatto (come in 06).
function hessiana_piena(t, p, x)
    n = length(x); H = zeros(n, n); gp = similar(x); gm = similar(x); ε = 1e-6
    y = copy(x)
    for j in 1:n
        y[j] = x[j] + ε; gradiente!(gp, t, p, y)
        y[j] = x[j] - ε; gradiente!(gm, t, p, y)
        y[j] = x[j]
        H[:, j] .= (gp .- gm) ./ 2ε
    end
    return H
end

function tessuto_disturbato(periodico)
    t = periodico ? foglio_periodico(fill(PRISMA, 16), 4, 4; lato = 5.0, altezza = 15.0) :
                    foglio_esagonale(fill(PRISMA, 19), 2; lato = 5.0, altezza = 15.0, bordo_fisso = true)
    rng = MersenneTwister(4)
    n = length(t.x) - (periodico ? 2 : 0)
    t.x[1:n] .+= 0.3 .* (rand(rng, n) .- 0.5)
    for c in 1:3
        t.intenzioni[c] = Intenzioni3D(1.0, 1.0, 1.0, 2.0, 900.0)
    end
    return t
end

function prova_hessiana_sparsa_uguale_alla_piena()
    for periodico in (false, true)
        t = tessuto_disturbato(periodico)
        H = hessiana_piena(t, P_PROGRAMMA, t.x)
        S = Morfogenetica._hessiana_sparsa(t, P_PROGRAMMA, t.x)
        maximum(abs, Matrix(S) .- H) <= 1e-6 * max(1, maximum(abs, H)) || return false
    end
    return true
end

function prova_stesso_minimo()
    for periodico in (false, true), ancorato in (false, true)
        a = tessuto_disturbato(periodico); b = deepcopy(a)
        kw = ancorato ? (ancora = copy(a.x), peso = 1.0) : NamedTuple()
        ra = minimizza!(a, P_PROGRAMMA; tolleranza = 1e-9, iterazioni_massime = 300, kw...)
        rb = minimizza_sparso!(b, P_PROGRAMMA; tolleranza = 1e-9, iterazioni_massime = 300, kw...)
        (ra.convergito && rb.convergito) || return false
        # un foglio periodico senza ancora si può traslare rigidamente senza cambiare energia: i
        # due minimi si confrontano a meno della traslazione media dei vertici (corretto il 4
        # novembre: prima si confrontavano le posizioni assolute, e passava per caso)
        nv = Morfogenetica._nvert(a)
        d = [a.x[3k-3+c] - b.x[3k-3+c] for k in 1:nv, c in 1:3]
        periodico && !ancorato && (d .-= sum(d; dims = 1) ./ nv)
        maximum(abs, d) < 1e-6 || return false
        maximum(abs, a.x[3nv+1:end] .- b.x[3nv+1:end]; init = 0.0) < 1e-6 || return false
        isapprox(energia(a, P_PROGRAMMA), energia(b, P_PROGRAMMA); rtol = 1e-10) || return false
    end
    return true
end

function prova_passata_sparsa_come_piena()
    a = caso_passata(); b = caso_passata()
    b.solutore = :sparso
    for _ in 1:3
        passata!(a); passata!(b)
    end
    return maximum(abs, a.tessuto.x .- b.tessuto.x) < 1e-6
end

function prova_scala_vera()
    # il foglio dell'organizzatore (1672 cellule, circa 1 × 1 mm) con una regione larga 300 µm che
    # costringe: il solutore sparso lo porta all'equilibrio. Il foglio parte all'equilibrio: con
    # σa = σb = 1 un foglio di prismi sta in piedi con altezza/lato = √3·(σa + σb)/σl, quindi per
    # lato e altezza di 15 µm serve σl = 2√3 (tensioni libere). Corretto nel verde: con σl = 1 il
    # foglio voleva cambiare proporzioni, la scatola sentiva migliaia di nN e, a 1 µm per
    # iterazione, non arrivava in 100 iterazioni.
    σl = 2sqrt(3)
    t = foglio_periodico(fill(Intenzioni3D(1.0, 1.0, σl, 0.0, 8768.0), 38 * 44), 38, 44;
                         lato = 15.0, altezza = 15.0)
    Lx, Ly = t.x[end-1], t.x[end]
    for c in eachindex(t.cellule)
        g = centro_cellula(t, c)
        hypot(g[1] - Lx / 2, g[2] - Ly / 2) < 150 && (t.intenzioni[c] = Intenzioni3D(1.0, 1.0, σl, 8.0, 8768.0))
    end
    r = minimizza_sparso!(t, P_PROGRAMMA; tolleranza = 1e-6, iterazioni_massime = 100)
    return r.convergito
end

# Newton modificato (4 novembre): costruire la matrice costa 2.7 s sul foglio grande, fattorizzarla
# 0.3 s, e un passo meccanico di 5 minuti con una regione che costringe faceva 5–10 iterazioni,
# ognuna con una matrice nuova (12 s). Riusando la matrice finché la forza scende abbastanza il
# passo deve costare meno di 5 s, arrivando allo stesso equilibrio (le prove sopra).
const _FOGLIO_GRANDE = Ref{Any}(nothing)
function foglio_grande_rilassato()
    if _FOGLIO_GRANDE[] === nothing
        σl = 2sqrt(3)
        t = foglio_periodico(fill(Intenzioni3D(1.0, 1.0, σl, 0.0, 8768.0), 38 * 44), 38, 44;
                             lato = 15.0, altezza = 15.0)
        minimizza_sparso!(t, P_PROGRAMMA; tolleranza = 1e-5, iterazioni_massime = 100)
        _FOGLIO_GRANDE[] = t
    end
    return deepcopy(_FOGLIO_GRANDE[])
end

function prova_passo_grande_veloce()
    t = foglio_grande_rilassato()
    Lx, Ly = t.x[end-1], t.x[end]
    for c in eachindex(t.cellule)
        g = centro_cellula(t, c)
        hypot(g[1] - Lx / 2, g[2] - Ly / 2) < 75 &&
            (t.intenzioni[c] = Intenzioni3D(1.0, 1.0, 2sqrt(3), 8.0, 8768.0))
    end
    passo() = minimizza_sparso!(t, P_PROGRAMMA; tolleranza = 1e-5, iterazioni_massime = 200,
                                ancora = copy(t.x), peso = P_PROGRAMMA.η / 5.0)
    passo(); passo()                                   # i primi due passi, i più lunghi
    tempo = @elapsed r = passo()
    return r.convergito && tempo < 5.0
end

@testset "R6 solutore sparso" begin
    @test prova_hessiana_sparsa_uguale_alla_piena()
    @test prova_stesso_minimo()
    @test prova_passata_sparsa_come_piena()
    @lenta @test prova_scala_vera()
    @lenta @test_broken prova_passo_grande_veloce()
end
