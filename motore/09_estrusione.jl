# Estrusione: rimozione di una cellula collassata (specifica, decisione 4, e repertorio:
# "morire: apoptosi, oppure estrusione dal foglio").
#
# Una cellula triangolare la cui faccia apicale e quella basale sono scese tutte e due sotto una
# soglia d'area esce dal foglio: i suoi tre vertici apicali si fondono in uno, così i tre basali,
# e le tre vicine si chiudono su quel vertice (perdono un lato ciascuna). È l'operazione T2 dei
# modelli a vertici, fatta sulle due superfici del guscio. La cellula arriva a essere un
# triangolo da sé, con gli scambi di vicine sulle sue pareti che si accorciano.
#
# Si estrude solo una cellula piccola tutta intera: nel guscio a prismi apice e base cambiano
# insieme, e togliere una cellula che ha chiuso l'apice ma ha ancora la base larga farebbe
# sparire di colpo una base grande (risultato negativo #15, come il #14 per gli scambi).
#
# Le cellule si rinumerano: chi viene dopo la cellula estrusa scala di un posto. Il registro
# dice chi è uscito e chi si è chiuso.

export Estrusione, estrudi!, estrusioni!, rilassa_con_topologia!

"""
    Estrusione(cellula, vicine)

Registro di un'estrusione: è uscita la cellula `cellula` (numerazione di prima) e le `vicine`
(numerazione di dopo) si sono chiuse sul suo posto.
"""
struct Estrusione
    cellula::Int
    vicine::NTuple{3,Int}
end

# Toglie i vertici del piano che nessun poligono usa più e rinumera gli altri.
function _compatta!(t::Tessuto3D, P)
    usati = sort!(unique(j for w in P for (j, _) in w))
    nuovo = Dict(j => k for (k, j) in enumerate(usati))
    idx = reduce(vcat, [[3(2j - 1) - 2, 3(2j - 1) - 1, 3(2j - 1), 3(2j) - 2, 3(2j) - 1, 3(2j)]
                        for j in usati])
    coda = t.periodico ? t.x[end-1:end] : Float64[]
    x = vcat(t.x[idx], coda)
    ancore = vcat(t.ancore[idx], t.periodico ? t.ancore[end-1:end] : Float64[])
    fissi = reduce(vcat, [[t.fissi[2j - 1], t.fissi[2j]] for j in usati])
    Q = [[(nuovo[j], im) for (j, im) in w] for w in P]
    return Q, x, BitVector(fissi), ancore
end

"""
    estrudi!(t, c)

Estrude la cellula `c`, che deve essere un triangolo con i tre vertici interni e non
bloccati. Restituisce il registro, o `nothing` se l'estrusione non è ammessa.
"""
function estrudi!(t::Tessuto3D, c::Int)
    P = poligoni(t)
    w = P[c]
    length(w) == 3 || return nothing
    vicine = Int[]
    for (j, _) in w
        t.fissi[2j - 1] && return nothing
        attorno = cellule_attorno(P, j)
        length(attorno) == 3 || return nothing
        append!(vicine, setdiff(attorno, [c]))
    end
    vicine = unique(vicine)
    length(vicine) == 3 || return nothing
    all(length(P[n]) > 3 for n in vicine) || return nothing      # le vicine restano poligoni
    # il vertice che resta: il primo del triangolo, al baricentro (apicale e basale)
    (j1, i1), (j2, i2), (j3, i3) = w
    Lx, Ly = _periodo(t, t.x)
    for s in (0, 1)                                              # 0 apicale, 1 basale
        v(j) = 2j - 1 + s
        a, b, d = _pos(t, t.x, v(j1), i1), _pos(t, t.x, v(j2), i2), _pos(t, t.x, v(j3), i3)
        g = ((a[1] + b[1] + d[1]) / 3, (a[2] + b[2] + d[2]) / 3, (a[3] + b[3] + d[3]) / 3)
        k = v(j1)
        t.x[3k-2] = g[1] - i1[1] * Lx; t.x[3k-1] = g[2] - i1[2] * Ly; t.x[3k] = g[3]
    end
    # nelle vicine, la coppia di vertici del triangolo diventa il solo j1, con l'immagine giusta
    scarto = Dict(j1 => (0, 0), j2 => _im_meno(i1, i2), j3 => _im_meno(i1, i3))
    for n in vicine
        u = P[n]
        k = findfirst(e -> haskey(scarto, e[1]), u)
        # i due vertici del triangolo in u sono consecutivi (ciclicamente)
        k2 = haskey(scarto, u[mod1(k + 1, length(u))][1]) ? mod1(k + 1, length(u)) : k
        k1 = k2 == k ? mod1(k - 1, length(u)) : k
        jv, iv = u[k1]
        nuovo = (j1, _im_piu(iv, scarto[jv]))
        u[k1] = nuovo
        deleteat!(u, k2)
    end
    deleteat!(P, c)
    intenzioni = deleteat!(copy(t.intenzioni), c)
    Q, x, fissi, ancore = _compatta!(t, P)
    nuovo = tessuto_da_poligoni(intenzioni, Q, x; periodico = t.periodico, fissi, ancore)
    t.x = nuovo.x; t.facce = nuovo.facce; t.cellule = nuovo.cellule
    t.intenzioni = nuovo.intenzioni; t.fissi = nuovo.fissi; t.ancore = nuovo.ancore
    dopo = Tuple(n > c ? n - 1 : n for n in vicine)
    return Estrusione(c, dopo)
end

_area_basale(t::Tessuto3D, c::Int) = _area(_vertici(t, t.x, _faccia_di(t, c, basale)))

"""
    estrusioni!(t; area_soglia)

Estrude ogni cellula triangolare la cui faccia apicale e quella basale sono tutte e due sotto
`area_soglia` (µm²). Restituisce il registro.
"""
function estrusioni!(t::Tessuto3D; area_soglia)
    registro = Estrusione[]
    while true
        c = findfirst(c -> length(t.cellule[c].facce) == 5 &&       # apicale, basale, 3 laterali
                           area_apicale(t, c) < area_soglia && _area_basale(t, c) < area_soglia,
                      eachindex(t.cellule))
        c === nothing && return registro
        r = estrudi!(t, c)
        r === nothing && return registro
        push!(registro, r)
    end
end

"""
    rilassa_con_topologia!(t, p; soglia, lunghezza_nuova, area_soglia, tolleranza,
                           cicli_massimi = 50, iterazioni_per_ciclo = 200)

Alterna Newton, scambi di vicine ed estrusioni finché il tessuto è all'equilibrio e nessuna
operazione topologica è più dovuta. Restituisce `(convergito, cicli, scambi, estrusioni)`.
"""
function rilassa_con_topologia!(t::Tessuto3D, p::ParametriMeccanici; soglia, lunghezza_nuova,
                                area_soglia, tolleranza, cicli_massimi = 50,
                                iterazioni_per_ciclo = 200)
    scambi = ScambioDiVicine[]
    uscite = Estrusione[]
    for ciclo in 1:cicli_massimi
        r = minimizza!(t, p; tolleranza, iterazioni_massime = iterazioni_per_ciclo)
        s = scambi_di_vicine!(t, p; soglia, lunghezza_nuova)
        e = estrusioni!(t; area_soglia)
        append!(scambi, s); append!(uscite, e)
        r.convergito && isempty(s) && isempty(e) &&
            return (convergito = true, cicli = ciclo, scambi, estrusioni = uscite)
    end
    return (convergito = false, cicli = cicli_massimi, scambi, estrusioni = uscite)
end
