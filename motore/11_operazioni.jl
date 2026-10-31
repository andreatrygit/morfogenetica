# Le operazioni topologiche elementari dei modelli 3D a vertici, in forma generale
# (docs/roadmap.md, passo R1; Okuda et al., "reversible network reconnection").
#
# Regioni. Attorno a ogni spigolo ci sono tre regioni e attorno a ogni vertice quattro, dove una
# regione è una cellula oppure il mezzo esterno: una pseudo-cellula per ogni etichetta di faccia
# libera (apicale, basale, laterale). Al bordo di un pezzo di foglio un vertice d'angolo può avere
# solo tre regioni (la cellula, il mezzo apicale e quello laterale).
#
#   spigolo_a_triangolo!   (I→H)  lo spigolo i–j, con attorno A, B, E e agli estremi in più C e
#                                 D, diventa un triangolo fra C e D; i suoi vertici stanno fra
#                                 C, D e le coppie (A,B), (B,E), (E,A)
#   triangolo_a_spigolo!   (H→I)  l'inversa: un triangolo fra C e D diventa uno spigolo con
#                                 attorno le altre tre regioni
#   contrai_spigolo!              uno spigolo del contorno libero di una cellula di bordo (la
#                                 cellula e due mezzi), con almeno un estremo d'angolo, si
#                                 contrae in un punto
#
# Le operazioni non creano mai una faccia fra due mezzi (aprirebbe un buco nel foglio) e, con
# `monostrato = true`, non tolgono a una cellula la sua faccia apicale o basale. Ogni operazione
# è registrata. Per ora solo su tessuti non periodici.

export Operazione, problemi_topologici, spigolo_a_triangolo!, triangolo_a_spigolo!
export contrai_spigolo!, rilassa_con_operazioni!

"""
    Operazione(tipo, vertici, perdono, guadagnano, triangolo)

Registro di un'operazione: `tipo` è :spigolo_a_triangolo, :triangolo_a_spigolo o :contrazione;
`perdono` sono le regioni che non si toccano più, `guadagnano` quelle che si toccano ora;
`triangolo` è l'indice della faccia triangolare creata (0 se nessuna).
"""
struct Operazione
    tipo::Symbol
    vertici::Vector{Int}
    perdono::Vector{Int}
    guadagnano::Vector{Int}
    triangolo::Int
end

# Le pseudo-cellule dei mezzi: numeri negativi.
const _MEZZO = Dict(apicale => -1, basale => -2, laterale => -3)
const _ETICHETTA_MEZZO = Dict(-1 => apicale, -2 => basale, -3 => laterale)
_mezzo(r) = r < 0

# Le due regioni ai lati di una faccia.
_lati(f::Faccia) = length(f.cellule) == 2 ? (f.cellule[1], f.cellule[2]) :
                                            (f.cellule[1], _MEZZO[f.etichetta])

_chiave(a, b) = minmax(a, b)

# Mappe di adiacenza, ricalcolate da capo (i tessuti sono piccoli).
function _adiacenze(t::Tessuto3D)
    facce_di_spigolo = Dict{Tuple{Int,Int},Vector{Int}}()
    regioni_di_vertice = Dict{Int,Set{Int}}()
    facce_di_vertice = Dict{Int,Vector{Int}}()
    for (fi, f) in enumerate(t.facce)
        w = f.vertici
        r1, r2 = _lati(f)
        for k in eachindex(w)
            push!(get!(facce_di_spigolo, _chiave(w[k], w[mod1(k + 1, length(w))]), Int[]), fi)
            union!(get!(regioni_di_vertice, w[k], Set{Int}()), (r1, r2))
            push!(get!(facce_di_vertice, w[k], Int[]), fi)
        end
    end
    return facce_di_spigolo, regioni_di_vertice, facce_di_vertice
end

function _regioni_di_spigolo(t, facce_di_spigolo, a, b)
    R = Set{Int}()
    for fi in get(facce_di_spigolo, _chiave(a, b), Int[])
        union!(R, _lati(t.facce[fi]))
    end
    return R
end

"""
    problemi_topologici(t; monostrato = true)

Elenco (vuoto se va tutto bene) delle violazioni degli invarianti: ogni cellula è un poliedro
chiuso (V − E + F = 2); ogni faccia ha almeno tre vertici distinti; nessuna faccia è doppia;
attorno a ogni spigolo ci sono tre regioni (due dove la superficie libera di una cellula fa
angolo); attorno a ogni vertice quattro, o tre se due sono
mezzi (l'angolo di una cellula di bordo); con `monostrato`, ogni cellula ha una faccia apicale e
una basale. Nei fogli periodici gli spigoli si distinguono anche per l'immagine.
"""
function problemi_topologici(t::Tessuto3D; monostrato = true)
    problemi = String[]
    for c in eachindex(t.cellule)
        χ = caratteristica_di_eulero(t, c)
        χ == 2 || push!(problemi, "cellula $c: V − E + F = $χ")
        if monostrato
            n_ap = count(fi -> t.facce[fi].etichetta == apicale, t.cellule[c].facce)
            n_ba = count(fi -> t.facce[fi].etichetta == basale, t.cellule[c].facce)
            (n_ap == 1 && n_ba == 1) ||
                push!(problemi, "cellula $c: $n_ap facce apicali e $n_ba basali")
        end
    end
    insiemi = Set{Vector{Int}}()
    regioni_spigolo = Dict{Tuple{Int,Int,Int,Int},Set{Int}}()
    regioni_vertice = Dict{Int,Set{Int}}()
    for (fi, f) in enumerate(t.facce)
        w = f.vertici
        length(w) >= 3 || push!(problemi, "faccia $fi con $(length(w)) vertici")
        length(unique(w)) == length(w) || push!(problemi, "faccia $fi con vertici ripetuti")
        s = sort(w)
        s in insiemi && push!(problemi, "faccia $fi doppia")
        push!(insiemi, s)
        r1, r2 = _lati(f)
        for k in eachindex(w)
            k2 = mod1(k + 1, length(w))
            (a, ia), (b, ib) = (w[k], f.immagini[k]), (w[k2], f.immagini[k2])
            d = a < b ? (a, b, ib[1] - ia[1], ib[2] - ia[2]) : (b, a, ia[1] - ib[1], ia[2] - ib[2])
            union!(get!(regioni_spigolo, d, Set{Int}()), (r1, r2))
            union!(get!(regioni_vertice, a, Set{Int}()), (r1, r2))
        end
    end
    for (d, R) in regioni_spigolo
        # tre regioni; due solo dove la superficie libera di una cellula si piega (lo spigolo
        # d'angolo di una cellula di bordo, fra due sue facce verso lo stesso mezzo)
        (length(R) == 3 || (length(R) == 2 && count(_mezzo, R) == 1)) ||
            push!(problemi, "spigolo $(d[1])–$(d[2]) con $(length(R)) regioni attorno")
    end
    for (v, R) in regioni_vertice
        n = length(R)
        mezzi = count(_mezzo, R)
        (n == 4 || (n == 3 && mezzi == 2)) ||
            push!(problemi, "vertice $v con $n regioni attorno ($mezzi mezzi)")
    end
    return problemi
end

# --- aiuti per le modifiche -------------------------------------------------------------------

function _nuova_faccia(f::Faccia, w)
    return Faccia(collect(w), fill((0, 0), length(w)), f.etichetta, copy(f.cellule), copy(f.segni),
                  f.etichetta == apicale ? zeros(Int, length(w)) : Int[])
end

# Sostituisce nel ciclo w il vertice v con la sequenza `nuovi`, scelta in base ai vicini: la
# funzione `per_vicino(u)` dice quale vertice nuovo eredita lo spigolo v–u.
function _sostituisci(w, v, per_vicino)
    i = findfirst(==(v), w)
    n = length(w)
    prec, succ = w[mod1(i - 1, n)], w[mod1(i + 1, n)]
    a, b = per_vicino(prec), per_vicino(succ)
    nuovi = a == b ? [a] : [a, b]
    return vcat(w[1:i-1], nuovi, w[i+1:end])
end

# Sostituisce in w la coppia consecutiva (a, b), in un verso o nell'altro, con il vertice m.
function _fondi_coppia(w, a, b, m)
    n = length(w)
    for i in 1:n
        j = mod1(i + 1, n)
        if (w[i] == a && w[j] == b) || (w[i] == b && w[j] == a)
            u = copy(w)
            u[i] = m
            deleteat!(u, j)
            return u
        end
    end
    throw(ArgumentError("$a e $b non sono consecutivi"))
end

_P(t, k) = _pos(t, t.x, k, (0, 0))

function _nuovo_vertice!(t::Tessuto3D, pos)
    append!(t.x, pos); append!(t.ancore, pos); push!(t.fissi, false)
    return length(t.x) ÷ 3
end

# Toglie il vertice v (che nessuna faccia usa più) e rinumera quelli dopo.
function _togli_vertice!(t::Tessuto3D, v)
    deleteat!(t.x, 3v-2:3v); deleteat!(t.ancore, 3v-2:3v); deleteat!(t.fissi, v)
    for (fi, f) in enumerate(t.facce)
        any(>(v), f.vertici) || continue
        t.facce[fi] = Faccia([k > v ? k - 1 : k for k in f.vertici], f.immagini, f.etichetta,
                             f.cellule, f.segni, f.condivisione)
    end
    return t
end

# Aggiunge una faccia nuova fra le regioni r1 e r2 (non entrambe mezzi), orientata verso
# l'esterno della prima cellula vera.
function _aggiungi_faccia!(t::Tessuto3D, w, r1, r2)
    _mezzo(r1) && ((r1, r2) = (r2, r1))
    etichetta = _mezzo(r2) ? _ETICHETTA_MEZZO[r2] : laterale
    P = [_P(t, k) for k in w]
    n = _croce(_meno(P[2], P[1]), _meno(P[3], P[1]))
    cel = t.cellule[r1]
    cc = _baricentro([_P(t, k) for fi in cel.facce for k in t.facce[fi].vertici])
    w = collect(w)
    _punto(n, _meno(_baricentro(P), cc)) < 0 && reverse!(w)
    cellule = _mezzo(r2) ? [r1] : [r1, r2]
    segni = _mezzo(r2) ? [1] : [1, -1]
    push!(t.facce, Faccia(w, fill((0, 0), length(w)), etichetta, cellule, segni,
                          etichetta == apicale ? zeros(Int, length(w)) : Int[]))
    fi = length(t.facce)
    for (c, s) in zip(cellule, segni)
        cc = t.cellule[c]
        t.cellule[c] = Cellula3D(vcat(cc.facce, fi), vcat(cc.segni, s), vcat(cc.traslazioni, (0, 0)))
    end
    return fi
end

# Toglie la faccia fi (dalle cellule e dall'elenco) e rinumera le facce dopo.
function _togli_faccia!(t::Tessuto3D, fi)
    deleteat!(t.facce, fi)
    for (c, cel) in enumerate(t.cellule)
        tieni = [k for k in eachindex(cel.facce) if cel.facce[k] != fi]
        t.cellule[c] = Cellula3D([g > fi ? g - 1 : g for g in cel.facce[tieni]],
                                 cel.segni[tieni], cel.traslazioni[tieni])
    end
    return t
end

# --- le operazioni ----------------------------------------------------------------------------

"""
    spigolo_a_triangolo!(t, i, j; lunghezza_nuova, monostrato = true)

Lo spigolo i–j diventa un triangolo di lato circa `lunghezza_nuova` (µm), nel piano
perpendicolare allo spigolo e centrato nel suo punto medio. Restituisce il registro, o `nothing`
se l'operazione non è ammessa.
"""
function spigolo_a_triangolo!(t::Tessuto3D, i::Int, j::Int; lunghezza_nuova, monostrato = true)
    t.periodico && return nothing
    (t.fissi[i] || t.fissi[j]) && return nothing
    fs, rv, _ = _adiacenze(t)
    Re = _regioni_di_spigolo(t, fs, i, j)
    length(Re) == 3 || return nothing
    Ri, Rj = rv[i], rv[j]
    (length(Ri) == 4 && length(Rj) == 4) || return nothing
    (issubset(Re, Ri) && issubset(Re, Rj)) || return nothing
    C = only(setdiff(Ri, Re)); D = only(setdiff(Rj, Re))
    C == D && return nothing
    (_mezzo(C) && _mezzo(D)) && return nothing
    # E: un mezzo, se c'è; le tre coppie di (A, B, E) danno i tre vertici nuovi
    tre = sort(collect(Re); by = r -> (_mezzo(r) ? 0 : 1, r))
    E = tre[1]; A, B = tre[2], tre[3]
    coppie = [(A, B), (B, E), (E, A)]

    # posizioni: triangolo nel piano perpendicolare allo spigolo, ogni vertice verso la faccia
    # che lo ha generato
    pi_, pj = _P(t, i), _P(t, j)
    m = (pi_ .+ pj) ./ 2
    e = _meno(pj, pi_); e = e ./ _norma(e)
    r = lunghezza_nuova / sqrt(3)
    direzioni = Dict{Tuple{Int,Int},NTuple{3,Float64}}()
    for fi in fs[_chiave(i, j)]
        x, y = _lati(t.facce[fi])
        g = _meno(_baricentro([_P(t, k) for k in t.facce[fi].vertici]), m)
        g = _meno(g, e .* _punto(g, e))
        direzioni[minmax(x, y)] = g ./ _norma(g)
    end
    for (x, y) in coppie                      # la faccia fra due mezzi non esiste
        if !haskey(direzioni, minmax(x, y))
            altre = [d for (k, d) in direzioni]
            s = reduce((a, b) -> a .+ b, altre)
            direzioni[minmax(x, y)] = .-(s ./ _norma(s))
        end
    end
    pos(x, y) = m .+ direzioni[minmax(x, y)] .* r
    # indici: i e j restano, come vertici (E,A) e (B,E); (A,B) è nuovo
    vAB = _nuovo_vertice!(t, pos(A, B))
    t.x[3i-2:3i] .= pos(E, A)
    t.x[3j-2:3j] .= pos(B, E)
    vert = Dict(minmax(A, B) => vAB, minmax(B, E) => j, minmax(E, A) => i)
    nuovo_di(x, y) = vert[minmax(x, y)]

    # facce che contengono lo spigolo: la coppia (i, j) diventa il vertice della loro coppia
    for fi in fs[_chiave(i, j)]
        x, y = _lati(t.facce[fi])
        t.facce[fi] = _nuova_faccia(t.facce[fi], _fondi_coppia(t.facce[fi].vertici, i, j, nuovo_di(x, y)))
    end
    # facce che contengono uno solo degli estremi: il vertice diventa i due vertici nuovi che
    # ereditano i suoi spigoli, ciascuno dal lato giusto
    _, _, fv = _adiacenze(t)
    for (v, X) in ((i, C), (j, D))
        for fi in copy(get(fv, v, Int[]))
            w = t.facce[fi].vertici
            v in w || continue
            (i in w && j in w) && continue
            per_vicino(u) = begin
                Ru = _regioni_di_spigolo(t, fs, v, u)       # regioni dello spigolo prima
                ab = setdiff(Ru, (X,))
                length(ab) == 2 || return v
                nuovo_di(collect(ab)...)
            end
            t.facce[fi] = _nuova_faccia(t.facce[fi], _sostituisci(w, v, per_vicino))
        end
    end
    # il triangolo nuovo
    tri = _aggiungi_faccia!(t, [vert[minmax(A, B)], vert[minmax(B, E)], vert[minmax(E, A)]], C, D)
    _ricalcola_condivisione!(t)
    if monostrato && !isempty(problemi_topologici(t; monostrato))
        throw(ErrorException("spigolo_a_triangolo!: invarianti violati, " *
                             join(problemi_topologici(t; monostrato), "; ")))
    end
    return Operazione(:spigolo_a_triangolo, [i, j, vAB], [A, B], [C, D], tri)
end

"""
    triangolo_a_spigolo!(t, f; lunghezza_nuova, monostrato = true)

La faccia triangolare `f` fra C e D diventa uno spigolo di lunghezza `lunghezza_nuova` (µm),
perpendicolare al triangolo e centrato nel suo baricentro, con attorno le altre tre regioni.
Restituisce il registro, o `nothing` se l'operazione non è ammessa.
"""
function triangolo_a_spigolo!(t::Tessuto3D, f::Int; lunghezza_nuova, monostrato = true)
    t.periodico && return nothing
    T = t.facce[f]
    length(T.vertici) == 3 || return nothing
    any(k -> t.fissi[k], T.vertici) && return nothing
    monostrato && T.etichetta != laterale && return nothing   # la faccia apicale o basale resta
    C, D = _lati(T)
    fs, rv, fv = _adiacenze(t)
    v1, v2, v3 = T.vertici
    all(length(rv[v]) == 4 for v in (v1, v2, v3)) || return nothing
    altre = setdiff(union(rv[v1], rv[v2], rv[v3]), (C, D))
    length(altre) == 3 || return nothing
    count(_mezzo, altre) >= 2 && return nothing          # lo spigolo starebbe fra due mezzi

    # posizioni: spigolo lungo la normale del triangolo, i dalla parte di C
    P = [_P(t, k) for k in T.vertici]
    g = _baricentro(P)
    n = _croce(_meno(P[2], P[1]), _meno(P[3], P[1])); n = n ./ _norma(n)
    cel = t.cellule[_mezzo(C) ? D : C]
    cc = _baricentro([_P(t, k) for fi in cel.facce for k in t.facce[fi].vertici])
    verso_C = _punto(n, _meno(cc, g)) > 0
    lato_c = _mezzo(C) ? !verso_C : verso_C           # n punta verso C?
    pc = g .+ n .* (lato_c ? 1 : -1) .* (lunghezza_nuova / 2)
    pd = g .- n .* (lato_c ? 1 : -1) .* (lunghezza_nuova / 2)
    i, j, morto = v1, v2, v3
    t.x[3i-2:3i] .= pc
    t.x[3j-2:3j] .= pd

    # le facce che contengono due vertici del triangolo (fuori dal triangolo) stanno fra C (o D)
    # e una delle altre regioni: la coppia diventa i (o j)
    for (fi, F) in enumerate(t.facce)
        fi == f && continue
        dentro = [k for k in F.vertici if k in (v1, v2, v3)]
        length(dentro) == 2 || continue
        x, y = _lati(F)
        nuovo = (x == C || y == C) ? i : j
        t.facce[fi] = _nuova_faccia(F, _fondi_coppia(F.vertici, dentro..., nuovo))
    end
    # le facce che contengono un solo vertice del triangolo stanno fra due delle altre regioni:
    # il vertice diventa i e j, ciascuno dal lato giusto
    for (fi, F) in enumerate(t.facce)
        fi == f && continue
        dentro = [k for k in F.vertici if k in (v1, v2, v3)]
        length(dentro) == 1 || continue
        v = only(dentro)
        per_vicino(u) = begin
            Ru = _regioni_di_spigolo(t, fs, v, u)
            C in Ru ? i : D in Ru ? j : v
        end
        nuovo_w = _sostituisci(F.vertici, v, per_vicino)
        t.facce[fi] = _nuova_faccia(F, nuovo_w)
    end
    _togli_faccia!(t, f)
    any(k -> k == morto, (k for F in t.facce for k in F.vertici)) &&
        throw(ErrorException("triangolo_a_spigolo!: il terzo vertice è ancora usato"))
    _togli_vertice!(t, morto)
    _ricalcola_condivisione!(t)
    if monostrato && !isempty(problemi_topologici(t; monostrato))
        throw(ErrorException("triangolo_a_spigolo!: invarianti violati, " *
                             join(problemi_topologici(t; monostrato), "; ")))
    end
    i2, j2 = i > morto ? i - 1 : i, j > morto ? j - 1 : j
    return Operazione(:triangolo_a_spigolo, [i2, j2], [C, D], collect(altre), 0)
end

"""
    contrai_spigolo!(t, i, j; monostrato = true)

Contrae in un punto lo spigolo i–j del contorno libero di una cellula di bordo: le tre regioni
dello spigolo sono la cellula e due mezzi, e almeno un estremo (l'angolo) non ha altre regioni.
L'angolo si fonde nell'altro estremo, che non si muove; se sono angoli tutti e due, il punto
che resta è quello medio. Le facce che contengono lo spigolo perdono un vertice. Restituisce il
registro, o `nothing`.
"""
function contrai_spigolo!(t::Tessuto3D, i::Int, j::Int; monostrato = true)
    t.periodico && return nothing
    (t.fissi[i] || t.fissi[j]) && return nothing
    fs, rv, _ = _adiacenze(t)
    Re = _regioni_di_spigolo(t, fs, i, j)
    (length(Re) == 3 && count(_mezzo, Re) == 2) || return nothing
    rv[j] == Re || ((i, j) = (j, i))                 # j è un angolo e sparisce
    rv[j] == Re || return nothing
    all(length(t.facce[fi].vertici) > 3 for fi in fs[_chiave(i, j)]) || return nothing
    rv[i] == Re && (t.x[3i-2:3i] .= (_P(t, i) .+ _P(t, j)) ./ 2)
    for fi in fs[_chiave(i, j)]
        t.facce[fi] = _nuova_faccia(t.facce[fi], _fondi_coppia(t.facce[fi].vertici, i, j, i))
    end
    for (fi, F) in enumerate(t.facce)
        j in F.vertici && (t.facce[fi] = _nuova_faccia(F, [k == j ? i : k for k in F.vertici]))
    end
    _togli_vertice!(t, j)
    _ricalcola_condivisione!(t)
    if monostrato && !isempty(problemi_topologici(t; monostrato))
        throw(ErrorException("contrai_spigolo!: invarianti violati, " *
                             join(problemi_topologici(t; monostrato), "; ")))
    end
    return Operazione(:contrazione, [i > j ? i - 1 : i], Int[], Int[], 0)
end

"""
    rilassa_con_operazioni!(t, p; soglia, area_soglia, lunghezza_nuova, tolleranza,
                            cicli_massimi = 30, iterazioni_per_ciclo = 200, monostrato = true)

Alterna Newton e le operazioni elementari finché il tessuto è all'equilibrio e nessuna
operazione è più dovuta: i triangoli sotto `area_soglia` (µm²) diventano spigoli, gli spigoli
sotto `soglia` (µm) diventano triangoli o si contraggono. Gli elementi appena creati sono più
grandi delle soglie (`lunghezza_nuova` deve superare `soglia` e dare un triangolo d'area maggiore
di `area_soglia`), così non si annullano subito. Restituisce `(convergito, cicli, registro)`.
"""
function rilassa_con_operazioni!(t::Tessuto3D, p::ParametriMeccanici; soglia, area_soglia,
                                 lunghezza_nuova, tolleranza, cicli_massimi = 30,
                                 iterazioni_per_ciclo = 200, monostrato = true)
    lunghezza_nuova > soglia || throw(ArgumentError("lunghezza_nuova deve superare soglia"))
    sqrt(3) / 4 * lunghezza_nuova^2 > area_soglia ||
        throw(ArgumentError("il triangolo nuovo deve superare area_soglia"))
    registro = Operazione[]
    for ciclo in 1:cicli_massimi
        r = minimizza!(t, p; tolleranza, iterazioni_massime = iterazioni_per_ciclo)
        nuovi = _operazioni_dovute!(t; soglia, area_soglia, lunghezza_nuova, monostrato)
        append!(registro, nuovi)
        r.convergito && isempty(nuovi) && return (convergito = true, cicli = ciclo, registro)
    end
    return (convergito = false, cicli = cicli_massimi, registro)
end

function _operazioni_dovute!(t; soglia, area_soglia, lunghezza_nuova, monostrato)
    fatti = Operazione[]
    tentati = Set{Any}()
    while true
        fatto = false
        # triangoli piccoli → spigoli
        for (fi, F) in enumerate(t.facce)
            length(F.vertici) == 3 || continue
            chiave = (:t, sort(F.vertici))
            chiave in tentati && continue
            area_faccia(t, F) < area_soglia || continue
            push!(tentati, chiave)
            r = triangolo_a_spigolo!(t, fi; lunghezza_nuova, monostrato)
            r === nothing && continue
            push!(fatti, r); fatto = true
            break
        end
        fatto && continue
        # spigoli corti → triangoli, o contrazioni al bordo
        corto = nothing; Lmin = soglia
        for F in t.facce, k in eachindex(F.vertici)
            a, b = F.vertici[k], F.vertici[mod1(k + 1, length(F.vertici))]
            chiave = (:s, minmax(a, b))
            chiave in tentati && continue
            L = _norma(_meno(_P(t, b), _P(t, a)))
            L < Lmin && (corto = (a, b); Lmin = L)
        end
        corto === nothing && return fatti
        push!(tentati, (:s, minmax(corto...)))
        r = spigolo_a_triangolo!(t, corto...; lunghezza_nuova, monostrato)
        r === nothing && (r = contrai_spigolo!(t, corto...; monostrato))
        r === nothing || push!(fatti, r)
    end
end
