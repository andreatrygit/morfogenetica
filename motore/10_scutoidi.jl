# Scutoidi: cellule con vicini diversi all'apice e alla base (docs/piano.md, sezione 12).
#
# Nel guscio una cellula resta un monostrato (va dalla superficie apicale a quella basale), ma
# non deve più essere un prisma: è un poliedro qualunque. Energia, aree, volumi e forze di
# 07_guscio.jl lavorano già su facce poligonali qualsiasi; qui ci sono le operazioni
# topologiche che producono poliedri non prismatici.
#
# Scambio apicale con scutoide. Quando si chiude solo la giunzione apicale fra le cellule A e B
# (la basale resta lunga):
#
#   - in alto lo spigolo passa a C e D, come in uno scambio di vicine;
#   - la parete fra A e B resta, ma in alto finisce in un vertice nuovo v, interno al foglio,
#     dove si incontrano A, B, C e D;
#   - sopra v nasce una parete triangolare fra C e D: lo spigolo apicale nuovo e v.
#
# È il cambio di vicini a metà altezza che fa uno scutoide. Le giunzioni che scendevano dagli
# estremi dello spigolo ora partono da v. Ogni scambio è registrato. Per ora solo su tessuti non
# periodici (le immagini dei vertici nuovi sono da aggiungere).

export ScambioApicale, scambio_apicale!, scutoidi_apicali!, rilassa_con_scutoidi!
export caratteristica_di_eulero

"""
    ScambioApicale(a1, a2, v, perdono, guadagnano)

Registro di uno scambio apicale con scutoide: lo spigolo apicale fra i vertici `a1` e `a2` è
passato dalle cellule `perdono` = (A, B) a `guadagnano` = (C, D), e la parete fra A e B finisce
ora nel vertice interno `v`.
"""
struct ScambioApicale
    a1::Int
    a2::Int
    v::Int
    perdono::NTuple{2,Int}
    guadagnano::NTuple{2,Int}
end

_faccia_apicale_indice(t::Tessuto3D, c) =
    only(fi for fi in t.cellule[c].facce if t.facce[fi].etichetta == apicale)

# Posizione di k in un ciclo e i suoi due vicini nel ciclo.
function _vicini_nel_ciclo(w, k)
    i = findfirst(==(k), w)
    i === nothing && return nothing
    return w[mod1(i - 1, length(w))], w[mod1(i + 1, length(w))]
end

# Inserisce `nuovo` fra a e b, che devono essere consecutivi nel ciclo (in un verso o nell'altro).
function _inserisci_fra(w, a, b, nuovo)
    n = length(w)
    for i in 1:n
        j = mod1(i + 1, n)
        if (w[i] == a && w[j] == b) || (w[i] == b && w[j] == a)
            return insert!(copy(w), i + 1, nuovo)
        end
    end
    throw(ArgumentError("$a e $b non sono consecutivi"))
end

_rinomina(w, da, a) = [k == da ? a : k for k in w]

function _con_vertici(f::Faccia, w)
    return Faccia(w, fill((0, 0), length(w)), f.etichetta, copy(f.cellule), copy(f.segni),
                  f.etichetta == apicale ? zeros(Int, length(w)) : Int[])
end

# Conta, per ogni spigolo apicale, quante facce apicali lo hanno (la tensione di linea va
# divisa fra loro) e lo scrive nelle facce apicali.
function _ricalcola_condivisione!(t::Tessuto3D)
    conta = Dict{Tuple{Int,Int},Int}()
    for f in t.facce
        f.etichetta == apicale || continue
        n = length(f.vertici)
        for k in 1:n
            e = minmax(f.vertici[k], f.vertici[mod1(k + 1, n)])
            conta[e] = get(conta, e, 0) + 1
        end
    end
    for (fi, f) in enumerate(t.facce)
        f.etichetta == apicale || continue
        n = length(f.vertici)
        f.condivisione .= [conta[minmax(f.vertici[k], f.vertici[mod1(k + 1, n)])] for k in 1:n]
    end
    return t
end

"""
    caratteristica_di_eulero(t, c)

V − E + F della superficie della cellula `c`: vale 2 per ogni poliedro chiuso senza buchi.
"""
function caratteristica_di_eulero(t::Tessuto3D, c::Int)
    V = Set{Int}(); E = Set{Tuple{Int,Int}}()
    F = length(t.cellule[c].facce)
    for fi in t.cellule[c].facce
        w = t.facce[fi].vertici
        union!(V, w)
        for k in eachindex(w)
            push!(E, minmax(w[k], w[mod1(k + 1, length(w))]))
        end
    end
    return length(V) - length(E) + F
end

"""
    scambio_apicale!(t, a1, a2; lunghezza_nuova, profondita = lunghezza_nuova)

Scambio apicale con scutoide sullo spigolo apicale fra i vertici `a1` e `a2`. I due vertici
vanno a distanza `lunghezza_nuova` (µm) nella direzione perpendicolare al vecchio spigolo; il
vertice interno nuovo sta `profondita` (µm) sotto il punto medio, verso la parete fra A e B.
Restituisce il registro, o `nothing` se lo scambio non è ammesso.
"""
function scambio_apicale!(t::Tessuto3D, a1::Int, a2::Int; lunghezza_nuova,
                          profondita = lunghezza_nuova)
    t.periodico && return nothing
    (t.fissi[a1] || t.fissi[a2]) && return nothing
    ap = [_faccia_apicale_indice(t, c) for c in eachindex(t.cellule)]
    W = [t.facce[fi].vertici for fi in ap]
    c1 = [c for c in eachindex(W) if a1 in W[c]]
    c2 = [c for c in eachindex(W) if a2 in W[c]]
    (length(c1) == 3 && length(c2) == 3) || return nothing
    AB = intersect(c1, c2)
    length(AB) == 2 || return nothing
    A = B = 0
    for c in AB
        _, dopo = _vicini_nel_ciclo(W[c], a1)
        dopo == a2 ? (A = c) : (B = c)
    end
    (A == 0 || B == 0) && return nothing
    C = only(setdiff(c1, AB)); D = only(setdiff(c2, AB))
    C == D && return nothing
    (length(W[A]) > 3 && length(W[B]) > 3) || return nothing

    # facce laterali coinvolte
    lat(fi) = t.facce[fi].etichetta == laterale
    fAB = findfirst(fi -> lat(fi) && A in t.facce[fi].cellule && B in t.facce[fi].cellule &&
                          a1 in t.facce[fi].vertici && a2 in t.facce[fi].vertici, eachindex(t.facce))
    fAB === nothing && return nothing
    trova(x, y, a) = findfirst(fi -> lat(fi) && fi != fAB && x in t.facce[fi].cellule &&
                                     y in t.facce[fi].cellule && a in t.facce[fi].vertici,
                               eachindex(t.facce))
    fAC, fBC, fAD, fBD = trova(A, C, a1), trova(B, C, a1), trova(A, D, a2), trova(B, D, a2)
    any(isnothing, (fAC, fBC, fAD, fBD)) && return nothing
    wAB = t.facce[fAB].vertici
    d1 = only(setdiff(collect(_vicini_nel_ciclo(wAB, a1)), [a2]))
    d2 = only(setdiff(collect(_vicini_nel_ciclo(wAB, a2)), [a1]))

    # posizioni
    P(k) = _pos(t, t.x, k, (0, 0))
    pa1, pa2 = P(a1), P(a2)
    m = ((pa1[1] + pa2[1]) / 2, (pa1[2] + pa2[2]) / 2, (pa1[3] + pa2[3]) / 2)
    e = _meno(pa2, pa1)
    giu = (P(d1) .+ P(d2)) ./ 2
    nrm = _meno(m, giu)                       # dal basso verso l'apice
    d = _croce(nrm, e)
    d = d ./ _norma(d) .* (lunghezza_nuova / 2)
    cA = _baricentro(_vertici(t, t.x, t.facce[ap[A]]))
    verso = _punto(d, _meno(cA, m)) >= 0 ? 1.0 : -1.0
    np, nq = m .+ verso .* d, m .- verso .* d
    dir = _meno(giu, m)
    nv = m .+ dir ./ _norma(dir) .* profondita
    t.x[3a1-2:3a1] .= np
    t.x[3a2-2:3a2] .= nq
    v = length(t.x) ÷ 3 + 1
    append!(t.x, nv); append!(t.ancore, nv); push!(t.fissi, false)

    # topologia: facce apicali (prima C e D, che usano i poligoni vecchi di A e B)
    xB = only(k for k in _vicini_nel_ciclo(W[C], a1) if k in W[B])
    yA = only(k for k in _vicini_nel_ciclo(W[D], a2) if k in W[A])
    F = t.facce
    F[ap[C]] = _con_vertici(F[ap[C]], _inserisci_fra(W[C], xB, a1, a2))
    F[ap[D]] = _con_vertici(F[ap[D]], _inserisci_fra(W[D], yA, a2, a1))
    F[ap[A]] = _con_vertici(F[ap[A]], filter(!=(a2), W[A]))
    F[ap[B]] = _con_vertici(F[ap[B]], filter(!=(a1), W[B]))
    # facce laterali
    i = findfirst(==(a1), wAB)
    nuovaAB = copy(wAB); nuovaAB[i] = v
    deleteat!(nuovaAB, findfirst(==(a2), nuovaAB))
    F[fAB] = _con_vertici(F[fAB], nuovaAB)
    F[fAC] = _con_vertici(F[fAC], _inserisci_fra(F[fAC].vertici, a1, d1, v))
    F[fBC] = _con_vertici(F[fBC], _inserisci_fra(_rinomina(F[fBC].vertici, a1, a2), a2, d1, v))
    F[fAD] = _con_vertici(F[fAD], _inserisci_fra(_rinomina(F[fAD].vertici, a2, a1), a1, d2, v))
    F[fBD] = _con_vertici(F[fBD], _inserisci_fra(F[fBD].vertici, a2, d2, v))
    # parete nuova fra C e D, orientata verso l'esterno di C
    tri = [a1, a2, v]
    Pt = [P(k) for k in tri]
    nt = _croce(_meno(Pt[2], Pt[1]), _meno(Pt[3], Pt[1]))
    cC = _baricentro([P(k) for fi in t.cellule[C].facce for k in F[fi].vertici])
    _punto(nt, _meno(_baricentro(Pt), cC)) < 0 && reverse!(tri)
    push!(F, Faccia(tri, fill((0, 0), 3), laterale, [C, D], [1, -1], Int[]))
    fCD = length(F)
    for (c, s) in ((C, 1), (D, -1))
        cel = t.cellule[c]
        t.cellule[c] = Cellula3D(vcat(cel.facce, fCD), vcat(cel.segni, s), vcat(cel.traslazioni, (0, 0)))
    end
    _ricalcola_condivisione!(t)
    return ScambioApicale(a1, a2, v, (A, B), (C, D))
end

"""
    scutoidi_apicali!(t, p; soglia, lunghezza_nuova)

Fa uno scambio apicale con scutoide su ogni spigolo apicale interno più corto di `soglia`
(µm), dal più corto. Restituisce il registro.
"""
function scutoidi_apicali!(t::Tessuto3D, p::ParametriMeccanici; soglia, lunghezza_nuova)
    lunghezza_nuova > soglia || throw(ArgumentError("lunghezza_nuova deve superare soglia"))
    soglia > p.lunghezza_minima || throw(ArgumentError("soglia deve superare lunghezza_minima"))
    registro = ScambioApicale[]
    fatti = Set{Tuple{Int,Int}}()
    while true
        corto = nothing; Lmin = soglia
        for f in t.facce
            f.etichetta == apicale || continue
            w = f.vertici
            for k in eachindex(w)
                e = minmax(w[k], w[mod1(k + 1, length(w))])
                e in fatti && continue
                L = _norma(_meno(_pos(t, t.x, e[2], (0, 0)), _pos(t, t.x, e[1], (0, 0))))
                L < Lmin && (corto = e; Lmin = L)
            end
        end
        corto === nothing && return registro
        push!(fatti, corto)
        r = scambio_apicale!(t, corto...; lunghezza_nuova)
        r === nothing || push!(registro, r)
    end
end

"""
    rilassa_con_scutoidi!(t, p; soglia, lunghezza_nuova, tolleranza, cicli_massimi = 30,
                          iterazioni_per_ciclo = 200)

Alterna Newton e scambi apicali con scutoide finché il tessuto è all'equilibrio e nessuno
spigolo apicale è sotto soglia. Restituisce `(convergito, cicli, registro)`.
"""
function rilassa_con_scutoidi!(t::Tessuto3D, p::ParametriMeccanici; soglia, lunghezza_nuova,
                               tolleranza, cicli_massimi = 30, iterazioni_per_ciclo = 200)
    registro = ScambioApicale[]
    for ciclo in 1:cicli_massimi
        r = minimizza!(t, p; tolleranza, iterazioni_massime = iterazioni_per_ciclo)
        nuovi = scutoidi_apicali!(t, p; soglia, lunghezza_nuova)
        append!(registro, nuovi)
        r.convergito && isempty(nuovi) && return (convergito = true, cicli = ciclo, registro)
    end
    return (convergito = false, cicli = cicli_massimi, registro)
end
