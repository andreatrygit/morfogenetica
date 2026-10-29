# Scambi di vicine (T1) nel guscio 3D: specifica, decisione 4 — cambio topologico controllato.
#
# Quando lo spigolo apicale fra due cellule A e B si accorcia sotto una soglia, lo si sostituisce
# con uno spigolo nuovo, perpendicolare, fra le altre due cellule che lo toccavano, C e D:
#
#         C                    C
#       \   /                  |
#     A  p─q  B     →     A    p    B       (e la stessa cosa sulla superficie basale)
#       /   \                  q
#         D                    |
#                              D
#
# A e B perdono un lato, C e D ne guadagnano uno; la faccia laterale fra A e B diventa una faccia
# fra C e D. Si fa insieme sulla superficie apicale e su quella basale, così ogni cellula resta un
# prisma. È un'operazione discreta, esplicita, con le sue condizioni, e ogni scambio viene
# registrato: niente cambi di topologia impliciti.
#
# Il tessuto scorre invece di strapparsi: è la lacerazione del solco del motore vecchio
# (specifica, sezione 2), che va evitata per costruzione.

export ScambioDiVicine, scambia!, scambi_di_vicine!, rilassa_con_scambi!, cellule_attorno

"""
    ScambioDiVicine(p, q, perdono, guadagnano)

Registro di uno scambio: lo spigolo fra i vertici del piano `p` e `q` è passato dalle cellule
`perdono` = (A, B) alle cellule `guadagnano` = (C, D).
"""
struct ScambioDiVicine
    p::Int
    q::Int
    perdono::NTuple{2,Int}
    guadagnano::NTuple{2,Int}
end

_im_meno(a, b) = (a[1] - b[1], a[2] - b[2])
_im_piu(a, b) = (a[1] + b[1], a[2] + b[2])

# Dove sta il vertice del piano j nel poligono w (indice), o nothing.
_dove(w, j) = findfirst(e -> e[1] == j, w)

"""
    cellule_attorno(P, j)

Le cellule il cui poligono contiene il vertice del piano j.
"""
cellule_attorno(P, j) = [c for c in eachindex(P) if _dove(P[c], j) !== nothing]

# Le quattro cellule di uno spigolo interno p–q, o nothing se lo scambio non è ammesso.
# A è la cellula in cui q segue p (senso antiorario), B quella in cui p segue q.
function _quartetto(t::Tessuto3D, P, p, q)
    (t.fissi[2p-1] || t.fissi[2q-1]) && return nothing
    cp, cq = cellule_attorno(P, p), cellule_attorno(P, q)
    (length(cp) == 3 && length(cq) == 3) || return nothing   # solo spigoli interni
    AB = intersect(cp, cq)
    length(AB) == 2 || return nothing
    A = B = 0
    for c in AB
        w = P[c]; k = _dove(w, p)
        w[mod1(k + 1, length(w))][1] == q ? (A = c) : (B = c)
    end
    (A == 0 || B == 0) && return nothing
    C = only(setdiff(cp, AB))
    D = only(setdiff(cq, AB))
    C == D && return nothing
    (length(P[A]) > 3 && length(P[B]) > 3) || return nothing  # restano almeno triangoli
    return A, B, C, D
end

"""
    scambia!(t, p, q; lunghezza_nuova)

Esegue lo scambio di vicine sullo spigolo fra i vertici del piano `p` e `q` (apicali 2p-1, 2q-1;
basali 2p, 2q). I due vertici vengono rimessi a distanza `lunghezza_nuova` (µm), attorno al
punto medio, nella direzione perpendicolare al vecchio spigolo. Restituisce il registro, o
`nothing` se lo spigolo non ammette lo scambio (bordo, vertici bloccati, cellule che
diventerebbero meno che triangoli).
"""
function scambia!(t::Tessuto3D, p::Int, q::Int; lunghezza_nuova)
    P = poligoni(t)
    quattro = _quartetto(t, P, p, q)
    quattro === nothing && return nothing
    A, B, C, D = quattro
    # immagini nel riferimento di A e scarto fra q e p
    ip = P[A][_dove(P[A], p)][2]
    iq = P[A][_dove(P[A], q)][2]
    δ = _im_meno(iq, ip)

    # nuove posizioni, superficie per superficie (apicale e basale)
    centroA = _baricentro(_vertici(t, t.x, _faccia_di(t, A, apicale)))
    for (vp, vq) in ((2p - 1, 2q - 1), (2p, 2q))
        a = _pos(t, t.x, vp, ip)
        b = _pos(t, t.x, vq, iq)
        m = ((a[1] + b[1]) / 2, (a[2] + b[2]) / 2, (a[3] + b[3]) / 2)
        e = _meno(b, a)
        # normale locale al foglio: dal basale all'apicale, media dei due vertici
        ax = _meno(_pos(t, t.x, 2p - 1, ip), _pos(t, t.x, 2p, ip))
        bx = _meno(_pos(t, t.x, 2q - 1, iq), _pos(t, t.x, 2q, iq))
        nrm = (ax[1] + bx[1], ax[2] + bx[2], ax[3] + bx[3])
        d = _croce(nrm, e)
        L = _norma(d)
        d = (d[1] / L * lunghezza_nuova / 2, d[2] / L * lunghezza_nuova / 2, d[3] / L * lunghezza_nuova / 2)
        # p va dalla parte di A (resta in A), q dalla parte di B
        verso = _punto(d, _meno(centroA, m)) >= 0 ? 1.0 : -1.0
        np = (m[1] + verso * d[1], m[2] + verso * d[2], m[3] + verso * d[3])
        nq = (m[1] - verso * d[1], m[2] - verso * d[2], m[3] - verso * d[3])
        Lx, Ly = _periodo(t, t.x)
        t.x[3vp-2] = np[1] - ip[1] * Lx; t.x[3vp-1] = np[2] - ip[2] * Ly; t.x[3vp] = np[3]
        t.x[3vq-2] = nq[1] - iq[1] * Lx; t.x[3vq-1] = nq[2] - iq[2] * Ly; t.x[3vq] = nq[3]
    end

    # nuova topologia
    wA = P[A]; deleteat!(wA, _dove(wA, q))                       # A: …, a_prec, p, a_succ, …
    wB = P[B]; deleteat!(wB, _dove(wB, p))                       # B: …, b_prec, q, b_succ, …
    wC = P[C]; kC = _dove(wC, p)                                 # C: …, b_succ, p, a_prec, …
    insert!(wC, kC, (q, _im_piu(wC[kC][2], δ)))                  #  → …, b_succ, q, p, a_prec, …
    wD = P[D]; kD = _dove(wD, q)                                 # D: …, a_succ, q, b_prec, …
    insert!(wD, kD, (p, _im_meno(wD[kD][2], δ)))                 #  → …, a_succ, p, q, b_prec, …

    nuovo = tessuto_da_poligoni(t.intenzioni, P, t.x; periodico = t.periodico,
                                fissi = t.fissi, ancore = t.ancore)
    t.facce = nuovo.facce
    t.cellule = nuovo.cellule
    return ScambioDiVicine(p, q, (A, B), (C, D))
end

"""
    scambi_di_vicine!(t, p; soglia, lunghezza_nuova)

Esegue uno scambio su ogni parete laterale interna i cui due spigoli, apicale **e** basale, sono
più corti di `soglia` (µm), dalla più corta, e restituisce il registro degli scambi fatti.
Servono tutti e due: nel guscio a prismi una cellula ha gli stessi vicini sopra e sotto, quindi
lo scambio cambia insieme le due superfici, e farlo quando solo lo spigolo apicale è corto
schiaccerebbe uno spigolo basale lungo (risultato negativo #14). `lunghezza_nuova` deve superare `soglia`, così uno
spigolo appena creato non viene scambiato di nuovo nella stessa chiamata; `soglia` deve
superare la lunghezza minima ammessa, così lo scambio arriva prima del rifiuto.
"""
function scambi_di_vicine!(t::Tessuto3D, p::ParametriMeccanici; soglia, lunghezza_nuova)
    lunghezza_nuova > soglia || throw(ArgumentError("lunghezza_nuova deve superare soglia"))
    soglia > p.lunghezza_minima || throw(ArgumentError("soglia deve superare lunghezza_minima"))
    registro = ScambioDiVicine[]
    fatti = Set{Tuple{Int,Int}}()
    while true
        corto = nothing
        Lmin = soglia
        for w in poligoni(t), k in eachindex(w)
            (j1, i1), (j2, i2) = w[k], w[mod1(k + 1, length(w))]
            minmax(j1, j2) in fatti && continue
            La = _norma(_meno(_pos(t, t.x, 2j2 - 1, i2), _pos(t, t.x, 2j1 - 1, i1)))
            Lb = _norma(_meno(_pos(t, t.x, 2j2, i2), _pos(t, t.x, 2j1, i1)))
            L = max(La, Lb)
            if L < Lmin && _quartetto(t, poligoni(t), j1, j2) !== nothing
                corto, Lmin = (j1, j2), L
            end
        end
        corto === nothing && return registro
        r = scambia!(t, corto...; lunghezza_nuova)
        push!(fatti, minmax(corto...))
        r === nothing || push!(registro, r)
    end
end

"""
    rilassa_con_scambi!(t, p; soglia, lunghezza_nuova, tolleranza, cicli_massimi = 50,
                        iterazioni_per_ciclo = 200)

Alterna la minimizzazione dell'energia (Newton, `minimizza!`) e gli scambi di vicine finché il
tessuto è all'equilibrio e nessuno spigolo è sotto soglia. Restituisce
`(convergito, cicli, registro)`.
"""
function rilassa_con_scambi!(t::Tessuto3D, p::ParametriMeccanici; soglia, lunghezza_nuova,
                             tolleranza, cicli_massimi = 50, iterazioni_per_ciclo = 200)
    registro = ScambioDiVicine[]
    for ciclo in 1:cicli_massimi
        r = minimizza!(t, p; tolleranza, iterazioni_massime = iterazioni_per_ciclo)
        nuovi = scambi_di_vicine!(t, p; soglia, lunghezza_nuova)
        append!(registro, nuovi)
        r.convergito && isempty(nuovi) && return (convergito = true, cicli = ciclo, registro)
    end
    return (convergito = false, cicli = cicli_massimi, registro)
end
