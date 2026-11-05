# La matrice delle derivate seconde esatta (docs/roadmap.md, R6): per il solutore sparso
# (15_sparso.jl) al posto della matrice per differenze del gradiente.
#
# L'energia del guscio (07_guscio.jl) è una somma di termini locali:
#
#   faccia per faccia     σ_f · A_f (+ Σ Λ/condivisione · lunghezza dello spigolo, se apicale)
#   cellula per cellula   K/2 · (V_c − V0_c)²,   V_c = Σ_facce ± v_f
#
# e ogni termine dipende solo dai vertici della sua faccia o della sua cellula (più le due voci
# della scatola, nei fogli periodici). Il gradiente di ogni termine di faccia (di σA e del
# contributo di volume v_f) è scritto in forma chiusa qui sotto, come in `gradiente!` di
# 07_guscio.jl ma sulle sole posizioni della faccia; la sua matrice jacobiana, cioè la matrice
# delle derivate seconde, si calcola per differenziazione automatica in avanti (ForwardDiff),
# esatta a precisione di macchina. Le derivate seconde scritte a mano sarebbero un campo minato:
# il centro di ogni faccia è il baricentro del contorno, pesato con le lunghezze dei lati.
# Differenziare due volte l'energia (ForwardDiff.hessian) costerebbe dieci volte tanto.
#
# Per la cellula: K·(∇V ∇Vᵀ + (V − V0)·Σ ±∇²v_f).
#
# Le derivate sono prese rispetto alle posizioni dei vertici della faccia, P = q + immagine·L,
# e portate sulle incognite con la regola della catena (lineare). I vertici bloccati sono
# costanti: le loro righe e colonne restano vuote, come nella matrice per differenze. I conti di
# faccia sono su tuple di lunghezza fissa (niente memoria temporanea); ogni cellula raccoglie in
# una matrice piena il suo volume e le facce di cui è la prima cellula, e le cellule si calcolano
# in parallelo, a blocchi raccolti in ordine fisso: la matrice è la stessa qualunque sia il
# numero di fili. Alla fine si simmetrizza.

using ForwardDiff

_più(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
_per(s, a) = (s * a[1], s * a[2], s * a[3])

# Il gradiente rispetto al centro c, `gc`, ripartito sul vertice k (come `_ripartisci_centro!`):
# Σ sui due spigoli di k di [u ((m − c)·gc) + (l/2) gc] / Ltot.
function _centro_su(P::NTuple{n}, c, gc, Ltot, k) where {n}
    Ltot == 0 && return _per(1 / n, gc)
    tot = _per(zero(Ltot), gc)
    for altro in (P[mod1(k - 1, n)], P[mod1(k + 1, n)])
        d = _meno(P[k], altro)
        l = _norma(d)
        l == 0 && continue
        m = _per(1 / 2, _più(P[k], altro))
        tot = _più(tot, _più(_per(_punto(_meno(m, c), gc) / l, d), _per(l / 2, gc)))
    end
    return _per(1 / Ltot, tot)
end

_lunghezza_contorno(P::NTuple{n}) where {n} = sum(_norma(_meno(P[mod1(k + 1, n)], P[k])) for k in 1:n)

_appiattisci(G::NTuple{n}) where {n} = ntuple(i -> G[(i - 1) ÷ 3 + 1][(i - 1) % 3 + 1], Val(3n))

# Il gradiente di σ·A (+ l'anello apicale, con `Λc[k]` = Λ/condivisione dello spigolo k, oppure
# `nothing`) rispetto alle posizioni P della faccia, appiattito.
function _gradiente_area(P::NTuple{n,NTuple{3,T}}, σ, Λc) where {n,T}
    c = _baricentro(P)
    z3 = (zero(T), zero(T), zero(T))
    # triangolo (a, b, c): ∂A/∂a = ½ (b−c)×n̂, ∂A/∂b = ½ n̂×(a−c)
    tri = ntuple(Val(n)) do k
        a, b = P[k], P[mod1(k + 1, n)]
        m = _croce(_meno(a, c), _meno(b, c))
        L = _norma(m)
        L > 0 || return (z3, z3)
        nh = _per(1 / L, m)
        return (_per(σ / 2, _croce(_meno(b, c), nh)), _per(σ / 2, _croce(nh, _meno(a, c))))
    end
    lin = ntuple(Val(n)) do k
        Λc === nothing && return z3
        d = _meno(P[mod1(k + 1, n)], P[k])
        return _per(Λc[k] / _norma(d), d)
    end
    # (una variabile riassegnata e catturata da una funzione anonima finirebbe in una scatola,
    # in memoria: per questo le somme sono in funzioni a parte)
    gc = _somma_tri(tri, z3)
    Ltot = _lunghezza_contorno(P)
    G = ntuple(Val(n)) do k
        km = mod1(k - 1, n)
        g = _più(_meno(_più(tri[k][1], tri[km][2]), lin[k]), lin[km])
        return _più(g, _centro_su(P, c, gc, Ltot, k))
    end
    return _appiattisci(G)
end

function _somma_tri(tri, gc)
    for (ga, gb) in tri
        gc = _meno(gc, _più(ga, gb))
    end
    return gc
end

function _volume_e_centro(P::NTuple{n,NTuple{3,T}}, c) where {n,T}
    v = zero(T)
    gc = (zero(T), zero(T), zero(T))
    for k in 1:n
        a, b = P[k], P[mod1(k + 1, n)]
        v += _punto(a, _croce(b, c))
        gc = _più(gc, _croce(a, b))
    end
    return v, gc
end

# Il contributo di volume della faccia, v = Σ P_k·(P_k+1 × c)/6, e il suo gradiente appiattito.
function _volume_faccia(P::NTuple{n,NTuple{3,T}}) where {n,T}
    c = _baricentro(P)
    v, gc = _volume_e_centro(P, c)
    Ltot = _lunghezza_contorno(P)
    G = ntuple(Val(n)) do k
        g = _più(_croce(P[mod1(k + 1, n)], c), _croce(c, P[mod1(k - 1, n)]))
        return _per(1 / 6, _più(g, _centro_su(P, c, gc, Ltot, k)))
    end
    return v / 6, _appiattisci(G)
end

struct _EtichettaHessiana end

# La matrice jacobiana di z ↦ g(punti(z)), in avanti, tutte le derivate in un passaggio.
function _jacobiana(g, P::NTuple{n,NTuple{3,Float64}}) where {n}
    D = ForwardDiff.Dual{_EtichettaHessiana}
    seme(i) = D(P[(i - 1) ÷ 3 + 1][(i - 1) % 3 + 1],
                ForwardDiff.Partials(ntuple(j -> Float64(i == j), Val(3n))))
    Pd = ntuple(k -> (seme(3k - 2), seme(3k - 1), seme(3k)), Val(n))
    r = g(Pd)
    H = Matrix{Float64}(undef, 3n, 3n)
    for i in 1:3n, j in 1:3n
        H[i, j] = ForwardDiff.partials(r[i], j)
    end
    return H
end

_hessiana_area(P, σ, Λc) = _jacobiana(Pd -> _gradiente_area(Pd, σ, Λc), P)
_hessiana_volume(P) = _jacobiana(Pd -> _volume_faccia(Pd)[2], P)

# La regola della catena per le coordinate dei vertici di una faccia traslate di `tr` periodi:
# ognuna dipende dalla coordinata del vertice (`vert`, 0 se è bloccato) e, nei fogli periodici,
# dalla voce della scatola (`scat`) moltiplicata per l'immagine (`coef`, 0 se nessuna).
function _catena(t::Tessuto3D, x, f::Faccia, tr)
    n = length(x)
    vert = Int[]; scat = Int[]; coef = Float64[]
    for (v, im) in zip(f.vertici, f.immagini), a in 1:3
        i = _somma(im, tr)
        push!(vert, t.fissi[v] ? 0 : 3v - 3 + a)
        b = t.periodico && a <= 2 && i[a] != 0
        push!(scat, b ? n - 2 + a : 0)
        push!(coef, b ? Float64(i[a]) : 0.0)
    end
    return (vert, scat, coef)
end

# Aggiunge w·Cᵀ H C (C la regola della catena) alla matrice piena Hc della cellula, che ha per
# incognite `incognite` (ordinate).
function _accumula!(Hc, incognite, H, (vert, scat, coef), w)
    loc(i) = searchsortedfirst(incognite, i)
    for β in axes(H, 2), α in axes(H, 1)
        h = w * H[α, β]
        h == 0 && continue
        for (i, ci) in ((vert[α], 1.0), (scat[α], coef[α])), (j, cj) in ((vert[β], 1.0), (scat[β], coef[β]))
            (i == 0 || j == 0 || ci == 0 || cj == 0) && continue
            Hc[loc(i), loc(j)] += ci * cj * h
        end
    end
end

# Una faccia della cellula, con le posizioni P (una tupla, la cui lunghezza decide quale versione
# compilata dei conti si usa): la tensione in Ht (se `area`), ±∇²v in Hv, ±∇v in G; restituisce ±v.
function _faccia_in_cellula!(Ht, Hv, G, incognite, catena, P, s, area, σ, Λc)
    area && _accumula!(Ht, incognite, _hessiana_area(P, σ, Λc), catena, 1.0)
    _accumula!(Hv, incognite, _hessiana_volume(P), catena, Float64(s))
    v, g = _volume_faccia(P)
    vert, scat, coef = catena
    for α in eachindex(g)
        vert[α] > 0 && (G[searchsortedfirst(incognite, vert[α])] += s * g[α])
        coef[α] != 0 && (G[searchsortedfirst(incognite, scat[α])] += s * coef[α] * g[α])
    end
    return s * v
end

function _hessiana_cellula!(I, J, V, t::Tessuto3D, p::ParametriMeccanici, x, c::Int)
    cel = t.cellule[c]
    catene = [_catena(t, x, t.facce[fi], tr) for (fi, tr) in zip(cel.facce, cel.traslazioni)]
    incognite = Int[]
    for (vert, scat, _) in catene
        append!(incognite, vert, scat)
    end
    incognite = filter!(!=(0), sort!(unique!(incognite)))
    m = length(incognite)
    Ht = zeros(m, m)
    Hv = zeros(m, m)
    G = zeros(m)
    Vc = 0.0
    for (q, (fi, s, tr)) in enumerate(zip(cel.facce, cel.segni, cel.traslazioni))
        f = t.facce[fi]
        P = Tuple(_pos(t, x, v, _somma(im, tr)) for (v, im) in zip(f.vertici, f.immagini))
        # la tensione della faccia la raccoglie la sua prima cellula
        area = first(f.cellule) == c
        Λc = area && f.etichetta == apicale ?
             Tuple(t.intenzioni[c].Λa / f.condivisione[k] for k in eachindex(f.vertici)) : nothing
        Vc += _faccia_in_cellula!(Ht, Hv, G, incognite, catene[q], P, s, area, _tensione_faccia(t, f), Λc)
    end
    w = p.K * (Vc - t.intenzioni[c].V0)
    for b in 1:m, a in 1:m
        h = Ht[a, b] + w * Hv[a, b] + p.K * G[a] * G[b]
        h == 0 && continue
        push!(I, incognite[a]); push!(J, incognite[b]); push!(V, h)
    end
end

"""
    _hessiana_esatta(t, p, x; ancora = nothing, peso = 0.0)

La matrice delle derivate seconde di E(x) + (peso/2)·|x − ancora|² in `x`, sparsa e simmetrica,
esatta (gradienti dei termini locali in forma chiusa, differenziati automaticamente).
"""
function _hessiana_esatta(t::Tessuto3D, p::ParametriMeccanici, x; ancora = nothing, peso = 0.0)
    _senza_contatto(p)
    n = length(x)
    blocchi = collect(Iterators.partition(eachindex(t.cellule), 16))
    pezzi = Vector{Tuple{Vector{Int},Vector{Int},Vector{Float64}}}(undef, length(blocchi))
    Threads.@threads for b in eachindex(blocchi)
        local Ib = Int[]
        local Jb = Int[]
        local Vb = Float64[]
        for c in blocchi[b]
            _hessiana_cellula!(Ib, Jb, Vb, t, p, x, c)
        end
        pezzi[b] = (Ib, Jb, Vb)
    end
    I = reduce(vcat, [pz[1] for pz in pezzi])
    J = reduce(vcat, [pz[2] for pz in pezzi])
    V = reduce(vcat, [pz[3] for pz in pezzi])
    if ancora !== nothing && peso != 0
        append!(I, 1:n); append!(J, 1:n); append!(V, fill(Float64(peso), n))
    end
    H = sparse(I, J, V, n, n)
    return (H + H') / 2
end
