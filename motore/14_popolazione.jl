# Popolazione: la divisione di una cellula (specifica, sezione 6, passo 7; repertorio: "dividersi,
# con un piano di divisione"). La morte è l'apoptosi seguita dall'estrusione di
# 09_estrusione.jl; la passata (13_programma.jl) la guida.
#
# Solo nel guscio a prismi: ogni cellula ha una faccia apicale e una basale con gli stessi
# vertici del piano (apicale 2j − 1, basale 2j) e pareti laterali di quattro vertici. Lì la
# divisione è un taglio del poligono della cellula, come gli scambi di 08_scambi.jl, e il
# tessuto si ricostruisce dai poligoni. Nel modello generale (scutoidi, 11_operazioni.jl)
# divisione e morte non ci sono ancora (docs/roadmap.md, R7).
#
# Il piano di divisione è verticale (contiene l'asse apico-basale) e passa per il centro della
# faccia apicale. La sua traccia nel piano del foglio ha l'angolo dato, oppure, con angolo NaN,
# segue la regola di Hertwig: la cellula si divide perpendicolarmente al proprio asse lungo.
# I punti nuovi stanno sugli spigoli tagliati, nello stesso punto relativo in alto e in basso:
# suddividere uno spigolo non cambia la geometria delle vicine (07_guscio.jl), quindi la
# divisione non muove nessun vertice e, in un foglio piatto, conserva il volume.
# La madre e la figlia hanno ciascuna metà del volume bersaglio della madre.

export dividi!, centro_basale

# Il tessuto è un guscio a prismi?
function _prismatico(t::Tessuto3D)
    for cel in t.cellule
        fa = [t.facce[fi] for fi in cel.facce]
        ap = [f for f in fa if f.etichetta == apicale]
        ba = [f for f in fa if f.etichetta == basale]
        la = [f for f in fa if f.etichetta == laterale]
        (length(ap) == 1 && length(ba) == 1) || return false
        (all(isodd, ap[1].vertici) && all(iseven, ba[1].vertici)) || return false
        Set((v + 1) ÷ 2 for v in ap[1].vertici) == Set(v ÷ 2 for v in ba[1].vertici) || return false
        (length(la) == length(ap[1].vertici) && all(length(f.vertici) == 4 for f in la)) || return false
    end
    return true
end

"""
    dividi!(t, c; angolo)

Divide la cellula `c` del guscio a prismi con un piano verticale per il centro della sua faccia
apicale; la traccia del piano ha angolo `angolo` (rad) rispetto all'asse x, oppure, con `NaN`,
è perpendicolare all'asse lungo della faccia apicale (Hertwig). La madre resta `c`, la figlia è
l'ultima cellula. Se il piano passa per un vertice lo si ruota di 0.01 rad. Restituisce
l'indice della figlia, o `nothing` se il tessuto non è a prismi o se un vertice della cellula è
bloccato.
"""
function dividi!(t::Tessuto3D, c::Int; angolo)
    _prismatico(t) || return nothing
    P = poligoni(t)
    w = P[c]
    m = length(w)
    any(t.fissi[2j - 1] for (j, _) in w) && return nothing
    A = [_pos(t, t.x, 2j - 1, im) for (j, im) in w]
    B = [_pos(t, t.x, 2j, im) for (j, im) in w]
    g = _baricentro(A)
    # direzione della traccia del piano nel piano del foglio
    if isnan(angolo)
        sxx = sum((p[1] - g[1])^2 for p in A); syy = sum((p[2] - g[2])^2 for p in A)
        sxy = sum((p[1] - g[1]) * (p[2] - g[2]) for p in A)
        θ = atan(2sxy, sxx - syy) / 2                 # l'asse lungo
        d = (-sin(θ), cos(θ))                         # la traccia gli è perpendicolare
    else
        d = (cos(angolo), sin(angolo))
    end
    # se la traccia passa per un vertice (nel reticolo regolare succede: la regola di Hertwig su
    # un esagono allungato lo taglia da vertice a vertice) la si ruota di 0.01 rad
    distanze(d) = [(p[1] - g[1]) * d[2] - (p[2] - g[2]) * d[1] for p in A]
    s = distanze(d)
    if any(x -> abs(x) < 1e-6 * maximum(abs, s), s)
        d = (d[1] * cos(0.01) - d[2] * sin(0.01), d[1] * sin(0.01) + d[2] * cos(0.01))
        s = distanze(d)
        any(x -> abs(x) < 1e-6 * maximum(abs, s), s) && return nothing
    end
    tagli = [k for k in 1:m if s[k] * s[mod1(k + 1, m)] < 0]
    length(tagli) == 2 || return nothing

    Lx, Ly = _periodo(t, t.x)
    nv = _nvert(t) ÷ 2                               # vertici del piano
    coda = t.periodico ? t.x[end-1:end] : Float64[]
    x = t.periodico ? t.x[1:end-2] : copy(t.x)
    ancore = t.periodico ? t.ancore[1:end-2] : copy(t.ancore)
    fissi = copy(t.fissi)
    nuovi = Tuple{Int,NTuple{2,Int}}[]
    for k in tagli
        k2 = mod1(k + 1, m)
        τ = s[k] / (s[k] - s[k2])
        im = w[k][2]
        for (Q, _) in ((A, :apicale), (B, :basale))
            p = Q[k] .+ τ .* (Q[k2] .- Q[k])
            q = (p[1] - im[1] * Lx, p[2] - im[2] * Ly, p[3])
            append!(x, q); append!(ancore, q); push!(fissi, false)
        end
        nv += 1
        push!(nuovi, (nv, im))
    end
    append!(x, coda); append!(ancore, t.periodico ? t.ancore[end-1:end] : Float64[])

    # le vicine che hanno uno spigolo tagliato prendono il punto nuovo sullo spigolo
    for (r, k) in enumerate(tagli)
        k2 = mod1(k + 1, m)
        (ja, ia), (jb, ib) = w[k], w[k2]
        for (ci, u) in enumerate(P)
            ci == c && continue
            L = length(u)
            for h in 1:L
                (jx, ix), (jy, iy) = u[h], u[mod1(h + 1, L)]
                if jx == jb && jy == ja && _im_meno(ix, iy) == _im_meno(ib, ia)
                    insert!(u, h + 1, (nuovi[r][1], iy))      # immagine come quella di a in u
                    break
                end
            end
        end
    end
    # madre e figlia
    k1, k2 = tagli
    giro(da, a) = [w[mod1(h, m)] for h in da:(a < da ? a + m : a)]
    madre = vcat([nuovi[1]], giro(k1 + 1, k2), [nuovi[2]])
    figlia = vcat([nuovi[2]], giro(k2 + 1, k1), [nuovi[1]])
    P[c] = madre
    push!(P, figlia)
    r = t.intenzioni[c]
    meta = Intenzioni3D(r.σa, r.σb, r.σl, r.Λa, r.V0 / 2)
    intenzioni = copy(t.intenzioni); intenzioni[c] = meta; push!(intenzioni, meta)
    nuovo = tessuto_da_poligoni(intenzioni, P, x; periodico = t.periodico, fissi = BitVector(fissi), ancore)
    t.x = nuovo.x; t.facce = nuovo.facce; t.cellule = nuovo.cellule
    t.intenzioni = nuovo.intenzioni; t.fissi = nuovo.fissi; t.ancore = nuovo.ancore
    return length(t.cellule)
end

"""
    centro_basale(t, c)

Dove la cellula `c` depone e legge la matrice: il centro della sua faccia basale (µm), riportato
nella scatola nei fogli periodici.
"""
function centro_basale(t::Tessuto3D, c::Int)
    cel = t.cellule[c]
    k = findfirst(fi -> t.facce[fi].etichetta == basale, cel.facce)
    g = _baricentro(_vertici(t, t.x, t.facce[cel.facce[k]], cel.traslazioni[k]))
    t.periodico || return g
    Lx, Ly = _periodo(t, t.x)
    riporta(x, L) = (y = mod(x, L); y == L ? 0.0 : y)
    return (riporta(g[1], Lx), riporta(g[2], Ly), g[3])
end
