# Il guscio 3D a due superfici (specifica, sezione 3).
#
# Una maglia apicale e la sua copia basale; ogni cellula è il prisma fra un poligono apicale e
# quello basale corrispondente, con facce laterali quadrilatere che li collegano. Le facce
# laterali sono condivise fra cellule vicine. Le facce non sono necessariamente piane: area e
# volume si calcolano a ventaglio dal baricentro dei vertici della faccia, la stessa
# triangolazione per le due cellule che la condividono, così i volumi si incastrano.
#
#   E = Σ_facce σ_f · A_f  +  Σ_spigoli apicali Λ · L  +  Σ_cellule (K/2) · (V_c − V0_c)²
#
# σ: tensione di superficie per faccia (nN/µm), con valori distinti per apicale, basale e
# laterale; Λ: tensione di linea sugli spigoli apicali (nN), l'anello apicale che `costringi`
# alza; K: incomprimibilità (nN/µm⁵), penalità sul volume. Qui K di `ParametriMeccanici` ha
# queste unità, non quelle d'area del modello laterale. Uno spigolo o una faccia condivisi
# prendono la media delle tensioni delle cellule, come nel modello laterale. Il contatto fra
# parti lontane non c'è ancora in 3D: con k_contatto > 0 il motore si rifiuta.
#
# Le posizioni stanno in x = [x₁, y₁, z₁, x₂, …]; il vertice apicale j è 2j-1, il basale 2j.
#
# Foglio periodico: le cellule al bordo destro toccano quelle al bordo sinistro, e così sopra
# e sotto, come su un toro piatto. Un vertice di una faccia è allora un vertice più un numero
# intero di periodi (la sua "immagine"), e le due dimensioni del periodo Lx, Ly sono gli ultimi
# due elementi di x: sono incognite come le posizioni, così il foglio trova da solo la sua
# estensione di equilibrio, senza bordi liberi.

export Faccia, Cellula3D, Intenzioni3D, Tessuto3D
export foglio_esagonale, foglio_periodico, area_faccia, volume_cellula
export area_apicale, altezza_prisma

"""
    Faccia(vertici, immagini, etichetta, cellule, segni, condivisione)

Poligono di vertici; `immagini[k]` è il numero di periodi (ix, iy) da aggiungere al vertice k.
`segni[c]` è +1 se l'ordine dei vertici gira in senso antiorario visto da fuori della cellula
`cellule[c]`. Per una faccia apicale, `condivisione[k]` è il numero di cellule che hanno lo
spigolo k → k+1 (serve alla tensione di linea).
"""
struct Faccia
    vertici::Vector{Int}
    immagini::Vector{NTuple{2,Int}}
    etichetta::Etichetta
    cellule::Vector{Int}
    segni::Vector{Int}
    condivisione::Vector{Int}
end

"""
    Cellula3D(facce, segni, traslazioni)

Le facce della cellula, con il segno dell'orientamento uscente e la traslazione (in periodi)
che porta la faccia nel riferimento della cellula: il volume si calcola con tutte le facce
nello stesso riferimento.
"""
struct Cellula3D
    facce::Vector{Int}
    segni::Vector{Int}
    traslazioni::Vector{NTuple{2,Int}}
end

"""
    Intenzioni3D(σa, σb, σl, Λa, V0)

Tensioni di superficie apicale, basale e laterale (nN/µm), tensione di linea dell'anello
apicale (nN) e volume bersaglio (µm³).
"""
struct Intenzioni3D
    σa::Float64
    σb::Float64
    σl::Float64
    Λa::Float64
    V0::Float64
end

mutable struct Tessuto3D
    x::Vector{Float64}
    facce::Vector{Faccia}
    cellule::Vector{Cellula3D}
    intenzioni::Vector{Intenzioni3D}
    periodico::Bool
end

_nvert(t::Tessuto3D) = (length(t.x) - (t.periodico ? 2 : 0)) ÷ 3
_periodo(t::Tessuto3D, x) = t.periodico ? (x[end-1], x[end]) : (0.0, 0.0)

function _pos(t::Tessuto3D, x, v, (ix, iy))
    Lx, Ly = _periodo(t, x)
    return (x[3v-2] + ix * Lx, x[3v-1] + iy * Ly, x[3v])
end

_somma(a, b) = (a[1] + b[1], a[2] + b[2])
_meno(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
_croce(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
_punto(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_norma(a) = sqrt(_punto(a, a))

# Posizioni dei vertici di una faccia, traslate di `tr` periodi.
_vertici(t::Tessuto3D, x, f::Faccia, tr = (0, 0)) =
    [_pos(t, x, v, _somma(im, tr)) for (v, im) in zip(f.vertici, f.immagini)]

function _baricentro(P)
    n = length(P)
    return (sum(p[1] for p in P) / n, sum(p[2] for p in P) / n, sum(p[3] for p in P) / n)
end

# Aggiunge w·G al gradiente del vertice v nell'immagine im (e del periodo, se periodico).
function _aggiungi!(g, t::Tessuto3D, v, im, G, w)
    g[3v-2] += w * G[1]; g[3v-1] += w * G[2]; g[3v] += w * G[3]
    if t.periodico
        g[end-1] += w * im[1] * G[1]
        g[end] += w * im[2] * G[2]
    end
    return g
end

"Area di una faccia (µm²), a ventaglio dal baricentro dei suoi vertici."
area_faccia(t::Tessuto3D, f::Faccia, x = t.x) = _area(_vertici(t, x, f))

function _area(P)
    c = _baricentro(P)
    n = length(P)
    return sum(_norma(_croce(_meno(P[k], c), _meno(P[mod1(k + 1, n)], c))) for k in 1:n) / 2
end

"Volume della cellula `c` (µm³), dal teorema della divergenza sulla stessa triangolazione."
function volume_cellula(t::Tessuto3D, c::Int, x = t.x)
    V = 0.0
    cel = t.cellule[c]
    for (fi, s, tr) in zip(cel.facce, cel.segni, cel.traslazioni)
        P = _vertici(t, x, t.facce[fi], tr)
        g = _baricentro(P)
        n = length(P)
        for k in 1:n
            V += s * _punto(P[k], _croce(P[mod1(k + 1, n)], g)) / 6
        end
    end
    return V
end

function _tensione_faccia(t::Tessuto3D, f::Faccia)
    σ = 0.0
    for c in f.cellule
        int = t.intenzioni[c]
        σ += f.etichetta == apicale ? int.σa : f.etichetta == basale ? int.σb : int.σl
    end
    return σ / length(f.cellule)
end

function _senza_contatto(p::ParametriMeccanici)
    p.k_contatto == 0 ||
        throw(ArgumentError("il contatto fra parti lontane non è ancora implementato in 3D"))
end

function energia(t::Tessuto3D, p::ParametriMeccanici, x = t.x)
    _senza_contatto(p)
    E = 0.0
    for f in t.facce
        P = _vertici(t, x, f)
        E += _tensione_faccia(t, f) * _area(P)
        if f.etichetta == apicale
            # anello apicale: ogni spigolo pesa per la quota di questa cellula
            Λ = t.intenzioni[only(f.cellule)].Λa
            n = length(P)
            for k in 1:n
                E += Λ / f.condivisione[k] * _norma(_meno(P[mod1(k + 1, n)], P[k]))
            end
        end
    end
    for c in eachindex(t.cellule)
        E += p.K / 2 * (volume_cellula(t, c, x) - t.intenzioni[c].V0)^2
    end
    return E
end

function gradiente!(g, t::Tessuto3D, p::ParametriMeccanici, x = t.x)
    _senza_contatto(p)
    fill!(g, 0)
    # tensione di superficie e anello apicale, faccia per faccia
    for f in t.facce
        σ = _tensione_faccia(t, f)
        P = _vertici(t, x, f)
        n = length(P)
        c = _baricentro(P)
        gc = (0.0, 0.0, 0.0)
        for k in 1:n
            k2 = mod1(k + 1, n)
            a, b = P[k], P[k2]
            # triangolo (a, b, c): ∂A/∂a = ½ (b−c)×n̂, ∂A/∂b = ½ n̂×(a−c), ∂A/∂c = −(somma)
            m = _croce(_meno(a, c), _meno(b, c))
            L = _norma(m)
            if L > 0
                nh = (m[1] / L, m[2] / L, m[3] / L)
                ga = _croce(_meno(b, c), nh)
                gb = _croce(nh, _meno(a, c))
                _aggiungi!(g, t, f.vertici[k], f.immagini[k], ga, σ / 2)
                _aggiungi!(g, t, f.vertici[k2], f.immagini[k2], gb, σ / 2)
                gc = (gc[1] - σ / 2 * (ga[1] + gb[1]), gc[2] - σ / 2 * (ga[2] + gb[2]),
                      gc[3] - σ / 2 * (ga[3] + gb[3]))
            end
            if f.etichetta == apicale
                d = _meno(b, a)
                w = t.intenzioni[only(f.cellule)].Λa / f.condivisione[k] / _norma(d)
                _aggiungi!(g, t, f.vertici[k2], f.immagini[k2], d, w)
                _aggiungi!(g, t, f.vertici[k], f.immagini[k], d, -w)
            end
        end
        for k in 1:n
            _aggiungi!(g, t, f.vertici[k], f.immagini[k], gc, 1 / n)
        end
    end
    # volume, cellula per cellula, con le facce nel riferimento della cellula
    for (ci, cel) in enumerate(t.cellule)
        Pc = p.K * (volume_cellula(t, ci, x) - t.intenzioni[ci].V0)
        for (fi, s, tr) in zip(cel.facce, cel.segni, cel.traslazioni)
            f = t.facce[fi]
            P = _vertici(t, x, f, tr)
            n = length(P)
            c = _baricentro(P)
            w = s * Pc / 6
            gc = (0.0, 0.0, 0.0)
            for k in 1:n
                k2 = mod1(k + 1, n)
                a, b = P[k], P[k2]
                # ∂(a·(b×c))/∂a = b×c, /∂b = c×a, /∂c = a×b
                _aggiungi!(g, t, f.vertici[k], _somma(f.immagini[k], tr), _croce(b, c), w)
                _aggiungi!(g, t, f.vertici[k2], _somma(f.immagini[k2], tr), _croce(c, a), w)
                ab = _croce(a, b)
                gc = (gc[1] + ab[1], gc[2] + ab[2], gc[3] + ab[3])
            end
            for k in 1:n
                _aggiungi!(g, t, f.vertici[k], _somma(f.immagini[k], tr), gc, w / n)
            end
        end
    end
    return g
end

"Forze sui vertici (e sulle due dimensioni del periodo, se periodico), −∇E."
forze(t::Tessuto3D, p::ParametriMeccanici, x = t.x) = -gradiente!(similar(x), t, p, x)

"Il modulo massimo fra le forze sui vertici e quelle sulle due dimensioni del periodo."
function forza_per_vertice(t::Tessuto3D, g)
    m = maximum(sqrt(g[3k-2]^2 + g[3k-1]^2 + g[3k]^2) for k in 1:_nvert(t))
    return t.periodico ? max(m, abs(g[end-1]), abs(g[end])) : m
end

function controlla_ammissibile(t::Tessuto3D, p::ParametriMeccanici, x = t.x)
    for f in t.facce
        P = _vertici(t, x, f)
        n = length(P)
        for k in 1:n
            L = _norma(_meno(P[mod1(k + 1, n)], P[k]))
            L < p.lunghezza_minima && throw(ConfigurazioneNonAmmissibile(
                "spigolo $(f.vertici[k])–$(f.vertici[mod1(k + 1, n)]) lungo $L µm"))
        end
    end
    for c in eachindex(t.cellule)
        V = volume_cellula(t, c, x)
        V <= 0 && throw(ConfigurazioneNonAmmissibile("cellula $c con volume $V µm³"))
    end
    return nothing
end

_controlla_contatto(t::Tessuto3D, p::ParametriMeccanici) = (_senza_contatto(p); nothing)

# Costruzione di un foglio di prismi esagonali regolari dai centri delle cellule.
# `avvolgi(px, py)` riporta un punto nella scatola e dice di quanti periodi lo si è spostato.
function _foglio(intenzioni, centri, lato, altezza, avvolgi, periodico, Lx, Ly)
    length(intenzioni) == length(centri) ||
        throw(ArgumentError("servono $(length(centri)) intenzioni, non $(length(intenzioni))"))
    indice = Dict{Tuple{Int,Int},Int}()
    piani = Tuple{Float64,Float64}[]
    function vertice(px, py)
        (qx, qy), im = avvolgi(px, py)
        chiave = (round(Int, qx * 1e6), round(Int, qy * 1e6))
        j = get!(indice, chiave) do
            push!(piani, (qx, qy))
            length(piani)
        end
        return j, im
    end
    ap(j) = 2j - 1
    ba(j) = 2j
    facce = Faccia[]
    cellule = Cellula3D[]
    laterali = Dict{Tuple{Int,Int,Int,Int},Tuple{Int,NTuple{2,Int}}}()  # → faccia, immagine di riferimento
    conta_apicali = Dict{Tuple{Int,Int,Int,Int},Int}()
    angoli = [[vertice(c[1] + lato * cos(π / 6 + k * π / 3), c[2] + lato * sin(π / 6 + k * π / 3))
               for k in 0:5] for c in centri]
    chiave_spigolo(a, b) = begin
        (ja, ia), (jb, ib) = a, b
        # lo stesso spigolo visto da due cellule differisce per una traslazione: la chiave
        # usa gli indici e lo scarto fra le due immagini
        ja < jb ? (ja, jb, ib[1] - ia[1], ib[2] - ia[2]) : (jb, ja, ia[1] - ib[1], ia[2] - ib[2])
    end
    for w in angoli, k in 1:6
        ch = chiave_spigolo(w[k], w[mod1(k + 1, 6)])
        conta_apicali[ch] = get(conta_apicali, ch, 0) + 1
    end
    for (ci, w) in enumerate(angoli)
        fs, ss, ts = Int[], Int[], NTuple{2,Int}[]
        cond = [conta_apicali[chiave_spigolo(w[k], w[mod1(k + 1, 6)])] for k in 1:6]
        push!(facce, Faccia([ap(j) for (j, _) in w], [im for (_, im) in w], apicale, [ci], [1], cond))
        push!(fs, length(facce)); push!(ss, 1); push!(ts, (0, 0))
        wr = reverse(w)
        push!(facce, Faccia([ba(j) for (j, _) in wr], [im for (_, im) in wr], basale, [ci], [1], Int[]))
        push!(fs, length(facce)); push!(ss, 1); push!(ts, (0, 0))
        for k in 1:6
            (jp, ip), (jq, iq) = w[k], w[mod1(k + 1, 6)]
            ch = chiave_spigolo(w[k], w[mod1(k + 1, 6)])
            if haskey(laterali, ch)
                fi, rif = laterali[ch]
                push!(facce[fi].cellule, ci); push!(facce[fi].segni, -1)
                # la faccia è salvata nel riferimento della prima cellula: traslazione per questa
                vq = facce[fi].vertici[1] == ba(jq) ? iq : ip
                push!(fs, fi); push!(ss, -1); push!(ts, (vq[1] - rif[1], vq[2] - rif[2]))
            else
                push!(facce, Faccia([ba(jp), ba(jq), ap(jq), ap(jp)], [ip, iq, iq, ip], laterale,
                                    [ci], [1], Int[]))
                laterali[ch] = (length(facce), ip)
                push!(fs, length(facce)); push!(ss, 1); push!(ts, (0, 0))
            end
        end
        push!(cellule, Cellula3D(fs, ss, ts))
    end
    x = zeros(3 * 2 * length(piani) + (periodico ? 2 : 0))
    for (j, (px, py)) in enumerate(piani)
        x[3ap(j)-2], x[3ap(j)-1], x[3ap(j)] = px, py, altezza
        x[3ba(j)-2], x[3ba(j)-1], x[3ba(j)] = px, py, 0.0
    end
    periodico && (x[end-1] = Lx; x[end] = Ly)
    return Tessuto3D(x, facce, cellule, copy(intenzioni), periodico)
end

"""
    foglio_esagonale(intenzioni, R; lato, altezza)

Pezzo di foglio piatto di prismi esagonali regolari: la cellula centrale e `R` anelli attorno
(1 + 3R(R+1) cellule, tante quante `intenzioni`), con bordo libero. Basale su z = 0, apicale su
z = `altezza`. Con R = 0 è un prisma isolato.
"""
function foglio_esagonale(intenzioni::Vector{Intenzioni3D}, R::Int; lato, altezza)
    centri = Tuple{Float64,Float64}[]
    for q in -R:R, r in -R:R
        abs(q + r) <= R || continue
        push!(centri, (lato * sqrt(3) * (q + r / 2), lato * 1.5 * r))
    end
    sort!(centri; by = c -> (round(hypot(c...); digits = 6), atan(c[2], c[1])))
    return _foglio(intenzioni, centri, lato, altezza, (px, py) -> ((px, py), (0, 0)), false, 0.0, 0.0)
end

"""
    foglio_periodico(intenzioni, nx, ny; lato, altezza)

Foglio piatto periodico di `nx` × `ny` prismi esagonali (ny pari), senza bordi. Le dimensioni
del periodo partono da Lx = nx·√3·lato, Ly = ny·1.5·lato e sono libere di cambiare. Le cellule
sono numerate per righe.
"""
function foglio_periodico(intenzioni::Vector{Intenzioni3D}, nx::Int, ny::Int; lato, altezza)
    iseven(ny) || throw(ArgumentError("ny deve essere pari per chiudere il reticolo esagonale"))
    (nx >= 3 && ny >= 4) || throw(ArgumentError("servono almeno 3 × 4 cellule"))
    Lx, Ly = nx * sqrt(3) * lato, ny * 1.5 * lato
    centri = [(sqrt(3) * lato * (q + (isodd(r) ? 0.5 : 0.0)), 1.5 * lato * r)
              for r in 0:ny-1 for q in 0:nx-1]
    function avvolgi(px, py)
        ix, iy = floor(Int, px / Lx + 1e-9), floor(Int, py / Ly + 1e-9)
        return (px - ix * Lx, py - iy * Ly), (ix, iy)
    end
    return _foglio(intenzioni, centri, lato, altezza, avvolgi, true, Lx, Ly)
end

# Misure per l'osservatore.

_faccia_di(t::Tessuto3D, c::Int, e::Etichetta) =
    t.facce[only(fi for fi in t.cellule[c].facce if t.facce[fi].etichetta == e)]

"Area della faccia apicale della cellula `c` (µm²)."
area_apicale(t::Tessuto3D, c::Int, x = t.x) = _area(_vertici(t, x, _faccia_di(t, c, apicale)))

"Distanza fra i baricentri delle facce apicale e basale della cellula `c` (µm)."
function altezza_prisma(t::Tessuto3D, c::Int, x = t.x)
    i = findfirst(fi -> t.facce[fi].etichetta == apicale, t.cellule[c].facce)
    j = findfirst(fi -> t.facce[fi].etichetta == basale, t.cellule[c].facce)
    cel = t.cellule[c]
    Pa = _baricentro(_vertici(t, x, t.facce[cel.facce[i]], cel.traslazioni[i]))
    Pb = _baricentro(_vertici(t, x, t.facce[cel.facce[j]], cel.traslazioni[j]))
    return _norma(_meno(Pa, Pb))
end
