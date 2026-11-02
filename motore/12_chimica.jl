# La chimica (specifica, decisione 5; docs/roadmap.md, passo R3): una sostanza che diffonde su
# una griglia cartesiana 3D di voxel, con degradazione del primo ordine, e le cellule che la
# secernono o la assorbono alla propria posizione.
#
#   ∂c/∂t = D ∇²c − k c + sorgenti − assorbimenti
#
# Unità: µm, min; la concentrazione è in quantità/µm³ con la quantità libera (molecole, o una
# unità arbitraria). La griglia ha un passo per asse; ogni asse ha il bordo chiuso (flusso nullo)
# o periodico (per i fogli periodici: il periodo è la scatola).
#
# Un passo di chimica, dt, fa nell'ordine:
#   1. diffusione, schema esplicito (differenze centrate); conserva la massa esattamente ed è
#      stabile per dt ≤ passo_massimo(campo) = 1 / (2 D Σ 1/h²);
#   2. degradazione e secrezione insieme, esatte nel passo: c ← c · exp(−k dt) + s (1 − exp(−k dt))/k,
#      dove s è la sorgente del voxel; la secrezione q di una cellula va sui voxel vicini alla
#      posizione con pesi trilineari (il baricentro di ciò che si aggiunge è la posizione). Così
#      allo stato stazionario ciò che entra esce per degradazione, senza errore di passo;
#   4. assorbimento del primo ordine, esatto nel voxel: c ← c · exp(−κ w dt / V), dove κ (µm³/min)
#      è il volume "ripulito" al minuto e w il peso del voxel; la concentrazione resta positiva.
# `avanza!` divide un intervallo in passi uguali non più lunghi di `passo_massimo`: la chimica
# ha così un passo proprio, più corto di quello meccanico (specifica, sezione 6).

export Bordo, bordo_chiuso, bordo_periodico, Sostanza, Campo, Scambio
export centro_voxel, massa, concentrazione, passo_massimo, avanza!, centro_cellula

@enum Bordo bordo_chiuso bordo_periodico

"""
    Sostanza(nome, D, k)

Una sostanza che diffonde: coefficiente di diffusione `D` (µm²/min) e costante di degradazione
`k` (1/min). La lunghezza di decadimento da una sorgente è √(D/k).
"""
struct Sostanza
    nome::Symbol
    D::Float64
    k::Float64
end

"""
    Campo(sostanza; origine, h, n, bordi = (bordo_chiuso, bordo_chiuso, bordo_chiuso))

Concentrazione della sostanza su `n` = (nx, ny, nz) voxel di lati `h` (µm, uno per asse o uno
solo), a partire dall'angolo `origine`. Il voxel (i, j, l) ha centro
origine + ((i − ½) hx, (j − ½) hy, (l − ½) hz). Le concentrazioni sono in `campo.c`, a zero.
"""
mutable struct Campo
    sostanza::Sostanza
    origine::NTuple{3,Float64}
    h::NTuple{3,Float64}
    c::Array{Float64,3}
    bordi::NTuple{3,Bordo}
end

function Campo(sostanza::Sostanza; origine, h, n, bordi = (bordo_chiuso, bordo_chiuso, bordo_chiuso))
    hh = h isa Number ? (Float64(h), Float64(h), Float64(h)) : Float64.(Tuple(h))
    all(>(0), hh) || throw(ArgumentError("il passo della griglia deve essere positivo"))
    return Campo(sostanza, Float64.(Tuple(origine)), hh, zeros(Tuple(n)...), Tuple(bordi))
end

"""
    Scambio(posizione, secrezione, assorbimento)

Ciò che una cellula fa alla sostanza nella sua posizione (µm): `secrezione` in quantità/min,
`assorbimento` come volume ripulito al minuto (µm³/min).
"""
struct Scambio
    posizione::NTuple{3,Float64}
    secrezione::Float64
    assorbimento::Float64
end

centro_voxel(campo::Campo, i, j, l) =
    ntuple(a -> campo.origine[a] + ((i, j, l)[a] - 0.5) * campo.h[a], 3)

_volume_voxel(campo::Campo) = campo.h[1] * campo.h[2] * campo.h[3]

massa(campo::Campo) = sum(campo.c) * _volume_voxel(campo)

"""
    passo_massimo(campo)

Il passo di chimica più lungo per cui lo schema esplicito resta stabile (min); `Inf` senza
diffusione.
"""
function passo_massimo(campo::Campo)
    D = campo.sostanza.D
    s = sum(1 / campo.h[a]^2 for a in 1:3 if size(campo.c, a) > 1; init = 0.0)
    (D == 0 || s == 0) && return Inf
    return 1 / (2D * s)
end

# Pesi lineari lungo un asse: (indice, peso), uno o due voxel.
function _pesi_asse(campo::Campo, a, x)
    n = size(campo.c, a)
    n == 1 && return ((1, 1.0),)
    s = (x - campo.origine[a]) / campo.h[a] + 0.5          # il centro del voxel i sta in s = i
    i0 = floor(Int, s); f = s - i0
    if campo.bordi[a] == bordo_periodico
        return ((mod1(i0, n), 1 - f), (mod1(i0 + 1, n), f))
    end
    s <= 1 && return ((1, 1.0),)
    s >= n && return ((n, 1.0),)
    return ((i0, 1 - f), (i0 + 1, f))
end

# Pesi trilineari di una posizione: [(CartesianIndex, peso)], con somma 1.
function _pesi(campo::Campo, p)
    px, py, pz = (_pesi_asse(campo, a, p[a]) for a in 1:3)
    return [(CartesianIndex(i, j, l), wi * wj * wl) for (i, wi) in px for (j, wj) in py
            for (l, wl) in pz if wi * wj * wl != 0]
end

"""
    concentrazione(campo, p)

La concentrazione alla posizione `p` (µm), per interpolazione trilineare fra i centri dei voxel;
esatta su un campo lineare. Fuori dalla griglia, sugli assi chiusi, vale quella del voxel di
bordo; sugli assi periodici la posizione si riporta nel periodo.
"""
concentrazione(campo::Campo, p) = sum(w * campo.c[I] for (I, w) in _pesi(campo, p))

# Un passo di diffusione esplicito: c ← c + dt D ∇²c, con il bordo chiuso (il vicino che manca
# è il voxel stesso, flusso nullo) o periodico.
function _diffondi!(campo::Campo, nuovo, dt)
    c = campo.c
    D = campo.sostanza.D
    n = size(c)
    copyto!(nuovo, c)
    for a in 1:3
        n[a] == 1 && continue
        α = dt * D / campo.h[a]^2
        periodico = campo.bordi[a] == bordo_periodico
        e = CartesianIndex(ntuple(b -> b == a ? 1 : 0, 3))
        for I in CartesianIndices(c)
            i = I[a]
            Ip = i < n[a] ? I + e : (periodico ? I - (n[a] - 1) * e : I)
            Im = i > 1 ? I - e : (periodico ? I + (n[a] - 1) * e : I)
            nuovo[I] += α * (c[Ip] - 2c[I] + c[Im])
        end
    end
    copyto!(c, nuovo)
    return campo
end

"""
    avanza!(campo, Δt, scambi = Scambio[]; passo = passo_massimo(campo))

Fa avanzare la chimica di `Δt` minuti, in passi uguali non più lunghi di `passo` (che non può
superare `passo_massimo`). Restituisce `(secreto, assorbito, passi)`: le quantità entrate e
uscite per le cellule. `secreto` è q·Δt: ciò che le cellule hanno emesso, anche la parte che si
è degradata nello stesso passo.
"""
function avanza!(campo::Campo, Δt, scambi = Scambio[]; passo = passo_massimo(campo))
    passo > passo_massimo(campo) * (1 + 1e-12) &&
        throw(ArgumentError("passo $passo oltre il limite di stabilità $(passo_massimo(campo))"))
    passi = isfinite(passo) ? max(1, ceil(Int, Δt / passo - 1e-9)) : 1
    dt = Δt / passi
    V = _volume_voxel(campo)
    pesi = [_pesi(campo, s.posizione) for s in scambi]
    nuovo = similar(campo.c)
    k = campo.sostanza.k
    fattore = exp(-k * dt)
    dt_sorgente = k == 0 ? dt : -expm1(-k * dt) / k     # ∫ exp(−k s) ds sul passo
    secreto = 0.0; assorbito = 0.0
    c = campo.c
    for _ in 1:passi
        campo.sostanza.D > 0 && _diffondi!(campo, nuovo, dt)
        fattore != 1 && (c .*= fattore)
        for (s, P) in zip(scambi, pesi)
            s.secrezione == 0 && continue
            for (I, w) in P
                c[I] += s.secrezione * dt_sorgente * w / V
            end
            secreto += s.secrezione * dt
        end
        for (s, P) in zip(scambi, pesi)
            s.assorbimento == 0 && continue
            for (I, w) in P
                prima = c[I]
                c[I] *= exp(-s.assorbimento * w * dt / V)
                assorbito += (prima - c[I]) * V
            end
        end
    end
    return (secreto = secreto, assorbito = assorbito, passi = passi)
end

"""
    centro_cellula(t, c)

La posizione da cui la cellula `c` del guscio secerne e legge: la media dei suoi vertici (µm),
riportata nella scatola nei fogli periodici.
"""
function centro_cellula(t::Tessuto3D, c::Int)
    cel = t.cellule[c]
    visti = Set{Tuple{Int,NTuple{2,Int}}}()
    S = (0.0, 0.0, 0.0)
    for (fi, tr) in zip(cel.facce, cel.traslazioni)
        f = t.facce[fi]
        for (v, im) in zip(f.vertici, f.immagini)
            chiave = (v, _somma(im, tr))
            chiave in visti && continue
            push!(visti, chiave)
            S = S .+ _pos(t, t.x, v, chiave[2])
        end
    end
    g = S ./ length(visti)
    t.periodico || return g
    Lx, Ly = _periodo(t, t.x)
    riporta(x, L) = (y = mod(x, L); y == L ? 0.0 : y)     # mod di −1e-16 dà L, arrotondato
    return (riporta(g[1], Lx), riporta(g[2], Ly), g[3])
end
