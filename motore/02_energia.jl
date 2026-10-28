# Il funzionale di energia e le forze (specifica, decisione 1).
#
#   E = Σ_spigoli Γ_s · L_s  +  Σ_cellule (K/2) · (A_c − A0_c)²  +  E_contatto
#
# Tensione di linea su ogni spigolo, incomprimibilità come penalità sull'area, esclusione di
# volume fra parti lontane del tessuto (05_contatto.jl). Le forze sono
# −∇E, scritte a mano in forma esatta: appaiate per costruzione, perché ogni termine dipende
# solo da differenze di posizioni. La batteria le confronta con le differenze finite.
# Unità: E in nN·µm, forze in nN, K in nN/µm³.

export ParametriMeccanici, tensione, energia, gradiente!, forze

"""
    ParametriMeccanici(; K, η, lunghezza_minima, k_contatto, d_contatto)

- `K`: incomprimibilità (nN/µm³), penalità sull'area
- `η`: attrito della dinamica sovrasmorzata (nN·min/µm)
- `lunghezza_minima`: sotto questa lunghezza uno spigolo rende la configurazione non
  ammissibile (µm)
- `k_contatto`: rigidezza del contatto fra parti lontane (nN/µm); 0 spegne il meccanismo
- `d_contatto`: distanza sotto cui vertice e spigolo lontani si respingono (µm)

Nessun valore ha un default: ogni corsa dichiara tutti i parametri, compreso se il contatto
è acceso.
"""
Base.@kwdef struct ParametriMeccanici
    K::Float64
    η::Float64
    lunghezza_minima::Float64
    k_contatto::Float64
    d_contatto::Float64
end

function _tensione_cellula(int::Intenzioni, e::Etichetta)
    e == apicale && return int.Γa
    e == basale && return int.Γb
    return int.Γl
end

"""
    tensione(t, s)

Tensione dello spigolo `s`: la media delle tensioni che le cellule adiacenti dichiarano per
quell'etichetta. Uno spigolo laterale condiviso pesa una volta sola, non due.
"""
function tensione(t::Tessuto, s::Spigolo)
    somma = 0.0
    for c in s.cellule
        somma += _tensione_cellula(t.intenzioni[c], s.etichetta)
    end
    return somma / length(s.cellule)
end

function energia(t::Tessuto, p::ParametriMeccanici, x = t.x)
    E = zero(eltype(x))
    for s in t.spigoli
        E += tensione(t, s) * lunghezza(x, s)
    end
    for (c, cel) in enumerate(t.cellule)
        dA = area(x, cel) - t.intenzioni[c].A0
        E += p.K / 2 * dA^2
    end
    return E + energia_contatto(t, p, x)
end

"""
    gradiente!(g, t, p, x = t.x)

Scrive ∇E in `g`.
"""
function gradiente!(g, t::Tessuto, p::ParametriMeccanici, x = t.x)
    fill!(g, 0)
    for s in t.spigoli
        xi, yi = posizione(x, s.i)
        xj, yj = posizione(x, s.j)
        L = hypot(xj - xi, yj - yi)
        Γ = tensione(t, s)
        ux, uy = Γ * (xj - xi) / L, Γ * (yj - yi) / L
        g[2s.j-1] += ux; g[2s.j] += uy
        g[2s.i-1] -= ux; g[2s.i] -= uy
    end
    for (c, cel) in enumerate(t.cellule)
        P = p.K * (area(x, cel) - t.intenzioni[c].A0)
        v = cel.vertici
        n = length(v)
        for k in 1:n
            xp, yp = posizione(x, v[mod1(k - 1, n)])
            xn, yn = posizione(x, v[mod1(k + 1, n)])
            g[2v[k]-1] += P * (yn - yp) / 2
            g[2v[k]] += P * (xp - xn) / 2
        end
    end
    return _gradiente_contatto!(g, t, p, x)
end

"Forze sui vertici, −∇E (nN)."
forze(t::Tessuto, p::ParametriMeccanici, x = t.x) = -gradiente!(similar(x), t, p, x)
