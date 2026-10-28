# Contatto fra parti lontane del tessuto: esclusione di volume.
#
# Un vertice e uno spigolo che non appartengono a cellule vicine si respingono quando la loro
# distanza d scende sotto d_contatto:
#
#   E_contatto = (k_contatto / 2) · (d_contatto − d)²   per d < d_contatto
#
# È un'energia come le altre, quindi le forze sono appaiate per costruzione. Non è un
# comportamento della cellula ma fisica: due cellule non occupano lo stesso posto.
# Interruttore: k_contatto = 0 spegne il meccanismo.
#
# Per non controllare tutte le coppie a ogni passo si tiene una lista delle coppie vicine
# (lista di Verlet): si ricostruisce quando qualche vertice si è spostato di più di metà del
# margine dall'ultima ricostruzione. A ogni ricostruzione si controlla anche che il tessuto non
# si sia attraversato: se è successo, la configurazione non è ammissibile.

export aggiorna_contatti!, energia_contatto

# Punto dello spigolo p–q più vicino a v: parametro u in [0, 1] e distanza.
function _punto_vicino(v, p, q)
    ex, ey = q[1] - p[1], q[2] - p[2]
    u = clamp(((v[1] - p[1]) * ex + (v[2] - p[2]) * ey) / (ex^2 + ey^2), 0.0, 1.0)
    cx, cy = p[1] + u * ex, p[2] + u * ey
    return u, cx, cy, hypot(v[1] - cx, v[2] - cy)
end

function _distanza_coppia(t::Tessuto, x, (v, k))
    s = t.spigoli[k]
    return _punto_vicino(posizione(x, v), posizione(x, s.i), posizione(x, s.j))
end

"Energia di contatto su tutte le coppie ammissibili (esatta, lenta: per l'osservatore e le prove)."
function energia_contatto(t::Tessuto, p::ParametriMeccanici, x = t.x)
    p.k_contatto == 0 && return 0.0
    E = 0.0
    for coppia in t.contatti.ammissibili
        d = _distanza_coppia(t, x, coppia)[4]
        d < p.d_contatto && (E += p.k_contatto / 2 * (p.d_contatto - d)^2)
    end
    return E
end

"""
    aggiorna_contatti!(t, p, x = t.x)

Ricostruisce la lista delle coppie vicine se qualche vertice si è spostato di più di metà
margine. Il margine è uguale a `d_contatto`.
"""
function aggiorna_contatti!(t::Tessuto, p::ParametriMeccanici, x = t.x)
    p.k_contatto == 0 && return nothing
    C = t.contatti
    margine = p.d_contatto
    if C.valida
        spostamento = 0.0
        for k in 1:length(x)÷2
            spostamento = max(spostamento, hypot(x[2k-1] - C.x_riferimento[2k-1],
                                                 x[2k] - C.x_riferimento[2k]))
        end
        spostamento < margine / 2 && return nothing
    end
    autointersezioni(t, x) > 0 &&
        throw(ConfigurazioneNonAmmissibile("il tessuto si è attraversato"))
    empty!(C.vicine)
    for coppia in C.ammissibili
        _distanza_coppia(t, x, coppia)[4] < p.d_contatto + margine && push!(C.vicine, coppia)
    end
    C.x_riferimento = copy(x)
    C.valida = true
    return nothing
end

# Aggiunge a g il gradiente dell'energia di contatto, sulle sole coppie vicine.
function _gradiente_contatto!(g, t::Tessuto, p::ParametriMeccanici, x)
    p.k_contatto == 0 && return g
    aggiorna_contatti!(t, p, x)
    for (v, k) in t.contatti.vicine
        s = t.spigoli[k]
        pv = posizione(x, v)
        u, cx, cy, d = _punto_vicino(pv, posizione(x, s.i), posizione(x, s.j))
        d >= p.d_contatto && continue
        # dE/dd = −k (d0 − d); ∂d/∂v = n, ∂d/∂p = −(1−u) n, ∂d/∂q = −u n
        f = -p.k_contatto * (p.d_contatto - d)
        nx, ny = (pv[1] - cx) / d, (pv[2] - cy) / d
        g[2v-1] += f * nx; g[2v] += f * ny
        g[2s.i-1] -= f * (1 - u) * nx; g[2s.i] -= f * (1 - u) * ny
        g[2s.j-1] -= f * u * nx; g[2s.j] -= f * u * ny
    end
    return g
end
