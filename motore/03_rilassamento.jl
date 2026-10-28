# Dinamica sovrasmorzata: η · dx/dt = −∇E, integrata con Eulero esplicito.
# Il tempo è in minuti.

export ConfigurazioneNonAmmissibile, controlla_ammissibile, passo!, rilassa!

"""
    ConfigurazioneNonAmmissibile(motivo)

Il motore si rifiuta di proseguire: uno spigolo è sceso sotto la lunghezza minima o una cellula
si è rovesciata. Non si simula in silenzio una geometria che il modello non descrive.
"""
struct ConfigurazioneNonAmmissibile <: Exception
    motivo::String
end

function controlla_ammissibile(t::Tessuto, p::ParametriMeccanici, x = t.x)
    for s in t.spigoli
        L = lunghezza(x, s)
        L < p.lunghezza_minima &&
            throw(ConfigurazioneNonAmmissibile("spigolo $(s.etichetta) $(s.i)–$(s.j) lungo $L µm"))
    end
    for (c, cel) in enumerate(t.cellule)
        a = area(x, cel)
        a <= 0 && throw(ConfigurazioneNonAmmissibile("cellula $c con area $a µm²"))
    end
    return nothing
end

"""
    passo!(t, p, dt, g = similar(t.x))

Un passo meccanico di durata `dt` (min). Restituisce la forza massima su un vertice (nN)
prima del passo.
"""
function passo!(t::Tessuto, p::ParametriMeccanici, dt, g = similar(t.x))
    gradiente!(g, t, p)
    t.x .-= (dt / p.η) .* g
    controlla_ammissibile(t, p)
    return _forza_massima(g)
end

function _forza_massima(g)
    m = 0.0
    for k in 1:length(g)÷2
        m = max(m, hypot(g[2k-1], g[2k]))
    end
    return m
end

"""
    rilassa!(t, p; dt, tolleranza, passi_massimi)

Ripete `passo!` finché la forza massima su un vertice scende sotto `tolleranza` (nN).
Restituisce `(convergito, passi, forza_massima)`.
"""
function rilassa!(t::Tessuto, p::ParametriMeccanici; dt, tolleranza, passi_massimi)
    g = similar(t.x)
    f = Inf
    for n in 1:passi_massimi
        f = passo!(t, p, dt, g)
        f < tolleranza && return (convergito = true, passi = n, forza_massima = f)
    end
    return (convergito = false, passi = passi_massimi, forza_massima = f)
end
