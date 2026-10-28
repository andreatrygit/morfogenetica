# Misure per l'osservatore. Nessuna cellula le legge (specifica, sezione 6, passo 9).

export semiangolo, altezza, larghezza_media

# Gli spigoli laterali sinistro e destro di una cellula: dal basale all'apicale.
function _laterali(x, cel::Cellula)
    b1, b2, a2, a1 = cel.vertici
    xb1, yb1 = posizione(x, b1); xa1, ya1 = posizione(x, a1)
    xb2, yb2 = posizione(x, b2); xa2, ya2 = posizione(x, a2)
    return (xa1 - xb1, ya1 - yb1), (xa2 - xb2, ya2 - yb2),
           ((xa1 + xb1) / 2, (ya1 + yb1) / 2), ((xa2 + xb2) / 2, (ya2 + yb2) / 2)
end

"""
    semiangolo(t, c)

Metà dell'angolo fra i due spigoli laterali della cellula `c` (radianti). Positivo se i
laterali convergono verso l'apice, cioè se la cellula è un cuneo stretto in alto.
"""
function semiangolo(t::Tessuto, c::Int, x = t.x)
    u, v, _, _ = _laterali(x, t.cellule[c])
    croce = u[1] * v[2] - u[2] * v[1]
    return atan(croce, u[1] * v[1] + u[2] * v[2]) / 2
end

"Lunghezza media dei due spigoli laterali della cellula `c` (µm)."
function altezza(t::Tessuto, c::Int, x = t.x)
    u, v, _, _ = _laterali(x, t.cellule[c])
    return (hypot(u...) + hypot(v...)) / 2
end

"Distanza fra i punti medi dei due spigoli laterali della cellula `c` (µm)."
function larghezza_media(t::Tessuto, c::Int, x = t.x)
    _, _, m1, m2 = _laterali(x, t.cellule[c])
    return hypot(m2[1] - m1[1], m2[2] - m1[2])
end

export autointersezioni

# I segmenti p1–p2 e q1–q2 si incrociano in senso proprio (non solo toccandosi agli estremi)?
function _si_incrociano(p1, p2, q1, q2)
    o(a, b, c) = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
    d1, d2 = o(q1, q2, p1), o(q1, q2, p2)
    d3, d4 = o(p1, p2, q1), o(p1, p2, q2)
    return d1 * d2 < 0 && d3 * d4 < 0
end

"""
    autointersezioni(t)

Numero di coppie di spigoli senza vertici in comune che si incrociano. Il motore non modella
ancora il contatto fra parti lontane del tessuto: una configurazione con autointersezioni
non è fisica, anche se le forze locali sono in equilibrio.
"""
function autointersezioni(t::Tessuto, x = t.x)
    n = 0
    S = t.spigoli
    for a in 1:length(S), b in a+1:length(S)
        sa, sb = S[a], S[b]
        (sa.i == sb.i || sa.i == sb.j || sa.j == sb.i || sa.j == sb.j) && continue
        _si_incrociano(posizione(x, sa.i), posizione(x, sa.j),
                       posizione(x, sb.i), posizione(x, sb.j)) && (n += 1)
    end
    return n
end
