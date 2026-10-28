# Geometria del modello a vertici laterale (sezione di un monostrato).
#
# I gradi di libertà sono i vertici. Ogni cellula è un poligono con spigoli etichettati
# apicale, basale o laterale; gli spigoli laterali sono condivisi fra cellule vicine.
# Unità: lunghezze in µm, aree in µm².
#
# Le posizioni stanno in un vettore piatto x = [x₁, y₁, x₂, y₂, …], così energia e forze
# lavorano su un solo array.

export Etichetta, apicale, basale, laterale
export Spigolo, Cellula, Intenzioni, Contatti, Tessuto
export striscia, anello, numero_vertici, posizione, lunghezza, area

@enum Etichetta apicale basale laterale

"""
    Spigolo(i, j, etichetta, cellule)

Segmento fra i vertici `i` e `j`. `cellule` sono le cellule che lo condividono (una o due).
"""
struct Spigolo
    i::Int
    j::Int
    etichetta::Etichetta
    cellule::Vector{Int}
end

"""
    Cellula(vertici)

Poligono con i vertici in ordine antiorario (area positiva).
"""
struct Cellula
    vertici::Vector{Int}
end

"""
    Intenzioni(Γa, Γb, Γl, A0)

Ciò che una cellula chiede alla meccanica: tensioni di linea apicale, basale e laterale (nN) e
area bersaglio (µm²). Sono le "intenzioni meccaniche" della specifica, sezione 6, passo 4.
"""
struct Intenzioni
    Γa::Float64
    Γb::Float64
    Γl::Float64
    A0::Float64
end

"""
    Contatti

Le coppie (vertice, spigolo) che possono toccarsi, e la lista di quelle abbastanza vicine.
Dipende solo dalla topologia: va ricostruita se la topologia cambia.
"""
mutable struct Contatti
    ammissibili::Vector{Tuple{Int,Int}}
    vicine::Vector{Tuple{Int,Int}}
    x_riferimento::Vector{Float64}
    valida::Bool
end

# Un vertice e uno spigolo sono "vicini" (e non si respingono) se una cellula del vertice
# condivide un vertice con una cellula dello spigolo: così le cellule adiacenti restano
# governate solo dalle loro tensioni e dalla loro area.
function Contatti(cellule::Vector{Cellula}, spigoli::Vector{Spigolo}, nv::Int)
    cellule_di = [Int[] for _ in 1:nv]
    for (c, cel) in enumerate(cellule), v in cel.vertici
        push!(cellule_di[v], c)
    end
    # per ogni cellula, l'insieme dei vertici delle cellule che la toccano
    vicinato = [Set{Int}() for _ in cellule]
    for (c, cel) in enumerate(cellule), v in cel.vertici, c2 in cellule_di[v]
        union!(vicinato[c], cellule[c2].vertici)
    end
    ammissibili = Tuple{Int,Int}[]
    for v in 1:nv, (k, s) in enumerate(spigoli)
        (v == s.i || v == s.j) && continue
        vicino = any(c -> v in vicinato[c], s.cellule)
        vicino || push!(ammissibili, (v, k))
    end
    return Contatti(ammissibili, Tuple{Int,Int}[], Float64[], false)
end

"""
    Tessuto

Stato meccanico: posizioni dei vertici, spigoli, cellule, le intenzioni di ogni cellula e la
lista dei contatti possibili fra parti lontane del tessuto.
"""
mutable struct Tessuto
    x::Vector{Float64}
    spigoli::Vector{Spigolo}
    cellule::Vector{Cellula}
    intenzioni::Vector{Intenzioni}
    contatti::Contatti
end

numero_vertici(t::Tessuto) = length(t.x) ÷ 2
posizione(x, k) = (x[2k-1], x[2k])

function lunghezza(x, s::Spigolo)
    xi, yi = posizione(x, s.i)
    xj, yj = posizione(x, s.j)
    return hypot(xj - xi, yj - yi)
end

"Area con segno (formula di Gauss): positiva se i vertici girano in senso antiorario."
function area(x, c::Cellula)
    v = c.vertici
    n = length(v)
    a = zero(eltype(x))
    for k in 1:n
        xk, yk = posizione(x, v[k])
        xn, yn = posizione(x, v[mod1(k + 1, n)])
        a += xk * yn - xn * yk
    end
    return a / 2
end

# Numerazione: il vertice apicale j è 2j-1, quello basale è 2j.
_ap(j) = 2j - 1
_ba(j) = 2j

"""
    striscia(intenzioni; larghezza, altezza)

Striscia aperta e piatta di `length(intenzioni)` cellule rettangolari, basale su y = 0 e
apicale su y = `altezza`, da sinistra a destra.
"""
function striscia(intenzioni::Vector{Intenzioni}; larghezza, altezza)
    n = length(intenzioni)
    x = zeros(2 * 2 * (n + 1))
    for j in 1:n+1
        xa = (j - 1) * larghezza
        x[2_ap(j)-1], x[2_ap(j)] = xa, altezza
        x[2_ba(j)-1], x[2_ba(j)] = xa, 0.0
    end
    return _collega(x, n, n + 1, intenzioni, false)
end

"""
    anello(intenzioni; r_apicale, r_basale)

Anello chiuso di cellule a cuneo, apicale all'interno (raggio `r_apicale`) e basale
all'esterno (raggio `r_basale`).
"""
function anello(intenzioni::Vector{Intenzioni}; r_apicale, r_basale)
    n = length(intenzioni)
    x = zeros(2 * 2n)
    for j in 1:n
        θ = 2π * (j - 1) / n
        x[2_ap(j)-1], x[2_ap(j)] = r_apicale * cos(θ), r_apicale * sin(θ)
        x[2_ba(j)-1], x[2_ba(j)] = r_basale * cos(θ), r_basale * sin(θ)
    end
    return _collega(x, n, n, intenzioni, true)
end

# Costruisce cellule e spigoli per una fila di n cellule su m coppie di vertici.
function _collega(x, n, m, intenzioni, chiuso)
    succ(j) = chiuso ? mod1(j + 1, m) : j + 1
    cellule = Cellula[]
    spigoli = Spigolo[]
    for i in 1:n
        j, k = i, succ(i)
        # antiorario: basale, basale successivo, apicale successivo, apicale.
        # Vale per la striscia (basale in basso) e per l'anello (basale all'esterno).
        push!(cellule, Cellula([_ba(j), _ba(k), _ap(k), _ap(j)]))
        push!(spigoli, Spigolo(_ap(j), _ap(k), apicale, [i]))
        push!(spigoli, Spigolo(_ba(j), _ba(k), basale, [i]))
    end
    for j in 1:m
        vicine = Int[]
        chiuso ? push!(vicine, mod1(j - 1, n)) : (j > 1 && push!(vicine, j - 1))
        j <= n && push!(vicine, j)
        push!(spigoli, Spigolo(_ap(j), _ba(j), laterale, unique(vicine)))
    end
    nv = length(x) ÷ 2
    return Tessuto(x, spigoli, cellule, copy(intenzioni), Contatti(cellule, spigoli, nv))
end
