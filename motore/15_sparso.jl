# Il solutore sparso (docs/roadmap.md, R6): `minimizza!` di 06_integratore.jl per tessuti di
# migliaia di cellule.
#
# `minimizza!` costruisce la matrice delle derivate seconde piena, per differenze centrali del
# gradiente esatto: 2n valutazioni del gradiente e n² numeri, con n il numero di incognite. Con
# 1700 cellule n è circa 20 000, e non regge. Ma la matrice è quasi tutta zeri: l'energia è una
# somma di termini di cellula (volume, facce, anello apicale), quindi due vertici sono accoppiati
# solo se stanno in una stessa cellula. Si sfrutta in due modi:
#
#   1. colorazione a distanza due: due vertici che non hanno nessun vicino in comune si possono
#      perturbare insieme, perché le righe che cambiano non si sovrappongono. Con una colorazione
#      greedy bastano qualche decina di colori, e la matrice costa 2 · 3 · colori valutazioni del
#      gradiente invece di 2n. E ognuna è parziale: il gradiente è una somma su facce e cellule,
#      e perturbando un gruppo cambiano solo i termini delle cellule che contengono i suoi
#      vertici; lo si calcola su un sotto-tessuto con quelle cellule e le loro facce (circa un
#      quarto del foglio). Le differenze sono le stesse, a meno dell'arrotondamento;
#   2. fattorizzazione di Cholesky sparsa (CHOLMOD) al posto di quella piena.
#
# I gruppi si calcolano in parallelo sui fili di Julia (`julia --threads=auto`; la batteria li
# usa), con lo stesso risultato bit per bit qualunque sia il numero di fili.
#
# Nei fogli periodici le due voci della scatola (Lx, Ly) sono accoppiate a tutto: si perturbano da
# sole, e le loro righe si ricavano per simmetria dalle loro colonne.
#
# Il metodo è lo stesso di `minimizza!`: Newton smorzato (Levenberg–Marquardt), passo accettato
# solo se l'energia scende, spostamento limitato per vertice. Solo per `Tessuto3D`.

using SparseArrays: sparse, SparseMatrixCSC
using LinearAlgebra: cholesky!

export minimizza_sparso!

# Chi è accoppiato con chi: per ogni vertice, i vertici delle cellule che lo contengono (lui
# compreso). Più la colorazione a distanza due dei vertici.
function _struttura(t::Tessuto3D)
    nv = _nvert(t)
    vicini = [Set{Int}() for _ in 1:nv]
    for cel in t.cellule
        vv = unique(v for fi in cel.facce for v in t.facce[fi].vertici)
        for v in vv
            union!(vicini[v], vv)
        end
    end
    vicini = [sort!(collect(s)) for s in vicini]
    colore = zeros(Int, nv)
    vietati = Set{Int}()
    for v in 1:nv
        empty!(vietati)
        for u in vicini[v], w in vicini[u]
            colore[w] > 0 && push!(vietati, colore[w])
        end
        k = 1
        while k in vietati
            k += 1
        end
        colore[v] = k
    end
    gruppi = [findall(==(k), colore) for k in 1:maximum(colore; init = 0)]
    return (vicini = vicini, gruppi = gruppi, sotto = [_sotto_tessuto(t, g) for g in gruppi])
end

# Il sotto-tessuto di un gruppo di vertici: le cellule che ne contengono almeno uno, con tutte le
# loro facce (rinumerate); le altre cellule restano senza facce, così la numerazione delle cellule
# e le intenzioni non cambiano. Posizioni, scatola e vertici bloccati sono quelli del tessuto.
function _sotto_tessuto(t::Tessuto3D, gruppo)
    nel = falses(_nvert(t)); nel[gruppo] .= true
    toccate = [any(nel[v] for fi in cel.facce for v in t.facce[fi].vertici) for cel in t.cellule]
    facce = sort!(unique(fi for (c, cel) in enumerate(t.cellule) if toccate[c] for fi in cel.facce))
    nuovo = Dict(fi => k for (k, fi) in enumerate(facce))
    cellule = [toccate[c] ? Cellula3D([nuovo[fi] for fi in cel.facce], cel.segni, cel.traslazioni) :
                            Cellula3D(Int[], Int[], NTuple{2,Int}[]) for (c, cel) in enumerate(t.cellule)]
    return Tessuto3D(t.x, t.facce[facce], cellule, t.intenzioni, t.periodico, t.fissi, t.ancore)
end

"""
    _hessiana_sparsa(t, p, x; struttura = _struttura(t), ancora = nothing, peso = 0.0)

La matrice delle derivate seconde di E(x) + (peso/2)·|x − ancora|² in `x`, sparsa, per
differenze centrali del gradiente esatto (passo 1e-6, come `minimizza!`), simmetrizzata.
"""
function _hessiana_sparsa(t::Tessuto3D, p::ParametriMeccanici, x;
                          struttura = _struttura(t), ancora = nothing, peso = 0.0)
    function grad!(g, T, y)
        gradiente!(g, T, p, y)
        ancora === nothing || (g .+= peso .* (y .- ancora))
        return g
    end
    nv = _nvert(t)
    n = length(x)
    ε = 1e-6
    # i gruppi (per ogni coordinata) sono indipendenti: si calcolano in parallelo, ognuno con le
    # sue copie di lavoro, e si raccolgono nel loro ordine, così la matrice è la stessa bit per
    # bit qualunque sia il numero di fili
    compiti = [(k, a) for k in eachindex(struttura.gruppi) for a in 1:3]
    pezzi = Vector{Tuple{Vector{Int},Vector{Int},Vector{Float64}}}(undef, length(compiti))
    Threads.@threads for q in eachindex(compiti)
        k, a = compiti[q]
        gruppo, T = struttura.gruppi[k], struttura.sotto[k]
        # `local`: senza, Julia le confonderebbe con le variabili omonime della funzione, uniche
        # e condivise da tutti i fili
        local y = copy(x)
        local gp = similar(x)
        local gm = similar(x)
        for v in gruppo
            y[3v-3+a] = x[3v-3+a] + ε
        end
        grad!(gp, T, y)
        for v in gruppo
            y[3v-3+a] = x[3v-3+a] - ε
        end
        grad!(gm, T, y)
        local Iq = Int[]
        local Jq = Int[]
        local Vq = Float64[]
        for v in gruppo
            j = 3v - 3 + a
            for w in struttura.vicini[v], b in 1:3
                i = 3w - 3 + b
                push!(Iq, i); push!(Jq, j); push!(Vq, (gp[i] - gm[i]) / 2ε)
            end
        end
        pezzi[q] = (Iq, Jq, Vq)
    end
    # una sola concatenazione per vettore, nell'ordine dei compiti (reduce con un generatore
    # ricopierebbe tutto a ogni aggiunta)
    I = reduce(vcat, [pz[1] for pz in pezzi])
    J = reduce(vcat, [pz[2] for pz in pezzi])
    V = reduce(vcat, [pz[3] for pz in pezzi])
    y = copy(x); gp = similar(x); gm = similar(x)
    # le voci della scatola: colonne piene, e le righe per simmetria
    for j in 3nv+1:n
        y[j] = x[j] + ε; grad!(gp, t, y)
        y[j] = x[j] - ε; grad!(gm, t, y)
        y[j] = x[j]
        for i in 1:n
            d = (gp[i] - gm[i]) / 2ε
            push!(I, i); push!(J, j); push!(V, d)
            i <= 3nv && (push!(I, j); push!(J, i); push!(V, d))
        end
    end
    H = sparse(I, J, V, n, n)
    return (H + H') / 2
end

"""
    minimizza_sparso!(t, p; tolleranza, iterazioni_massime, spostamento_massimo = 1.0,
                      ancora = nothing, peso = 0.0)

Come `minimizza!` (06_integratore.jl), con la matrice delle derivate seconde sparsa: porta il
tessuto 3D al minimo di E(x) + (peso/2)·|x − ancora|², finché la forza massima su un vertice
scende sotto `tolleranza` (nN). Restituisce `(convergito, iterazioni, forza_massima)`.
"""
function minimizza_sparso!(t::Tessuto3D, p::ParametriMeccanici; tolleranza, iterazioni_massime,
                           spostamento_massimo = 1.0, ancora = nothing, peso = 0.0)
    p.k_contatto > 0 && (spostamento_massimo = min(spostamento_massimo, p.d_contatto / 4))
    x = t.x
    obiettivo(y) = energia(t, p, y) + (ancora === nothing ? 0.0 : peso / 2 * sum(abs2, y .- ancora))
    function grad!(g, y)
        gradiente!(g, t, p, y)
        ancora === nothing || (g .+= peso .* (y .- ancora))
        return g
    end
    struttura = _struttura(t)
    g = similar(x); gp = similar(x)
    simbolica = nothing
    μ = 1e-3
    f = Inf
    for it in 1:iterazioni_massime
        grad!(g, x)
        f = forza_per_vertice(t, g)
        if f < tolleranza
            _controlla_contatto(t, p)
            return (convergito = true, iterazioni = it, forza_massima = f)
        end
        H = _hessiana_sparsa(t, p, x; struttura, ancora, peso)
        grad!(g, x)
        E0 = obiettivo(x)
        accettato = false
        for _ in 1:60
            # l'analisi simbolica (l'ordinamento) dipende solo dalla struttura della matrice, che
            # non cambia: la prima fattorizzazione la fa, le altre la riusano
            A = Symmetric(H + μ * I)
            if simbolica === nothing
                C = cholesky(A; check = false)
                simbolica = C
            else
                C = cholesky!(simbolica, A; check = false)
            end
            if !issuccess(C)
                # una fattorizzazione fallita lascia l'oggetto in uno stato di errore: la prossima
                # ricomincia da capo (con il riuso cieco, μ cresceva senza fine)
                simbolica = nothing
                μ *= 4
                continue
            end
            d = -(C \ g)
            s = forza_per_vertice(t, d)
            s > spostamento_massimo && (d .*= spostamento_massimo / s)
            prova = x .+ d
            ammissibile = try
                controlla_ammissibile(t, p, prova)
                true
            catch e
                e isa ConfigurazioneNonAmmissibile || rethrow()
                false
            end
            Ep = ammissibile ? obiettivo(prova) : Inf
            scende = Ep < E0 || (abs(Ep - E0) <= 1e-12 * abs(E0) &&
                                 forza_per_vertice(t, grad!(gp, prova)) < f)
            if ammissibile && scende
                x .= prova
                μ = max(μ / 3, 1e-12)
                accettato = true
                break
            end
            μ *= 4
        end
        accettato || break
    end
    grad!(g, x)
    f = forza_per_vertice(t, g)
    _controlla_contatto(t, p)
    return (convergito = f < tolleranza, iterazioni = iterazioni_massime, forza_massima = f)
end
