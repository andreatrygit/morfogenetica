# Integratore: minimizzazione dell'energia e passo implicito della dinamica sovrasmorzata.
#
# Eulero esplicito (03_rilassamento.jl) è limitato dalla parte più rigida del tessuto
# (penalità d'area, contatto): il passo deve restare piccolo anche quando il tessuto si muove
# piano, come quando una striscia lunga si piega. Qui due strumenti che non hanno quel limite.
#
# 1. `minimizza!`: cerca direttamente il minimo dell'energia con il metodo di Newton smorzato
#    (Levenberg–Marquardt): a ogni iterazione risolve (H + μI)·d = −∇E, con H la matrice delle
#    derivate seconde, e accetta il passo solo se l'energia scende; altrimenti alza μ e
#    riprova. Vicino al minimo converge in poche iterazioni anche dove Eulero e FIRE ne
#    richiedono milioni (una striscia lunga che si piega). H si calcola per differenze del
#    gradiente esatto: va bene per i tessuti 2D di qualche centinaio di vertici, non per il 3D.
#    Per le prove di equilibrio, dove conta solo la forma finale.
#
# 2. `passo_implicito!`: Eulero implicito per η·dx/dt = −∇E. Per la dinamica sovrasmorzata il
#    passo implicito di durata dt è esattamente il minimo di
#
#        E(x) + (η / 2dt) · |x − x_n|²
#
#    (lo "spostamento minimizzante"): stabile per qualunque dt, l'energia non sale mai da un
#    passo all'altro, e per dt → 0 ritrova la traiettoria di Eulero esplicito. Quando conta il
#    tempo fisico.
#
# Le iterazioni di Newton non hanno un tempo; il tempo fisico è solo quello di
# `passo_implicito!` ed `evolvi!`, in minuti.

using LinearAlgebra: Symmetric, cholesky, issuccess, I

export minimizza!, passo_implicito!, evolvi!

"""
    minimizza!(t, p; tolleranza, iterazioni_massime, spostamento_massimo = 1.0,
               ancora = nothing, peso = 0.0)

Porta il tessuto al minimo di E(x) + (peso/2)·|x − ancora|², finché la forza massima su un
vertice (nN) scende sotto `tolleranza`. Nessun vertice si sposta più di `spostamento_massimo`
(µm) per iterazione; con il contatto acceso il limite scende a d_contatto/4, così il contatto
non può essere scavalcato.

Restituisce `(convergito, iterazioni, forza_massima)`.
"""
function minimizza!(t::Tessuto, p::ParametriMeccanici; tolleranza, iterazioni_massime,
                    spostamento_massimo = 1.0, ancora = nothing, peso = 0.0)
    p.k_contatto > 0 && (spostamento_massimo = min(spostamento_massimo, p.d_contatto / 4))
    x = t.x
    n = length(x)
    obiettivo(y) = energia(t, p, y) + (ancora === nothing ? 0.0 : peso / 2 * sum(abs2, y .- ancora))
    function grad!(g, y)
        gradiente!(g, t, p, y)
        ancora === nothing || (g .+= peso .* (y .- ancora))
        return g
    end
    g = similar(x); gp = similar(x); gm = similar(x)
    H = zeros(n, n)
    μ = 1e-3
    f = Inf
    for it in 1:iterazioni_massime
        grad!(g, x)
        f = _forza_massima(g)
        if f < tolleranza
            _controlla_contatto(t, p)
            return (convergito = true, iterazioni = it, forza_massima = f)
        end
        # matrice delle derivate seconde per differenze centrali del gradiente esatto
        ε = 1e-6
        for j in 1:n
            xj = x[j]
            x[j] = xj + ε; grad!(gp, x)
            x[j] = xj - ε; grad!(gm, x)
            x[j] = xj
            H[:, j] .= (gp .- gm) ./ 2ε
        end
        grad!(g, x)                    # ripristina la lista dei contatti su x
        E0 = obiettivo(x)
        accettato = false
        for _ in 1:60
            C = cholesky(Symmetric((H .+ H') ./ 2 + μ * I); check = false)
            if !issuccess(C)
                μ *= 4
                continue
            end
            d = -(C \ g)
            s = _forza_massima(d)
            s > spostamento_massimo && (d .*= spostamento_massimo / s)
            prova = x .+ d
            ammissibile = try
                controlla_ammissibile(t, p, prova)
                p.k_contatto == 0 || autointersezioni(t, prova) == 0
            catch e
                e isa ConfigurazioneNonAmmissibile || rethrow()
                false
            end
            # vicino al minimo l'energia cambia meno della precisione dei numeri: lì si
            # accetta il passo se riduce la forza
            Ep = ammissibile ? obiettivo(prova) : Inf
            scende = Ep < E0 || (abs(Ep - E0) <= 1e-12 * abs(E0) &&
                                 _forza_massima(grad!(gp, prova)) < f)
            if ammissibile && scende
                x .= prova
                μ = max(μ / 3, 1e-12)
                accettato = true
                break
            end
            μ *= 4
        end
        accettato || break             # nessun passo fa scendere l'energia: si è al minimo numerico
    end
    grad!(g, x)
    f = _forza_massima(g)
    _controlla_contatto(t, p)
    return (convergito = f < tolleranza, iterazioni = iterazioni_massime, forza_massima = f)
end

"""
    passo_implicito!(t, p, dt; tolleranza = 1e-9, iterazioni_massime = 1_000_000)

Un passo di Eulero implicito di durata `dt` (min): il minimo di E(x) + (η/2dt)·|x − x_n|².
`tolleranza` è il residuo ammesso sull'equazione del passo (nN).
"""
function passo_implicito!(t::Tessuto, p::ParametriMeccanici, dt; tolleranza = 1e-9,
                          iterazioni_massime = 200)
    r = minimizza!(t, p; tolleranza, iterazioni_massime, ancora = copy(t.x), peso = p.η / dt)
    r.convergito || error("passo implicito non risolto: residuo $(r.forza_massima) nN")
    return r
end

"""
    evolvi!(t, p; durata, dt, tolleranza = 1e-9)

Fa evolvere il tessuto per `durata` minuti con passi impliciti di `dt` minuti.
Restituisce il numero di passi.
"""
function evolvi!(t::Tessuto, p::ParametriMeccanici; durata, dt, tolleranza = 1e-9)
    n = max(1, round(Int, durata / dt))
    for _ in 1:n
        passo_implicito!(t, p, durata / n; tolleranza)
    end
    return n
end
