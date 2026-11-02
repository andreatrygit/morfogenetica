# R4 (docs/roadmap.md): percezione e programma (specifica, sezioni 5, 6 e 8).
#
#   Percezione                  ciò che una regola riceve, e nient'altro: i quattro canali
#     .fase .orologio .eta .stato .intenzioni .misure     dentro di sé
#     .vicine::Vector{Vicina}  (fase e dichiarazioni)     dalle vicine
#     .matrice                                            dalla matrice (vuota per ora)
#     .concentrazioni                                     misura assoluta
#   Uscita(stato, intenzioni, secrezioni, dichiara)       ciò che una regola produce
#   Fase(nome, regola, parametri, fine)                   una riga della tabella di fasi;
#                                                         regola(p, parametri) → Uscita,
#                                                         fine(p, parametri) → Bool
#   Interruttori(; chimica, onda, meccanica, topologia)
#   Embrione(t, programma; parametri, dt, campi, stati, interruttori)
#   passata!(e; ordine)          una passata (specifica, sezione 6); `ordine` è l'ordine in cui
#                                si visitano le cellule, che non deve contare
#   vicine_di(t, c)              le cellule che toccano c (una faccia laterale in comune)
#   e.cellule[c].fase, .stato, .dichiara, .intenzioni, .ingresso   lo stato di ogni cellula
#   e.registro                   l'osservatore: per passata, le fasi di tutte le cellule
#   costringi(; Λ, sostanza, soglia)   la prima primitiva (primitive/costringi.jl): alza la
#                                tensione dell'anello apicale dove la sostanza supera la soglia
#   adatta_alla_scatola!(campo, Lx, Ly)   la griglia periodica segue la scatola
#
# Ogni prova è una funzione che restituisce un Bool. Commit rosso: marcate @test_broken finché
# il codice (motore/13_programma.jl, primitive/costringi.jl) non c'è.

const P_PROGRAMMA = ParametriMeccanici(K = 0.1, η = 1.0, lunghezza_minima = 0.05,
                                       k_contatto = 0.0, d_contatto = 0.5)
const SPENTO = (chimica = false, onda = true, meccanica = false, topologia = false)

foglio_programma(R = 2) = foglio_esagonale(fill(PRISMA, 1 + 3R * (R + 1)), R; lato = 5.0,
                                           altezza = 15.0, bordo_fisso = true)

# Distanze sul grafo dei contatti, calcolate dall'osservatore.
function distanze(t, da)
    d = fill(typemax(Int), length(t.cellule)); d[da] = 0
    coda = [da]
    while !isempty(coda)
        c = popfirst!(coda)
        for v in vicine_di(t, c)
            d[v] == typemax(Int) || continue
            d[v] = d[c] + 1; push!(coda, v)
        end
    end
    return d
end

uscita_ferma(p; stato = p.stato, dichiara = Dict{Symbol,Float64}()) =
    Uscita(copy(stato), p.intenzioni, Dict{Symbol,Float64}(), dichiara)

# Nessun campo della percezione può portare lo stato globale: si percorrono i tipi.
function _tipi_raggiungibili(T, visti = Set{Any}())
    T in visti && return visti
    push!(visti, T)
    if T isa DataType
        for S in fieldtypes(T)
            _tipi_raggiungibili(S, visti)
        end
        for S in T.parameters
            S isa Type && _tipi_raggiungibili(S, visti)
        end
    elseif T isa Union
        _tipi_raggiungibili(T.a, visti); _tipi_raggiungibili(T.b, visti)
    end
    return visti
end

function prova_percezione_senza_stato_globale()
    tipi = _tipi_raggiungibili(Percezione)
    vietati = (Tessuto3D, Campo, Embrione, Faccia, Cellula3D, Any, Function)
    return !any(T -> T in tipi, vietati) && Vicina in tipi
end

function prova_la_regola_riceve_solo_la_percezione()
    visti = Set{Any}()
    regola(p, par) = (push!(visti, (typeof(p), typeof(par))); uscita_ferma(p))
    e = Embrione(foglio_programma(1), [Fase(:f, regola, (a = 1.0,), (p, par) -> false)];
                 parametri = P_PROGRAMMA, dt = 1.0, interruttori = Interruttori(; SPENTO...))
    passata!(e)
    return visti == Set([(Percezione, typeof((a = 1.0,)))])
end

# Due fasi: la cellula con stato[:inizio] = 1 finisce la prima da sé; le altre ci arrivano
# solo per contatto (l'onda).
function programma_onda()
    regola(p, par) = uscita_ferma(p; stato = merge(p.stato, Dict(:orologio => p.orologio,
                                                                :eta => p.eta)))
    return [Fase(:prima, regola, NamedTuple(), (p, par) -> get(p.stato, :inizio, 0.0) == 1),
            Fase(:seconda, regola, NamedTuple(), (p, par) -> false)]
end

function embrione_onda(onda)
    t = foglio_programma(2)
    stati = [Dict{Symbol,Float64}(:inizio => (c == 1 ? 1.0 : 0.0)) for c in eachindex(t.cellule)]
    return Embrione(t, programma_onda(); parametri = P_PROGRAMMA, dt = 1.0, stati,
                    interruttori = Interruttori(; SPENTO..., onda))
end

function prova_onda_monotona()
    e = embrione_onda(true)
    d = distanze(e.tessuto, 1)
    for _ in 1:6
        passata!(e)
    end
    # la cellula a distanza d entra nella seconda fase alla passata d + 2 (la prima passata
    # la cellula 1 si accorge della fine, la successiva è già nella seconda fase)
    ingressi = [e.cellule[c].ingresso for c in eachindex(e.cellule)]
    ingressi == d .+ 2 || return false
    # monotona: nel registro la fase di una cellula non scende mai
    fasi = [r.fasi for r in e.registro]
    all(all(fasi[k+1] .>= fasi[k]) for k in 1:length(fasi)-1) || return false
    # spenta l'onda, avanza solo la cellula che ha visto la fine da sé
    s = embrione_onda(false)
    for _ in 1:6
        passata!(s)
    end
    return [c.fase for c in s.cellule] == [c == 1 ? 2 : 1 for c in eachindex(s.cellule)]
end

function prova_orologio_proprio()
    # burattinaio 3: l'orologio è della cellula, da quando è entrata nella fase
    e = embrione_onda(true)
    P = 6
    for _ in 1:P
        passata!(e)
    end
    return all(e.cellule[c].stato[:orologio] == P - e.cellule[c].ingresso &&
               e.cellule[c].stato[:eta] == P - 1 for c in eachindex(e.cellule))
end

# Una regola che propaga il massimo di un valore per contatto.
function programma_massimo()
    function regola(p, par)
        v = max(get(p.stato, :v, 0.0), maximum((get(n.dichiara, :v, 0.0) for n in p.vicine); init = 0.0))
        return uscita_ferma(p; stato = Dict(:v => v), dichiara = Dict(:v => v))
    end
    return [Fase(:massimo, regola, NamedTuple(), (p, par) -> false)]
end

function embrione_massimo(X, valore)
    t = foglio_programma(2)
    stati = [Dict{Symbol,Float64}(:v => (c == X ? valore : 0.0)) for c in eachindex(t.cellule)]
    return Embrione(t, programma_massimo(); parametri = P_PROGRAMMA, dt = 1.0, stati,
                    interruttori = Interruttori(; SPENTO...))
end

function prova_cono_di_luce()
    # burattinai 1 e 2: nessuna decisione globale. Due embrioni che differiscono solo nello stato
    # della cellula X: dopo k passate le cellule a distanza > k sono identiche, quelle a
    # distanza ≤ k no (la regola è attiva).
    X = 12
    a = embrione_massimo(X, 1.0); b = embrione_massimo(X, 0.0)
    d = distanze(a.tessuto, X)
    for k in 1:3
        passata!(a); passata!(b)
        for c in eachindex(a.cellule)
            uguali = a.cellule[c].stato == b.cellule[c].stato &&
                     a.cellule[c].dichiara == b.cellule[c].dichiara &&
                     a.cellule[c].fase == b.cellule[c].fase
            # dopo k passate l'informazione è arrivata a distanza k - 1 nello stato (la
            # dichiarazione di X si legge alla passata successiva)
            uguali == (d[c] > k - 1) || return false
        end
    end
    return true
end

function prova_ordine_irrilevante()
    a = embrione_massimo(12, 1.0); b = embrione_massimo(12, 1.0)
    ordine = shuffle(MersenneTwister(3), collect(eachindex(a.cellule)))
    for _ in 1:4
        passata!(a); passata!(b; ordine)
    end
    return all(a.cellule[c].stato == b.cellule[c].stato for c in eachindex(a.cellule))
end

function prova_interruttori()
    t = foglio_programma(1)
    campo = Campo(Sostanza(:s, 6.0, 0.1); origine = (-15.0, -15.0, -5.0), h = 2.0, n = (15, 15, 13))
    regola(p, par) = Uscita(copy(p.stato), Intenzioni3D(1.0, 1.0, 1.0, 1.0, 1000.0),
                            Dict(:s => 1.0), Dict{Symbol,Float64}())
    programma = [Fase(:f, regola, NamedTuple(), (p, par) -> false)]
    e = Embrione(t, programma; parametri = P_PROGRAMMA, dt = 1.0, campi = Dict(:s => campo),
                 interruttori = Interruttori(chimica = false, onda = true, meccanica = false,
                                             topologia = false))
    x0 = copy(e.tessuto.x)
    passata!(e); passata!(e)
    (e.tessuto.x == x0 && massa(e.campi[:s]) == 0) || return false
    # accese: la sostanza c'è e il foglio si muove (la tensione di linea stringe gli apici)
    e = Embrione(foglio_programma(1), programma; parametri = P_PROGRAMMA, dt = 1.0,
                 campi = Dict(:s => deepcopy(campo)), interruttori = Interruttori())
    passata!(e); passata!(e)
    return e.tessuto.x != x0 && massa(e.campi[:s]) > 0
end

# costringi, con la sostanza secreta dalla cellula centrale finché il suo orologio lo dice.
function programma_costringi(Λ, soglia, durata_sorgente)
    base = costringi(; Λ, sostanza = :segnale, soglia)
    function regola(p, par)
        u = base.regola(p, base.parametri)
        q = get(p.stato, :sorgente, 0.0) == 1 && p.eta < durata_sorgente ? 20.0 : 0.0
        return Uscita(u.stato, u.intenzioni, merge(u.secrezioni, Dict(:segnale => q)), u.dichiara)
    end
    return [Fase(:costringi, regola, base.parametri, base.fine)]
end

function embrione_costringi(Λ, soglia; durata_sorgente = 10^6, meccanica = true, k = 0.1)
    t = foglio_programma(3)
    campo = Campo(Sostanza(:segnale, 6.0, k); origine = (-40.0, -40.0, -6.0), h = 2.0,
                  n = (40, 40, 14))
    stati = [Dict{Symbol,Float64}(:sorgente => (c == 1 ? 1.0 : 0.0)) for c in eachindex(t.cellule)]
    return Embrione(t, programma_costringi(Λ, soglia, durata_sorgente); parametri = P_PROGRAMMA,
                    dt = 1.0, campi = Dict(:segnale => campo), stati,
                    interruttori = Interruttori(chimica = true, onda = true, meccanica,
                                                topologia = false))
end

rapporto_apice_base(t, c) = area_apicale(t, c) / Morfogenetica._area_basale(t, c)

# La soglia fra il primo e il secondo anello, letta dall'osservatore a regime (prepara la prova,
# non entra nel programma).
function soglia_fra_anelli()
    e = embrione_costringi(0.0, Inf; meccanica = false)
    for _ in 1:40
        passata!(e)
    end
    c = [concentrazione(e.campi[:segnale], centro_cellula(e.tessuto, k)) for k in eachindex(e.cellule)]
    return (minimum(c[2:7]) + maximum(c[8:19])) / 2
end

function prova_costringi_fa_i_cunei()
    soglia = soglia_fra_anelli()
    e = embrione_costringi(3.0, soglia)
    for _ in 1:30
        passata!(e)
    end
    dentro = 1:7; fuori = 8:19
    attive = [e.cellule[c].intenzioni.Λa > 0 for c in eachindex(e.cellule)]
    (all(attive[dentro]) && !any(attive[fuori])) || return false
    r_dentro = maximum(rapporto_apice_base(e.tessuto, c) for c in dentro)
    r_fuori = minimum(rapporto_apice_base(e.tessuto, c) for c in fuori)
    r_dentro < 0.85 && r_fuori > r_dentro || return false
    # interruttore: con Λ = 0 nessun cuneo
    s = embrione_costringi(0.0, soglia)
    for _ in 1:30
        passata!(s)
    end
    return all(abs(rapporto_apice_base(s.tessuto, c) - 1) < 0.02 for c in 1:19)
end

function prova_la_marca_e_un_ipotesi_sul_presente()
    # la sorgente tace dopo 15 passate: la sostanza si degrada e le cellule smettono di
    # costringere da sé, perché la condizione si rivaluta a ogni passata
    soglia = soglia_fra_anelli()
    e = embrione_costringi(3.0, soglia; durata_sorgente = 15, meccanica = false)
    attive_a_metà = 0
    for k in 1:80
        passata!(e)
        k == 15 && (attive_a_metà = count(c -> c.intenzioni.Λa > 0, e.cellule))
    end
    return attive_a_metà == 7 && all(c.intenzioni.Λa == 0 for c in e.cellule)
end

function prova_griglia_segue_la_scatola()
    s = foglio_periodico(fill(PRISMA, 16), 4, 4; lato = 5.0, altezza = 15.0)
    Lx, Ly = s.x[end-1], s.x[end]
    campo = Campo(Sostanza(:s, 6.0, 0.0); origine = (0.0, 0.0, -5.0), h = (Lx / 16, Ly / 16, 2.0),
                  n = (16, 16, 13), bordi = (bordo_periodico, bordo_periodico, bordo_chiuso))
    campo.c .= rand(MersenneTwister(5), size(campo.c)...)
    M = massa(campo)
    c0 = copy(campo.c)
    adatta_alla_scatola!(campo, 1.1Lx, 0.95Ly)
    return isapprox(campo.h[1], 1.1Lx / 16; rtol = 1e-14) && isapprox(campo.h[2], 0.95Ly / 16; rtol = 1e-14) &&
           isapprox(massa(campo), M; rtol = 1e-12) &&
           all(isapprox.(campo.c, c0 ./ (1.1 * 0.95); rtol = 1e-12))
end

function prova_riproducibile()
    a = embrione_costringi(3.0, 0.5; meccanica = true)
    b = embrione_costringi(3.0, 0.5; meccanica = true)
    for _ in 1:3
        passata!(a); passata!(b)
    end
    return a.tessuto.x == b.tessuto.x && a.campi[:segnale].c == b.campi[:segnale].c
end

@testset "R4 percezione e programma" begin
    @test_broken prova_percezione_senza_stato_globale()
    @test_broken prova_la_regola_riceve_solo_la_percezione()
    @test_broken prova_onda_monotona()
    @test_broken prova_orologio_proprio()
    @test_broken prova_cono_di_luce()
    @test_broken prova_ordine_irrilevante()
    @test_broken prova_interruttori()
    @test_broken prova_costringi_fa_i_cunei()
    @test_broken prova_la_marca_e_un_ipotesi_sul_presente()
    @test_broken prova_griglia_segue_la_scatola()
    @test_broken prova_riproducibile()
end
