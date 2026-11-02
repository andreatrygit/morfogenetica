# Percezione e programma (specifica, sezioni 5, 6 e 8; docs/roadmap.md, passo R4).
#
# Un programma è una tabella di fasi copiata dentro ogni cellula, come un genoma. Ogni fase ha
# una regola, i suoi parametri e un predicato di fine. La regola è una funzione pura della sola
# percezione: riceve una `Percezione`, che porta i quattro canali e nient'altro, e restituisce
# un'`Uscita` (stato interno nuovo, intenzioni meccaniche, secrezioni, ciò che dichiara alle
# vicine). La percezione non ha modo di raggiungere il tessuto, i campi o le altre cellule: il
# tipo lo impedisce (nessun campo di tipo Tessuto3D, Campo, Any o Function).
#
# I quattro canali (specifica, sezione 8):
#   1. dentro di sé: fase, orologio (passate nella fase), età (passate di vita), stato interno,
#      intenzioni a riposo, misure del proprio corpo (volume, area apicale, area basale);
#   2. dalle vicine: per ogni cellula a contatto la sua fase e ciò che ha dichiarato alla
#      passata precedente;
#   3. dalla matrice: vuoto, finché la matrice non c'è;
#   4. misura assoluta: la concentrazione di ogni sostanza alla propria posizione.
#
# Avanzamento di fase, monotono e senza orologio centrale: una cellula passa da k a k+1 quando
# il predicato di fine di k è vero per lei, oppure (l'onda) quando una vicina è già oltre k.
# Mai indietro, mai più di una fase per passata.
#
# Una passata (specifica, sezione 6), qui nei passi che ci sono già:
#   1. vicinato: chi tocca chi (una faccia laterale in comune);
#   2. chimica: ogni campo avanza di dt, con le secrezioni decise alla passata precedente; nei
#      fogli periodici la griglia segue prima la scatola;
#   3. percezione di ogni cellula, dallo stato della passata precedente;
#   4. programma: regola e avanzamento di fase, per tutte le cellule sullo stesso stato
#      (sincrono: l'ordine delle cellule non conta);
#   5. meccanica: le intenzioni diventano le tensioni del tessuto, un passo implicito di dt;
#   6. topologia: le operazioni elementari dovute (11_operazioni.jl);
#   9. registro: l'osservatore annota le fasi. Può guardare tutto; nessuna cellula lo legge.
# Mancano la popolazione (divisione e morte, passo 7) e la matrice (passo 8).
#
# Ogni meccanismo ha un interruttore (`Interruttori`).

export Vicina, Percezione, Uscita, Fase, Interruttori, StatoCellula, Embrione
export passata!, vicine_di

"""
    Vicina(fase, dichiara)

Ciò che una cellula sa di una vicina: la sua fase e ciò che ha dichiarato.
"""
struct Vicina
    fase::Int
    dichiara::Dict{Symbol,Float64}
end

const Misure = NamedTuple{(:volume, :area_apicale, :area_basale),NTuple{3,Float64}}

"""
    Percezione

Ciò che la regola di una fase riceve: i quattro canali, e nient'altro (vedi l'intestazione di
`motore/13_programma.jl`). `intenzioni` sono quelle a riposo della cellula; `orologio` conta le
passate nella fase corrente (0 alla prima), `eta` quelle vissute prima di questa.
"""
struct Percezione
    fase::Int
    orologio::Int
    eta::Int
    stato::Dict{Symbol,Float64}
    intenzioni::Intenzioni3D
    misure::Misure
    vicine::Vector{Vicina}
    matrice::Dict{Symbol,Float64}
    concentrazioni::Dict{Symbol,Float64}
end

"""
    Uscita(stato, intenzioni, secrezioni, dichiara)

Ciò che la regola produce: il nuovo stato interno, le intenzioni meccaniche (tensioni e volume
bersaglio), le secrezioni (quantità/min per sostanza) e ciò che la cellula dichiara alle vicine.
"""
struct Uscita
    stato::Dict{Symbol,Float64}
    intenzioni::Intenzioni3D
    secrezioni::Dict{Symbol,Float64}
    dichiara::Dict{Symbol,Float64}
end

"""
    Fase(nome, regola, parametri, fine)

Una riga della tabella di fasi: `regola(p, parametri)` → `Uscita` e `fine(p, parametri)` →
Bool, funzioni pure della percezione `p`.
"""
struct Fase{R,F,P<:NamedTuple}
    nome::Symbol
    regola::R
    parametri::P
    fine::F
end

"""
    Interruttori(; chimica = true, onda = true, meccanica = true, topologia = true)

Un interruttore per meccanismo: spento, il passo corrispondente della passata non si fa
(`onda` spenta: si avanza di fase solo per la propria fine).
"""
Base.@kwdef struct Interruttori
    chimica::Bool = true
    onda::Bool = true
    meccanica::Bool = true
    topologia::Bool = true
end

"""
    StatoCellula

Ciò che il motore tiene per ogni cellula: la sua copia del programma, la fase e la passata in
cui vi è entrata (`ingresso`), la passata di nascita, lo stato interno, ciò che dichiara, le
intenzioni a riposo e quelle correnti, le secrezioni correnti.
"""
mutable struct StatoCellula
    programma::Tuple
    fase::Int
    ingresso::Int
    nascita::Int
    stato::Dict{Symbol,Float64}
    dichiara::Dict{Symbol,Float64}
    riposo::Intenzioni3D
    intenzioni::Intenzioni3D
    secrezioni::Dict{Symbol,Float64}
end

"""
    Embrione(t, programma; parametri, dt, campi = Dict(), stati = nothing,
             interruttori = Interruttori(), tolleranza = 1e-8,
             operazioni = (soglia = 0.3, area_soglia = 0.3, lunghezza_nuova = 1.0))

Un tessuto con un programma copiato in ogni cellula, i campi chimici per nome di sostanza, i
parametri meccanici e il passo meccanico `dt` (min). `stati[c]` è lo stato interno iniziale
della cellula c (vuoto se non dato). Le intenzioni a riposo sono quelle del tessuto.
"""
mutable struct Embrione
    tessuto::Tessuto3D
    cellule::Vector{StatoCellula}
    campi::Dict{Symbol,Campo}
    parametri::ParametriMeccanici
    dt::Float64
    passata::Int
    interruttori::Interruttori
    tolleranza::Float64
    operazioni::NamedTuple
    registro::Vector{NamedTuple{(:passata, :fasi),Tuple{Int,Vector{Int}}}}
end

function Embrione(t::Tessuto3D, programma; parametri, dt, campi = Dict{Symbol,Campo}(),
                  stati = nothing, interruttori = Interruttori(), tolleranza = 1e-8,
                  operazioni = (soglia = 0.3, area_soglia = 0.3, lunghezza_nuova = 1.0))
    tabella = Tuple(programma)
    isempty(tabella) && throw(ArgumentError("il programma non ha fasi"))
    cellule = [StatoCellula(tabella, 1, 1, 1,
                            stati === nothing ? Dict{Symbol,Float64}() : Dict{Symbol,Float64}(stati[c]),
                            Dict{Symbol,Float64}(), t.intenzioni[c], t.intenzioni[c],
                            Dict{Symbol,Float64}()) for c in eachindex(t.cellule)]
    return Embrione(t, cellule, Dict{Symbol,Campo}(campi), parametri, Float64(dt), 0,
                    interruttori, tolleranza, operazioni, NamedTuple{(:passata, :fasi),Tuple{Int,Vector{Int}}}[])
end

"""
    vicine_di(t, c)

Le cellule che toccano la cellula `c`: quelle con cui ha una faccia in comune.
"""
function vicine_di(t::Tessuto3D, c::Int)
    v = Int[]
    for fi in t.cellule[c].facce
        f = t.facce[fi]
        length(f.cellule) == 2 || continue
        push!(v, f.cellule[1] == c ? f.cellule[2] : f.cellule[1])
    end
    return sort!(unique!(v))
end

# Passo 3: la percezione della cellula c, dallo stato della passata precedente. Tutto è copiato:
# la regola non può toccare lo stato di nessuno.
function _percezione(e::Embrione, c::Int, vicine)
    t = e.tessuto
    s = e.cellule[c]
    misure = (volume = volume_cellula(t, c), area_apicale = area_apicale(t, c),
              area_basale = _area_basale(t, c))
    posizione = isempty(e.campi) ? (0.0, 0.0, 0.0) : centro_cellula(t, c)
    concentrazioni = Dict{Symbol,Float64}(nome => concentrazione(campo, posizione)
                                          for (nome, campo) in e.campi)
    return Percezione(s.fase, e.passata - s.ingresso, e.passata - s.nascita, copy(s.stato),
                      s.riposo, misure,
                      [Vicina(e.cellule[v].fase, copy(e.cellule[v].dichiara)) for v in vicine],
                      Dict{Symbol,Float64}(), concentrazioni)
end

"""
    passata!(e; ordine = eachindex(e.cellule))

Una passata del motore (vedi l'intestazione del file). `ordine` è l'ordine in cui si visitano
le cellule al passo 4: non conta, perché tutte leggono lo stato della passata precedente.
Restituisce l'embrione.
"""
function passata!(e::Embrione; ordine = eachindex(e.cellule))
    sort(collect(ordine)) == collect(eachindex(e.cellule)) ||
        throw(ArgumentError("l'ordine deve visitare ogni cellula una volta"))
    e.passata += 1
    t = e.tessuto
    # 1. vicinato
    vicine = [vicine_di(t, c) for c in eachindex(e.cellule)]
    # 2. chimica
    if e.interruttori.chimica
        for (nome, campo) in e.campi
            t.periodico && campo.bordi[1] == bordo_periodico && adatta_alla_scatola!(campo, _periodo(t, t.x)...)
            scambi = [Scambio(centro_cellula(t, c), get(s.secrezioni, nome, 0.0), 0.0)
                      for (c, s) in enumerate(e.cellule) if get(s.secrezioni, nome, 0.0) != 0]
            avanza!(campo, e.dt, scambi)
        end
    end
    # 3–4. percezione e programma, tutte sullo stato precedente
    uscite = Vector{Uscita}(undef, length(e.cellule))
    avanza = falses(length(e.cellule))
    for c in ordine
        s = e.cellule[c]
        p = _percezione(e, c, vicine[c])
        fase = s.programma[s.fase]
        uscite[c] = fase.regola(p, fase.parametri)
        if s.fase < length(s.programma)
            avanza[c] = fase.fine(p, fase.parametri) ||
                        (e.interruttori.onda && any(v.fase > s.fase for v in p.vicine))
        end
    end
    for (c, s) in enumerate(e.cellule)
        u = uscite[c]
        s.stato = copy(u.stato); s.dichiara = copy(u.dichiara)
        s.intenzioni = u.intenzioni; s.secrezioni = copy(u.secrezioni)
        if avanza[c]
            s.fase += 1
            s.ingresso = e.passata + 1
        end
    end
    # 5. meccanica
    if e.interruttori.meccanica
        for (c, s) in enumerate(e.cellule)
            t.intenzioni[c] = s.intenzioni
        end
        passo_implicito!(t, e.parametri, e.dt; tolleranza = e.tolleranza)
    end
    # 6. topologia
    if e.interruttori.topologia
        _operazioni_dovute!(t; e.operazioni..., monostrato = true)
    end
    # 9. registro dell'osservatore
    push!(e.registro, (passata = e.passata, fasi = [s.fase for s in e.cellule]))
    return e
end
