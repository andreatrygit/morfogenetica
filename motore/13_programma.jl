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
#   1. dentro di sé: fase, orologio (passate nella fase), età (passate di vita), il proprio
#      numero a caso, stato interno, intenzioni a riposo, misure del proprio corpo (volume, area
#      apicale, area basale);
#   2. dalle vicine: per ogni cellula a contatto la sua fase e ciò che ha dichiarato alla
#      passata precedente;
#   3. dalla matrice: la densità di ogni matrice al proprio lato basale;
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
#   5. meccanica: le intenzioni diventano le tensioni del tessuto, un passo implicito di dt
#      (se non converge si registra e si va avanti: vedi sotto);
#   6. topologia: le operazioni elementari dovute, nel modello generale (11_operazioni.jl),
#      oppure nel guscio a prismi uno scambio di vicine (08), il più corto, e le estrusioni
#      (09). Uno solo per passata, con la meccanica in mezzo: una cellula in apoptosi arriva
#      alla soglia con tutti i lati insieme, e più scambi nella stessa volta, ognuno con una
#      parete nuova più lunga della cellula, la rovesciano (volume negativo);
#   7. popolazione, ogni `ogni` passate: si accolgono le richieste del passo 4. Chi chiede di
#      morire entra in apoptosi (programma fermo, volume bersaglio che scende con tempo
#      caratteristico `tempo_apoptosi` fino a `volume_apoptosi`) e uscirà dal foglio con
#      un'estrusione; chi chiede di dividersi si divide (14_popolazione.jl), e la figlia
#      eredita fase, stato e programma, con un'identità sua. Solo nel guscio a prismi;
#   8. matrice: ogni cellula depone sul proprio lato basale; il deposto non diffonde e resta
#      anche dopo di lei;
#   9. registro: l'osservatore annota le fasi e gli eventi. Può guardare tutto; nessuna cellula
#      lo legge.
#
# Casualità: ogni cellula riceve a ogni passata un numero in [0, 1) che dipende solo dal seme
# dell'embrione, dalla sua identità e dalla passata. Le corse sono riproducibili e l'ordine delle
# cellule non conta.
#
# Ogni meccanismo ha un interruttore (`Interruttori`).

export Vicina, Percezione, Uscita, Fase, SCRIVIBILI, Interruttori, StatoCellula, Embrione
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
passate nella fase corrente (0 alla prima), `eta` quelle vissute prima di questa; `caso` è un
numero in [0, 1) dal seme, proprio della cellula e della passata.
"""
struct Percezione
    fase::Int
    orologio::Int
    eta::Int
    caso::Float64
    stato::Dict{Symbol,Float64}
    intenzioni::Intenzioni3D
    misure::Misure
    vicine::Vector{Vicina}
    matrice::Dict{Symbol,Float64}
    concentrazioni::Dict{Symbol,Float64}
end

"""
    Uscita(; intenzioni, stato, dichiara, secrezioni, deposita, richiesta = nothing, piano = NaN)

Ciò che la regola produce, come delta: contiene solo ciò che la regola scrive, e la presenza
di una voce è il registro delle scritture (una scrittura del valore a riposo si distingue da
un "non tocco"). `intenzioni`: i campi di `Intenzioni3D` scritti, per esempio `(Λa = 3.0,)`;
quelli non scritti tornano a riposo. `stato`: le voci dello stato interno scritte; le altre
restano com'erano. `dichiara`: il messaggio di questa passata alle vicine. `secrezioni` e
`deposita`: quantità/min per sostanza e per matrice. `richiesta`: `:dividiti` (con `piano`,
l'angolo in rad della linea di divisione nel piano del foglio, `NaN` per la regola di
Hertwig), `:muori`, `:nessuna`, o `nothing` se non scritta.
"""
struct Uscita
    intenzioni::Dict{Symbol,Float64}
    stato::Dict{Symbol,Float64}
    dichiara::Dict{Symbol,Float64}
    secrezioni::Dict{Symbol,Float64}
    deposita::Dict{Symbol,Float64}
    richiesta::Union{Nothing,Symbol}
    piano::Float64
end

function Uscita(; intenzioni = NamedTuple(), stato = Dict{Symbol,Float64}(),
                dichiara = Dict{Symbol,Float64}(), secrezioni = Dict{Symbol,Float64}(),
                deposita = Dict{Symbol,Float64}(), richiesta = nothing, piano = NaN)
    ints = Dict{Symbol,Float64}(Symbol(k) => Float64(v) for (k, v) in pairs(intenzioni))
    all(in(fieldnames(Intenzioni3D)), keys(ints)) ||
        throw(ArgumentError("intenzioni sconosciute: $(setdiff(keys(ints), fieldnames(Intenzioni3D)))"))
    (richiesta === nothing || richiesta in (:nessuna, :dividiti, :muori)) ||
        throw(ArgumentError("richiesta sconosciuta: $richiesta"))
    return Uscita(ints, Dict{Symbol,Float64}(stato), Dict{Symbol,Float64}(dichiara),
                  Dict{Symbol,Float64}(secrezioni), Dict{Symbol,Float64}(deposita), richiesta,
                  Float64(piano))
end

# Le intenzioni di una cellula: quelle a riposo, con i campi che l'uscita scrive.
_applica(riposo::Intenzioni3D, u::Uscita) =
    Intenzioni3D((get(u.intenzioni, f, getfield(riposo, f)) for f in fieldnames(Intenzioni3D))...)

_richiede(u::Uscita) = !(u.richiesta === nothing || u.richiesta == :nessuna)

"""
    Fase(nome, regola, parametri, fine; scrive = nothing)

Una riga della tabella di fasi: `regola(p, parametri)` → `Uscita` e `fine(p, parametri)` →
Bool, funzioni pure della percezione `p`. La fase dichiara che cosa la regola scrive:
`scrive`, fra le intenzioni (:σa, :σb, :σl, :Λa, :V0), :secrezioni, :deposita, :richiesta;
`stato`, le voci dello stato interno; `dichiara`, le voci del messaggio alle vicine. Una fase
dichiarata si compone con `insieme`, che rifiuta i conflitti quando compone (per ogni cellula,
qualunque sia la corsa), e la passata confronta l'uscita, che è un delta, con la
dichiarazione: una voce scritta e non dichiarata è un errore. Una fase non dichiarata
(`scrive = nothing`) può scrivere tutto ma non si compone.
"""
struct Fase{R,F,P<:NamedTuple}
    nome::Symbol
    regola::R
    parametri::P
    fine::F
    scrive::Union{Nothing,Tuple{Vararg{Symbol}}}
    scrive_stato::Union{Nothing,Tuple{Vararg{Symbol}}}
    scrive_dichiara::Union{Nothing,Tuple{Vararg{Symbol}}}
end

const SCRIVIBILI = (:σa, :σb, :σl, :Λa, :V0, :secrezioni, :deposita, :richiesta)

function Fase(nome, regola, parametri, fine; scrive = nothing, stato = (), dichiara = ())
    scrive === nothing && return Fase(Symbol(nome), regola, parametri, fine, nothing, nothing, nothing)
    d = Tuple(Symbol.(scrive))
    all(in(SCRIVIBILI), d) ||
        throw(ArgumentError("la fase $nome dichiara $(setdiff(d, SCRIVIBILI)), che non si scrive"))
    return Fase(Symbol(nome), regola, parametri, fine, d, Tuple(Symbol.(stato)),
                Tuple(Symbol.(dichiara)))
end

# La regola ha scritto soltanto ciò che la sua fase dichiara? L'uscita è un delta: ogni voce
# presente è una scrittura, anche se vale quanto prima.
function _controlla_dichiarazione(f::Fase, p::Percezione, u::Uscita)
    f.scrive === nothing && return nothing
    fuori(cosa) = throw(ArgumentError("la fase $(f.nome) scrive $cosa senza dichiararlo"))
    for campo in keys(u.intenzioni)
        campo in f.scrive || fuori(campo)
    end
    (:secrezioni in f.scrive || isempty(u.secrezioni)) || fuori(:secrezioni)
    (:deposita in f.scrive || isempty(u.deposita)) || fuori(:deposita)
    (:richiesta in f.scrive || u.richiesta === nothing) || fuori(:richiesta)
    for k in keys(u.stato)
        k in f.scrive_stato || fuori("stato[:$k]")
    end
    for k in keys(u.dichiara)
        k in f.scrive_dichiara || fuori("dichiara[:$k]")
    end
    return nothing
end

"""
    Interruttori(; chimica = true, onda = true, meccanica = true, topologia = true,
                   popolazione = true, matrice = true)

Un interruttore per meccanismo: spento, il passo corrispondente della passata non si fa
(`onda` spenta: si avanza di fase solo per la propria fine).
"""
Base.@kwdef struct Interruttori
    chimica::Bool = true
    onda::Bool = true
    meccanica::Bool = true
    topologia::Bool = true
    popolazione::Bool = true
    matrice::Bool = true
end

"""
    StatoCellula

Ciò che il motore tiene per ogni cellula: la sua copia del programma, l'identità (non cambia
quando le cellule si rinumerano), la fase e la passata in cui vi è entrata (`ingresso`), la
passata di nascita, se è in apoptosi, lo stato interno, ciò che dichiara, le intenzioni a
riposo e quelle correnti, le secrezioni correnti.
"""
mutable struct StatoCellula
    programma::Tuple
    identita::Int
    morente::Bool
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
    Embrione(t, programma; parametri, dt, campi = Dict(), matrici = Dict(), stati = nothing,
             interruttori = Interruttori(), seme = 0, ogni = 1, topologia = :generale,
             volume_apoptosi = 1.0, tempo_apoptosi = 5.0, tolleranza = 1e-8, solutore = :pieno,
             operazioni = (soglia = 0.3, area_soglia = 0.3, lunghezza_nuova = 1.0))

Un tessuto con un programma copiato in ogni cellula, i campi chimici e le matrici per nome, i
parametri meccanici e il passo meccanico `dt` (min). `stati[c]` è lo stato interno iniziale
della cellula c (vuoto se non dato). Le intenzioni a riposo sono quelle del tessuto. `seme`
fissa la casualità; la popolazione si aggiorna ogni `ogni` passate; `topologia` è `:generale`
(11_operazioni.jl) o `:prismi` (scambi di 08 ed estrusioni di 09: serve per dividersi e
morire). Una cellula in apoptosi ha il volume bersaglio che scende esponenzialmente, con tempo
`tempo_apoptosi` (min, libero), fino a `volume_apoptosi` (µm³, libero). `solutore` è il metodo
del passo meccanico: `:pieno` (`minimizza!`, 06_integratore.jl) o `:sparso`
(`minimizza_sparso!`, 15_sparso.jl, per i tessuti di migliaia di cellule).
"""
mutable struct Embrione
    tessuto::Tessuto3D
    cellule::Vector{StatoCellula}
    campi::Dict{Symbol,Campo}
    matrici::Dict{Symbol,Campo}
    parametri::ParametriMeccanici
    dt::Float64
    passata::Int
    interruttori::Interruttori
    seme::Int
    ogni::Int
    topologia::Symbol
    volume_apoptosi::Float64
    tempo_apoptosi::Float64
    tolleranza::Float64
    solutore::Symbol
    operazioni::NamedTuple
    prossima_identita::Int
    registro::Vector{NamedTuple{(:passata, :fasi),Tuple{Int,Vector{Int}}}}
    eventi::Vector{NamedTuple{(:passata, :tipo, :identita),Tuple{Int,Symbol,Int}}}
end

function Embrione(t::Tessuto3D, programma; parametri, dt, campi = Dict{Symbol,Campo}(),
                  matrici = Dict{Symbol,Campo}(), stati = nothing, interruttori = Interruttori(),
                  seme = 0, ogni = 1, topologia = :generale, volume_apoptosi = 1.0,
                  tempo_apoptosi = 5.0, tolleranza = 1e-8, solutore = :pieno,
                  operazioni = (soglia = 0.3, area_soglia = 0.3, lunghezza_nuova = 1.0))
    tabella = Tuple(programma)
    isempty(tabella) && throw(ArgumentError("il programma non ha fasi"))
    topologia in (:generale, :prismi) || throw(ArgumentError("topologia :generale o :prismi"))
    ogni >= 1 || throw(ArgumentError("ogni deve essere almeno 1"))
    solutore in (:pieno, :sparso) || throw(ArgumentError("solutore :pieno o :sparso"))
    cellule = [StatoCellula(tabella, c, false, 1, 1, 1,
                            stati === nothing ? Dict{Symbol,Float64}() : Dict{Symbol,Float64}(stati[c]),
                            Dict{Symbol,Float64}(), t.intenzioni[c], t.intenzioni[c],
                            Dict{Symbol,Float64}()) for c in eachindex(t.cellule)]
    return Embrione(t, cellule, Dict{Symbol,Campo}(campi), Dict{Symbol,Campo}(matrici), parametri,
                    Float64(dt), 0, interruttori, seme, ogni, topologia, Float64(volume_apoptosi),
                    Float64(tempo_apoptosi), tolleranza, solutore, operazioni, length(cellule) + 1,
                    NamedTuple{(:passata, :fasi),Tuple{Int,Vector{Int}}}[],
                    NamedTuple{(:passata, :tipo, :identita),Tuple{Int,Symbol,Int}}[])
end

# Un numero in [0, 1) da (seme, identità, passata): SplitMix64, stabile fra versioni di Julia.
function _mescola(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end
_caso(seme, identita, passata) =
    (_mescola(_mescola(_mescola(seme % UInt64) ⊻ (identita % UInt64)) ⊻ (passata % UInt64)) >> 11) * 2.0^-53

_con_volume(r::Intenzioni3D, V0) = Intenzioni3D(r.σa, r.σb, r.σl, r.Λa, V0)

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
    base = isempty(e.matrici) ? (0.0, 0.0, 0.0) : centro_basale(t, c)
    matrice = Dict{Symbol,Float64}(nome => concentrazione(campo, base) for (nome, campo) in e.matrici)
    return Percezione(s.fase, e.passata - s.ingresso, e.passata - s.nascita,
                      _caso(e.seme, s.identita, e.passata), copy(s.stato), s.riposo, misure,
                      [Vicina(e.cellule[v].fase, copy(e.cellule[v].dichiara)) for v in vicine],
                      matrice, concentrazioni)
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
        if s.morente                    # in apoptosi il programma è fermo e il volume scende
            V0 = max(e.volume_apoptosi, s.intenzioni.V0 * exp(-e.dt / e.tempo_apoptosi))
            r = s.intenzioni
            uscite[c] = Uscita(; intenzioni = (σa = r.σa, σb = r.σb, σl = r.σl, Λa = r.Λa, V0 = V0))
            continue
        end
        p = _percezione(e, c, vicine[c])
        fase = s.programma[s.fase]
        uscite[c] = fase.regola(p, fase.parametri)
        _controlla_dichiarazione(fase, p, uscite[c])
        if s.fase < length(s.programma)
            avanza[c] = fase.fine(p, fase.parametri) ||
                        (e.interruttori.onda && any(v.fase > s.fase for v in p.vicine))
        end
    end
    for (c, s) in enumerate(e.cellule)
        u = uscite[c]
        merge!(s.stato, u.stato); s.dichiara = copy(u.dichiara)
        s.intenzioni = _applica(s.riposo, u); s.secrezioni = copy(u.secrezioni)
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
        # il passo implicito di 06_integratore.jl, senza fermarsi se non converge: succede
        # quando la topologia è in ritardo (una cellula in apoptosi che stringe i lati sotto
        # la lunghezza minima prima degli scambi); come in rilassa_con_topologia!, si va avanti
        # e la topologia, al passo 6, sblocca. Resta registrato.
        minimo! = e.solutore == :sparso ? minimizza_sparso! : minimizza!
        r = minimo!(t, e.parametri; tolleranza = e.tolleranza, iterazioni_massime = 200,
                    ancora = copy(t.x), peso = e.parametri.η / e.dt)
        r.convergito || push!(e.eventi, (passata = e.passata, tipo = :meccanica_non_risolta,
                                         identita = 0))
    end
    # 6. topologia
    if e.interruttori.topologia
        if e.topologia == :generale
            _operazioni_dovute!(t; e.operazioni..., monostrato = true)
        else
            _uno_scambio!(t; soglia = e.operazioni.soglia,
                          lunghezza_nuova = e.operazioni.lunghezza_nuova)
            for r in estrusioni!(t; area_soglia = e.operazioni.area_soglia)
                push!(e.eventi, (passata = e.passata, tipo = :estrusione,
                                 identita = e.cellule[r.cellula].identita))
                deleteat!(e.cellule, r.cellula); deleteat!(uscite, r.cellula)
            end
        end
    end
    # 7. popolazione
    if e.interruttori.popolazione && e.passata % e.ogni == 0
        _popolazione!(e, uscite)
    end
    # 8. matrice
    if e.interruttori.matrice
        for (nome, campo) in e.matrici
            scambi = [Scambio(centro_basale(t, c), get(uscite[c].deposita, nome, 0.0), 0.0)
                      for c in eachindex(e.cellule) if get(uscite[c].deposita, nome, 0.0) != 0]
            avanza!(campo, e.dt, scambi)
        end
    end
    # 9. registro dell'osservatore
    push!(e.registro, (passata = e.passata, fasi = [s.fase for s in e.cellule]))
    return e
end

# Passo 7: le richieste accolte. Prima le morti, poi le divisioni (le figlie vanno in fondo, così
# gli indici delle madri non cambiano).
function _popolazione!(e::Embrione, uscite)
    t = e.tessuto
    richieste = [c for c in eachindex(e.cellule)
                 if _richiede(uscite[c]) && !e.cellule[c].morente]
    isempty(richieste) && return e
    e.topologia == :prismi ||
        throw(ArgumentError("dividersi e morire chiedono il guscio a prismi (topologia = :prismi)"))
    for c in richieste
        uscite[c].richiesta == :muori || continue
        s = e.cellule[c]
        s.morente = true
        s.riposo = _con_volume(s.riposo, e.volume_apoptosi)
        push!(e.eventi, (passata = e.passata, tipo = :apoptosi, identita = s.identita))
    end
    for c in richieste
        uscite[c].richiesta == :dividiti || continue
        s = e.cellule[c]
        t.intenzioni[c] = s.intenzioni
        f = dividi!(t, c; angolo = uscite[c].piano)
        if f === nothing
            push!(e.eventi, (passata = e.passata, tipo = :divisione_rifiutata, identita = s.identita))
            continue
        end
        s.riposo = _con_volume(s.riposo, s.riposo.V0 / 2)
        s.intenzioni = _con_volume(s.intenzioni, s.intenzioni.V0 / 2)
        figlia = StatoCellula(s.programma, e.prossima_identita, false, s.fase, s.ingresso,
                              e.passata + 1, copy(s.stato), copy(s.dichiara), s.riposo,
                              s.intenzioni, copy(s.secrezioni))
        e.prossima_identita += 1
        push!(e.cellule, figlia)
        t.intenzioni[c] = s.intenzioni; t.intenzioni[f] = figlia.intenzioni
        push!(uscite, Uscita())
        push!(e.eventi, (passata = e.passata, tipo = :divisione, identita = figlia.identita))
    end
    return e
end

# Nel guscio a prismi: lo scambio di vicine sulla parete più corta, se sta sotto soglia in alto
# e in basso (come `scambi_di_vicine!` di 08_scambi.jl, ma uno solo).
function _uno_scambio!(t::Tessuto3D; soglia, lunghezza_nuova)
    W = poligoni(t)
    corto = nothing
    Lmin = soglia
    for w in W, k in eachindex(w)
        (j1, i1), (j2, i2) = w[k], w[mod1(k + 1, length(w))]
        La = _norma(_meno(_pos(t, t.x, 2j2 - 1, i2), _pos(t, t.x, 2j1 - 1, i1)))
        Lb = _norma(_meno(_pos(t, t.x, 2j2, i2), _pos(t, t.x, 2j1, i1)))
        L = max(La, Lb)
        if L < Lmin && _quartetto(t, W, j1, j2) !== nothing
            corto, Lmin = (j1, j2), L
        end
    end
    corto === nothing && return nothing
    return scambia!(t, corto...; lunghezza_nuova)
end
