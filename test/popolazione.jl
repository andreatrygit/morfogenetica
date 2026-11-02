# R4, chiusura (docs/roadmap.md): casualità con seme, matrice, divisione e morte (specifica,
# sezione 6, passi 7 e 8; repertorio, sezione 7).
#
#   p.caso                          un numero in [0, 1) per cellula e per passata, dal seme
#                                   dell'embrione e dall'identità della cellula
#   Uscita(stato, intenzioni, secrezioni, dichiara; deposita, richiesta, piano)
#                                   richiesta :nessuna, :dividiti o :muori; piano è l'angolo
#                                   (rad, nel piano del foglio) della linea di divisione, NaN
#                                   per la regola di Hertwig (perpendicolare all'asse lungo)
#   Embrione(...; seme, ogni, topologia, matrici, volume_apoptosi)
#                                   `ogni`: la popolazione si aggiorna ogni `ogni` passate;
#                                   `topologia`: :generale (11_operazioni.jl) o :prismi
#                                   (08 scambi e 09 estrusioni, il guscio a prismi)
#   Interruttori(; ..., popolazione, matrice)
#   e.cellule[c].identita, .nascita, .morente;  e.eventi (passata, tipo, identita)
#   dividi!(t, c; angolo)           divide una cellula del guscio a prismi con un piano
#                                   verticale per il centro della faccia apicale; restituisce
#                                   l'indice della figlia, o `nothing`
#   centro_basale(t, c)             dove la cellula depone e legge la matrice
#
# Ogni prova è una funzione che restituisce un Bool. Commit rosso: marcate @test_broken finché
# il codice (motore/13_programma.jl, motore/14_popolazione.jl) non c'è.

volume_tessuto(t) = sum(volume_cellula(t, c) for c in eachindex(t.cellule))

function prova_dividi_conserva_e_rispetta_gli_invarianti()
    for t in (foglio_programma(2), foglio_periodico(fill(PRISMA, 16), 4, 4; lato = 5.0, altezza = 15.0))
        V = volume_tessuto(t); V1 = volume_cellula(t, 1); n = length(t.cellule)
        f = dividi!(t, 1; angolo = 0.3)
        f == n + 1 || return false
        isapprox(volume_tessuto(t), V; rtol = 1e-12) || return false
        # un esagono regolare tagliato per il centro: due metà uguali
        isapprox(volume_cellula(t, 1), V1 / 2; rtol = 1e-12) || return false
        isapprox(volume_cellula(t, f), V1 / 2; rtol = 1e-12) || return false
        isempty(problemi_topologici(t)) || return false
        f in vicine_di(t, 1) || return false
        errore_gradiente(t, P_PROGRAMMA) < 1e-7 || return false
    end
    return true
end

function prova_dividi_secondo_hertwig()
    # una cellula allungata in x si divide con una linea lungo y: le figlie stanno ai due lati
    # in x
    t = foglio_programma(1)
    for k in 1:length(t.x) ÷ 3
        t.x[3k-2] *= 1.6; t.ancore[3k-2] *= 1.6
    end
    f = dividi!(t, 1; angolo = NaN)
    f === nothing && return false
    a = centro_cellula(t, 1); b = centro_cellula(t, f)
    return abs(a[1] - b[1]) > 5 * abs(a[2] - b[2])
end

regola_linea(; richiesta = :dividiti) =
    (p, par) -> Uscita(copy(p.stato), p.intenzioni, Dict{Symbol,Float64}(), Dict{Symbol,Float64}();
                       richiesta = get(p.stato, :linea, 0.0) == 1 ? richiesta : :nessuna)

function embrione_linea(; ogni = 1, popolazione = true, meccanica = false, topologia = :prismi,
                        richiesta = :dividiti, t = foglio_programma(2))
    stati = [Dict{Symbol,Float64}(:linea => (c == 1 ? 1.0 : 0.0)) for c in eachindex(t.cellule)]
    return Embrione(t, [Fase(:f, regola_linea(; richiesta), NamedTuple(), (p, par) -> false)];
                    parametri = P_PROGRAMMA, dt = 1.0, stati, ogni, topologia,
                    interruttori = Interruttori(; chimica = false, onda = true, meccanica,
                                                topologia = meccanica, popolazione))
end

function prova_la_figlia_eredita()
    e = embrione_linea()
    n = length(e.cellule)
    V0 = e.cellule[1].riposo.V0
    e.cellule[1].stato[:altro] = 7.0
    passata!(e)
    length(e.cellule) == n + 1 || return false
    m = e.cellule[1]; f = e.cellule[end]
    return f.fase == m.fase && f.programma === m.programma && f.stato == m.stato &&
           f.identita != m.identita && f.nascita == 2 &&
           m.riposo.V0 == V0 / 2 && f.riposo.V0 == V0 / 2 &&
           length(unique(c.identita for c in e.cellule)) == length(e.cellule)
end

function prova_popolazione_ogni_m_passate()
    # la linea chiede di dividersi a ogni passata; si divide solo alle passate 3 e 6
    e = embrione_linea(; ogni = 3)
    n = length(e.cellule)
    conti = Int[]
    for _ in 1:7
        passata!(e)
        push!(conti, length(e.cellule))
    end
    conti == [n, n, n + 1, n + 1, n + 1, n + 3, n + 3] || return false
    # interruttore: spenta la popolazione, nessuno si divide
    s = embrione_linea(; popolazione = false)
    for _ in 1:3
        passata!(s)
    end
    return length(s.cellule) == n
end

function prova_la_morte_e_un_estrusione()
    # la cellula centrale chiede di morire: si restringe, diventa un triangolo con gli scambi di
    # vicine ed esce dal foglio; le altre restano, con le loro identità
    e = embrione_linea(; richiesta = :muori, meccanica = true)
    n = length(e.cellule)
    morta = e.cellule[1].identita
    for _ in 1:40
        passata!(e)
        length(e.cellule) < n && break
    end
    length(e.cellule) == n - 1 || return false
    any(ev -> ev.tipo == :estrusione && ev.identita == morta, e.eventi) || return false
    return !(morta in [c.identita for c in e.cellule]) && isempty(problemi_topologici(e.tessuto))
end

function prova_la_popolazione_chiede_i_prismi()
    e = embrione_linea(; topologia = :generale)
    return try
        passata!(e); false
    catch err
        err isa ArgumentError
    end
end

function prova_la_matrice_persiste()
    # la cellula centrale depone per 5 passate; poi smette, ma ciò che ha deposto resta lì
    t = foglio_programma(2)
    matrice = Campo(Sostanza(:collagene, 0.0, 0.0); origine = (-30.0, -30.0, -3.0), h = 2.0,
                    n = (30, 30, 4))
    function regola(p, par)
        q = get(p.stato, :linea, 0.0) == 1 && p.eta < 5 ? 1.0 : 0.0
        return Uscita(merge(p.stato, Dict(:letta => get(p.matrice, :collagene, 0.0))), p.intenzioni,
                      Dict{Symbol,Float64}(), Dict{Symbol,Float64}(); deposita = Dict(:collagene => q))
    end
    stati = [Dict{Symbol,Float64}(:linea => (c == 1 ? 1.0 : 0.0)) for c in eachindex(t.cellule)]
    e = Embrione(t, [Fase(:f, regola, NamedTuple(), (p, par) -> false)]; parametri = P_PROGRAMMA,
                 dt = 1.0, stati, matrici = Dict(:collagene => matrice),
                 interruttori = Interruttori(; SPENTO...))
    for _ in 1:5
        passata!(e)
    end
    M = massa(e.matrici[:collagene])
    isapprox(M, 5.0; rtol = 1e-12) || return false
    d = concentrazione(e.matrici[:collagene], centro_basale(e.tessuto, 1))
    for _ in 1:5
        passata!(e)
    end
    return massa(e.matrici[:collagene]) == M &&
           concentrazione(e.matrici[:collagene], centro_basale(e.tessuto, 1)) == d && d > 0 &&
           e.cellule[1].stato[:letta] == d && e.cellule[19].stato[:letta] == 0
end

function caso_registrato(; seme, ordine = nothing)
    t = foglio_programma(2)
    regola(p, par) = uscita_ferma(p; stato = Dict(:caso => p.caso))
    e = Embrione(t, [Fase(:f, regola, NamedTuple(), (p, par) -> false)]; parametri = P_PROGRAMMA,
                 dt = 1.0, seme, interruttori = Interruttori(; SPENTO...))
    valori = Float64[]
    for _ in 1:3
        ordine === nothing ? passata!(e) : passata!(e; ordine)
        append!(valori, [c.stato[:caso] for c in e.cellule])
    end
    return valori
end

function prova_casualita_con_seme()
    a = caso_registrato(; seme = 1)
    b = caso_registrato(; seme = 1, ordine = reverse(1:19))
    c = caso_registrato(; seme = 2)
    (a == b && a != c && all(0 .<= a .< 1) && length(unique(a)) == length(a)) || return false
    # uniforme: media di molti valori vicina a 1/2
    v = reduce(vcat, caso_registrato(; seme = s) for s in 1:200)
    return abs(sum(v) / length(v) - 0.5) < 0.01
end

@testset "R4 popolazione, matrice, caso" begin
    @test_broken prova_dividi_conserva_e_rispetta_gli_invarianti()
    @test_broken prova_dividi_secondo_hertwig()
    @test_broken prova_la_figlia_eredita()
    @test_broken prova_popolazione_ogni_m_passate()
    @test_broken prova_la_morte_e_un_estrusione()
    @test_broken prova_la_popolazione_chiede_i_prismi()
    @test_broken prova_la_matrice_persiste()
    @test_broken prova_casualita_con_seme()
end
