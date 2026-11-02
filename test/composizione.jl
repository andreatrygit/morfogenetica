# Il conflitto fra regole si vede alla composizione, non alla corsa (docs/roadmap.md, R5).
#
# Ogni fase dichiara che cosa scrive: le intenzioni (:σa, :σb, :σl, :Λa, :V0), :secrezioni,
# :deposita, :richiesta.
#
#   Fase(nome, regola, parametri, fine; scrive)   `scrive` è una tupla di simboli; senza, la
#                                                 fase non è dichiarata (può scrivere tutto)
#   f.scrive                                      la dichiarazione, `nothing` se manca
#   insieme(a, b; nome)                           rifiuta, quando compone, due fasi che
#                                                 scrivono la stessa intenzione o la stessa
#                                                 richiesta, e le fasi non dichiarate;
#                                                 secrezioni e deposizioni si sommano
#   passata!                                      rifiuta una regola che scrive fuori da ciò che
#                                                 la sua fase dichiara
#
# Il caso che ha fatto nascere la prova: `costringi(Λ = 0)` ‖ `costringi(Λ = 3)`. Confrontando
# i valori con quelli a riposo, la regola che scrive Λa = 0 sembrava non toccare niente, e
# vinceva in silenzio Λ = 3. Una scrittura del valore a riposo e un "non tocco" non si
# distinguono dai risultati: si distinguono solo dichiarandoli.
#
# Ogni prova è una funzione che restituisce un Bool. Scritte prima del codice (commit rosso, con
# @test_broken), passano dal commit verde. Le prove della seconda parte, più sotto, hanno
# seguito lo stesso percorso.

function prova_conflitto_alla_composizione()
    a = costringi(; Λ = 0.0, sostanza = :s, soglia = 0.0)
    b = costringi(; Λ = 3.0, sostanza = :s, soglia = 0.0)
    return try
        insieme(a, b; nome = :due); false
    catch err
        err isa ArgumentError
    end
end

function prova_insieme_dichiara_l_unione()
    g = insieme(secerni(; sostanza = :s, tasso = 1.0), costringi(; Λ = 1.0, sostanza = :s, soglia = 0.0);
                nome = :g)
    h = con_fine(g, (p, par) -> true)
    return costringi(; Λ = 1.0, sostanza = :s, soglia = 0.0).scrive == (:Λa,) &&
           secerni(; sostanza = :s, tasso = 1.0).scrive == (:secrezioni,) &&
           Set(g.scrive) == Set((:secrezioni, :Λa)) && h.scrive == g.scrive
end

function prova_non_dichiarata_non_si_compone()
    libera = Fase(:libera, (p, par) -> uscita_ferma(p), NamedTuple(), (p, par) -> false)
    libera.scrive === nothing || return false
    return try
        insieme(libera, secerni(; sostanza = :s, tasso = 1.0); nome = :x); false
    catch err
        err isa ArgumentError
    end
end

function prova_la_passata_controlla_le_dichiarazioni()
    # una regola che scrive Λa: dichiarata passa, non dichiarata no
    regola(p, par) = Uscita(; intenzioni = (Λa = 2.0,))
    fai(scrive) = Embrione(foglio_programma(1), [Fase(:f, regola, NamedTuple(), (p, par) -> false; scrive)];
                           parametri = P_PROGRAMMA, dt = 1.0, interruttori = Interruttori(; SPENTO...))
    passata!(fai((:Λa,)))
    return try
        passata!(fai((:V0,))); false
    catch err
        err isa ArgumentError
    end
end

# --- Seconda parte (2 novembre): uscite a delta e dichiarazioni di stato e vicine ----------------
#
#   Uscita(; intenzioni, stato, dichiara, secrezioni, deposita, richiesta, piano)
#       l'uscita è un delta: contiene solo ciò che la regola scrive. `intenzioni` sono i campi
#       scritti (per esempio (Λa = 3.0,)); ciò che non è scritto resta: lo stato interno quello
#       di prima, le intenzioni quelle a riposo; la dichiarazione alle vicine è il messaggio di
#       questa passata. La presenza di una voce è il registro delle scritture: una scrittura del
#       valore a riposo si distingue da un "non tocco".
#   Fase(...; scrive, stato, dichiara)
#       le voci dello stato interno e della dichiarazione che la regola scrive
#
# Prima: lo stato si univa per valori, e ogni regola restituiva lo stato intero; una regola che
# aggiorna un contatore, composta con una che lo ricopia, dava un conflitto falso.

regola_contatore(p, par) = Uscita(; stato = Dict(:n => get(p.stato, :n, 0.0) + 1))
contatore() = Fase(:contatore, regola_contatore, NamedTuple(), (p, par) -> false; scrive = (),
                   stato = (:n,))

embrione_di(fase; stati = nothing) =
    Embrione(foglio_programma(1), [fase]; parametri = P_PROGRAMMA, dt = 1.0, stati,
             interruttori = Interruttori(; SPENTO...))

function prova_stato_senza_falsi_conflitti()
    f = insieme(contatore(), costringi(; Λ = 1.0, sostanza = :s, soglia = 0.0); nome = :cc)
    e = embrione_di(f; stati = [Dict(:n => 0.0, :altro => 7.0) for _ in 1:7])
    passata!(e); passata!(e)
    return all(c.stato == Dict(:n => 2.0, :altro => 7.0) for c in e.cellule)
end

function prova_conflitto_su_stato_e_dichiara_alla_composizione()
    rifiuta(a, b) = try
        insieme(a, b; nome = :x); false
    catch err
        err isa ArgumentError
    end
    voce(d) = Fase(:voce, (p, par) -> Uscita(; dichiara = Dict(:v => 1.0)), NamedTuple(),
                   (p, par) -> false; scrive = (), dichiara = d)
    return rifiuta(contatore(), contatore()) && rifiuta(voce((:v,)), voce((:v,))) &&
           !rifiuta(contatore(), voce((:v,)))
end

function prova_la_scrittura_del_valore_a_riposo_si_vede()
    # Λa scritto uguale a quello a riposo, senza dichiararlo: prima passava, ora no
    f = Fase(:muta, (p, par) -> Uscita(; intenzioni = (Λa = p.intenzioni.Λa,)), NamedTuple(),
             (p, par) -> false; scrive = ())
    return try
        passata!(embrione_di(f)); false
    catch err
        err isa ArgumentError
    end
end

function prova_lo_stato_non_dichiarato_non_si_scrive()
    f = Fase(:muta, (p, par) -> Uscita(; stato = Dict(:k => 1.0)), NamedTuple(),
             (p, par) -> false; scrive = ())
    return try
        passata!(embrione_di(f)); false
    catch err
        err isa ArgumentError
    end
end

function prova_lo_stato_e_un_delta()
    # chi non scrive lo stato lo lascia com'è; chi scrive una voce non cancella le altre
    e = embrione_di(contatore(); stati = [Dict(:n => 0.0, :b => 5.0) for _ in 1:7])
    passata!(e)
    s = embrione_di(costringi(; Λ = 1.0, sostanza = :s, soglia = 0.0);
                    stati = [Dict(:b => 5.0) for _ in 1:7])
    passata!(s)
    return all(c.stato == Dict(:n => 1.0, :b => 5.0) for c in e.cellule) &&
           all(c.stato == Dict(:b => 5.0) for c in s.cellule)
end

@testset "composizione: il conflitto si vede alla composizione" begin
    @test prova_conflitto_alla_composizione()
    @test prova_insieme_dichiara_l_unione()
    @test prova_non_dichiarata_non_si_compone()
    @test prova_la_passata_controlla_le_dichiarazioni()
    @test prova_stato_senza_falsi_conflitti()
    @test prova_conflitto_su_stato_e_dichiara_alla_composizione()
    @test prova_la_scrittura_del_valore_a_riposo_si_vede()
    @test prova_lo_stato_non_dichiarato_non_si_scrive()
    @test prova_lo_stato_e_un_delta()
end
