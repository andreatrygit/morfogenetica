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
# Ogni prova è una funzione che restituisce un Bool. Commit rosso: marcate @test_broken finché
# il codice non c'è.

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
    regola(p, par) = Uscita(copy(p.stato), Intenzioni3D(p.intenzioni.σa, p.intenzioni.σb,
                                                         p.intenzioni.σl, 2.0, p.intenzioni.V0),
                            Dict{Symbol,Float64}(), Dict{Symbol,Float64}())
    fai(scrive) = Embrione(foglio_programma(1), [Fase(:f, regola, NamedTuple(), (p, par) -> false; scrive)];
                           parametri = P_PROGRAMMA, dt = 1.0, interruttori = Interruttori(; SPENTO...))
    passata!(fai((:Λa,)))
    return try
        passata!(fai((:V0,))); false
    catch err
        err isa ArgumentError
    end
end

@testset "composizione: il conflitto si vede alla composizione" begin
    @test_broken prova_conflitto_alla_composizione()
    @test_broken prova_insieme_dichiara_l_unione()
    @test_broken prova_non_dichiarata_non_si_compone()
    @test_broken prova_la_passata_controlla_le_dichiarazioni()
end
