# R5 (docs/roadmap.md): l'organizzatore che nasce senza seme, per attivazione locale e
# inibizione a lungo raggio (Gierer e Meinhardt), con i numeri misurati di Nodal e Lefty.
#
#   attiva_inibisci(; attivatore, inibitore, ka, kh, κ, rumore, c0)
#       la primitiva (primitive/attiva_inibisci.jl): ogni cellula produce l'attivatore in modo
#       autocatalitico e saturante, frenato dall'inibitore, e l'inibitore in proporzione al
#       quadrato dell'attivatore; la sua espressione ha un rumore (dal seme) di ampiezza `rumore`.
#       Le produzioni sono tali che lo stato uniforme sia a = h = c0.
#
# Parametri (Müller et al., Science 336, 721, 2012, doi:10.1126/science.1221920): Cyclops
# (attivatore) D = 0.7 µm²/s = 42 µm²/min, k = 1.22e-4/s = 0.0073/min; Lefty2 (inibitore)
# D = 18.9 µm²/s = 1134 µm²/min, k = 0.69e-4/s = 0.0041/min. Liberi: κ = 0.29 (saturazione
# dell'autoattivazione: senza, poiché Lefty si elimina più lentamente di Nodal, lo stato uniforme
# oscilla; 0.29 massimizza la crescita), rumore 20%, c0. Con questi numeri l'instabilità ha
# lunghezze d'onda fra 0.7 e 1.9 mm (massimo 1.13 mm) e cresce al più di un fattore e ogni 13
# ore: una sola macchia nasce in un tessuto di circa 1 mm, in giorni. Cellule di 15 µm di lato
# (26 µm di diametro, come le cellule profonde di una blastula).
#
# Ogni prova è una funzione che restituisce un Bool. Scritte prima del codice (commit rosso, con
# @test_broken), passano dal commit verde.

const NODAL_LEFTY = (Da = 42.0, ka = 0.0073, Dh = 1134.0, kh = 0.0041, κ = 0.29)

function embrione_organizzatore(; seme = 1, rumore = 0.2, Dh = NODAL_LEFTY.Dh)
    t = foglio_periodico(fill(Intenzioni3D(1.0, 1.0, 1.0, 0.0, 8768.0), 38 * 44), 38, 44;
                         lato = 15.0, altezza = 15.0)
    Lx, Ly = t.x[end-1], t.x[end]
    griglia(nome, D, k) = Campo(Sostanza(nome, D, k); origine = (0.0, 0.0, 0.0),
                                h = (Lx / 82, Ly / 82, 15.0), n = (82, 82, 1),
                                bordi = (bordo_periodico, bordo_periodico, bordo_chiuso))
    c0 = 1.0
    campi = Dict(:nodal => griglia(:nodal, NODAL_LEFTY.Da, NODAL_LEFTY.ka),
                 :lefty => griglia(:lefty, Dh, NODAL_LEFTY.kh))
    for campo in values(campi)
        campo.c .= c0
    end
    f = attiva_inibisci(; attivatore = :nodal, inibitore = :lefty, ka = NODAL_LEFTY.ka,
                        kh = NODAL_LEFTY.kh, κ = NODAL_LEFTY.κ, rumore, c0)
    return Embrione(t, [f]; parametri = P_PROGRAMMA, dt = 30.0, campi, seme,
                    interruttori = Interruttori(; chimica = true, meccanica = false,
                                                topologia = false, popolazione = false))
end

attivatore(e) = [concentrazione(e.campi[:nodal], centro_cellula(e.tessuto, c)) for c in eachindex(e.cellule)]

# Le macchie: componenti connesse (per contatto) delle cellule sopra la metà fra minimo e massimo.
function macchie(e)
    a = attivatore(e)
    s = (maximum(a) + minimum(a)) / 2
    sopra = Set(c for c in eachindex(a) if a[c] > s)
    n = 0
    while !isempty(sopra)
        n += 1
        coda = [pop!(sopra)]
        while !isempty(coda)
            c = pop!(coda)
            for v in vicine_di(e.tessuto, c)
                v in sopra && (delete!(sopra, v); push!(coda, v))
            end
        end
    end
    return n
end

contrasto(e) = (a = attivatore(e); maximum(a) / minimum(a))

const _ORGANIZZATORI = Dict{Any,Any}()
function corsa_organizzatore(; kw...)
    return get!(_ORGANIZZATORI, sort(collect(kw))) do
        e = embrione_organizzatore(; kw...)
        contrasti = Float64[]
        for _ in 1:200                                # 100 ore
            passata!(e)
            push!(contrasti, contrasto(e))
        end
        (e, contrasti)
    end
end

function prova_una_macchia_senza_seme()
    # una sola macchia, che nasce in giorni: nessun organizzatore nel primo giorno
    e, k = corsa_organizzatore(; seme = 1)
    return k[48] < 1.5 && k[end] > 3 && macchie(e) == 1
end

function prova_decide_il_rumore()
    # semi diversi, macchie in posti diversi; lo stesso seme, lo stesso posto
    a, _ = corsa_organizzatore(; seme = 1)
    b, _ = corsa_organizzatore(; seme = 2)
    c = embrione_organizzatore(; seme = 1)
    for _ in 1:3
        passata!(c)
    end
    d = embrione_organizzatore(; seme = 1)
    for _ in 1:3
        passata!(d)
    end
    return argmax(attivatore(a)) != argmax(attivatore(b)) && attivatore(c) == attivatore(d)
end

function prova_senza_diffusivita_diversa_niente_macchia()
    # interruttore: con l'inibitore lento quanto l'attivatore (Dh = Da) non c'è instabilità;
    # la diffusività diversa è la condizione che Müller et al. hanno misurato
    _, k = corsa_organizzatore(; seme = 1, Dh = NODAL_LEFTY.Da)
    return k[end] < 1.2
end

@testset "R5 organizzatore senza seme (Nodal/Lefty)" begin
    @lenta @test prova_una_macchia_senza_seme()
    @lenta @test prova_decide_il_rumore()
    @lenta @test prova_senza_diffusivita_diversa_niente_macchia()
end
