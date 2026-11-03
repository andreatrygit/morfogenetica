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
# D = 18.9 µm²/s = 1134 µm²/min, k = 0.69e-4/s = 0.0041/min. Liberi: κ = 0.35 (saturazione
# dell'autoattivazione: senza, poiché Lefty si elimina più lentamente di Nodal, lo stato uniforme
# oscilla), rumore 20%, c0. Con questi numeri lo stato uniforme è un oscillatore smorzato (periodo
# circa 16 ore, smorzamento in circa 57 ore) e l'instabilità ha lunghezze d'onda fra 0.8 e 1.6 mm
# e cresce di un fattore e ogni 23 ore: una sola macchia nasce in un tessuto di circa 1 mm, in
# giorni. Cellule di 15 µm di lato (26 µm di diametro, come le cellule profonde di una blastula).
#
# La passata è di 5 minuti. Le secrezioni sono decise all'inizio di ogni passata e restano fisse
# fino alla successiva: per il circuito di retroazione è un passo esplicito, stabile solo se la
# passata è più corta di 2δ/ω² (δ lo smorzamento, ω la frequenza dello stato uniforme), 12.9
# minuti con κ = 0.35. Risultato negativo #17: con κ = 0.29 il limite è 1.9 minuti, e con passate
# di 30 minuti lo stato uniforme oscilla per un difetto numerico (prova in fondo al file).
#
# Corretto il 3 novembre: le prove verdi della prima versione usavano κ = 0.29 e passate di 30
# minuti. Le loro conclusioni (una macchia, decide il rumore, niente macchia a diffusività uguali)
# reggono anche così, ma i numeri del contrasto erano gonfiati dall'oscillazione numerica.
#
# Ogni prova è una funzione che restituisce un Bool. Scritte prima del codice (commit rosso, con
# @test_broken), passano dal commit verde.

const NODAL_LEFTY = (Da = 42.0, ka = 0.0073, Dh = 1134.0, kh = 0.0041, κ = 0.35)

function embrione_organizzatore(; seme = 1, rumore = 0.2, Dh = NODAL_LEFTY.Dh, κ = NODAL_LEFTY.κ,
                                dt = 5.0)
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
                        kh = NODAL_LEFTY.kh, κ, rumore, c0)
    return Embrione(t, [f]; parametri = P_PROGRAMMA, dt, campi, seme,
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
        for _ in 1:round(Int, 168 * 60 / e.dt)       # 7 giorni
            passata!(e)
            push!(contrasti, contrasto(e))
        end
        (e, contrasti)
    end
end

function prova_una_macchia_senza_seme()
    # una sola macchia, che nasce in giorni: nessun organizzatore nel primo giorno
    e, k = corsa_organizzatore(; seme = 1)
    return k[round(Int, 24 * 60 / e.dt)] < 1.5 && k[end] > 3 && macchie(e) == 1
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

# Risultato negativo #17: con κ = 0.29 lo stato uniforme è appena smorzato (periodo circa 16
# ore, limite di stabilità della passata 1.9 minuti); con passate di 30 minuti il livello medio di
# Nodal oscilla, per un difetto numerico, del ±30% in 60 ore; con passate di 5 minuti resta
# entro pochi per cento.
function prova_risultato_negativo_17()
    escursione(dt) = begin
        e = embrione_organizzatore(; seme = 1, κ = 0.29, dt)
        medie = Float64[]
        for _ in 1:round(Int, 60 * 60 / dt)
            passata!(e)
            a = attivatore(e)
            push!(medie, sum(a) / length(a))
        end
        maximum(medie) - minimum(medie)
    end
    return escursione(30.0) > 0.4 && escursione(5.0) < 0.15
end

@testset "R5 organizzatore senza seme (Nodal/Lefty)" begin
    @lenta @test prova_una_macchia_senza_seme()
    @lenta @test prova_decide_il_rumore()
    @lenta @test prova_senza_diffusivita_diversa_niente_macchia()
    @lenta @test prova_risultato_negativo_17()
end
