# R1 (docs/roadmap.md): le due operazioni elementari dei modelli 3D a vertici, in forma generale.
#
#   spigolo_a_triangolo!(t, i, j; lunghezza_nuova)   (I→H)
#   triangolo_a_spigolo!(t, f; lunghezza_nuova)      (H→I, l'inversa)
#   problemi_topologici(t; monostrato = true)        → elenco delle violazioni (vuoto se va bene)
#   rilassa_con_operazioni!(t, p; soglia, area_soglia, lunghezza_nuova, tolleranza)
#
# Il mezzo esterno conta come una "cellula" per ogni etichetta di faccia libera (apicale, basale,
# laterale): così ogni spigolo ha tre cellule attorno e ogni vertice quattro.
#
# Ogni prova è una funzione che restituisce un Bool. Scritte prima del codice (commit rosso, con
# @test_broken), passano dal commit verde.

const P_OPERAZIONI = P_SCUTOIDI

foglio_operazioni() = foglio_scutoidi()

# Firma della topologia, indipendente da come sono numerati i vertici: per ogni cellula, chi sta
# dall'altra parte di ogni sua faccia (un'altra cellula o il mezzo, per etichetta) e quanti
# vertici ha la faccia.
function firma(t)
    return [sort([begin
                      f = t.facce[fi]
                      altra = length(f.cellule) == 2 ? only(setdiff(f.cellule, [c])) : Symbol(f.etichetta)
                      (string(altra), length(f.vertici))
                  end for fi in t.cellule[c].facce]) for c in eachindex(t.cellule)]
end

# chi tocca chi, senza contare i vertici delle facce
vicinato(t) = [sort([string(length(t.facce[fi].cellule) == 2 ?
                            only(setdiff(t.facce[fi].cellule, [c])) : Symbol(t.facce[fi].etichetta))
                     for fi in t.cellule[c].facce]) for c in eachindex(t.cellule)]

function prova_invarianti_dei_fogli()
    t = foglio_operazioni()
    isempty(problemi_topologici(t)) || return false
    q = foglio_periodico(fill(PRISMA, 16), 4, 4; lato = 5.0, altezza = 15.0)
    isempty(problemi_topologici(q)) || return false
    s = foglio_operazioni()
    scambio_apicale!(s, spigolo_apicale_12(s)...; lunghezza_nuova = 1.0)
    return isempty(problemi_topologici(s))
end

function prova_invarianti_violati_si_vedono()
    # una faccia tolta da una cellula: la sua superficie non è più chiusa
    t = foglio_operazioni()
    cel = t.cellule[1]
    t.cellule[1] = Cellula3D(cel.facce[2:end], cel.segni[2:end], cel.traslazioni[2:end])
    return !isempty(problemi_topologici(t))
end

function prova_spigolo_apicale_come_scambio_apicale()
    a = foglio_operazioni(); b = foglio_operazioni()
    i, j = spigolo_apicale_12(a)
    r = spigolo_a_triangolo!(a, i, j; lunghezza_nuova = 1.0)
    scambio_apicale!(b, i, j; lunghezza_nuova = 1.0)
    return r !== nothing && firma(a) == firma(b) && isempty(problemi_topologici(a))
end

function prova_triangolo_a_spigolo_inverte()
    t = foglio_operazioni()
    prima = firma(t)
    i, j = spigolo_apicale_12(t)
    r = spigolo_a_triangolo!(t, i, j; lunghezza_nuova = 1.0)
    r === nothing && return false
    firma(t) != prima || return false
    s = triangolo_a_spigolo!(t, r.triangolo; lunghezza_nuova = 1.0)
    return s !== nothing && firma(t) == prima && isempty(problemi_topologici(t))
end

function prova_continuita()
    salti = Float64[]
    for ε in (1e-2, 1e-3)
        t = accorcia_apicale!(foglio_operazioni(), spigolo_apicale_12(foglio_operazioni())..., ε)
        E = energia(t, P_OPERAZIONI)
        V = sum(volume_cellula(t, c) for c in eachindex(t.cellule))
        r = spigolo_a_triangolo!(t, spigolo_apicale_12(foglio_operazioni())...; lunghezza_nuova = ε)
        r === nothing && return false
        isapprox(sum(volume_cellula(t, c) for c in eachindex(t.cellule)), V; rtol = 1e-12) || return false
        E1 = energia(t, P_OPERAZIONI)
        # e indietro: il triangolo appena nato è lungo ε, torna a spigolo con un salto ∝ ε
        s = triangolo_a_spigolo!(t, r.triangolo; lunghezza_nuova = ε)
        s === nothing && return false
        push!(salti, abs(E1 - E) + abs(energia(t, P_OPERAZIONI) - E1))
    end
    return 9 < salti[1] / salti[2] < 11
end

function prova_forze_esatte_dopo_le_operazioni()
    t = foglio_operazioni()
    r = spigolo_a_triangolo!(t, spigolo_apicale_12(t)...; lunghezza_nuova = 1.0)
    r === nothing && return false
    errore_gradiente(t, P_OPERAZIONI) < 1e-7 || return false
    triangolo_a_spigolo!(t, r.triangolo; lunghezza_nuova = 1.0) === nothing && return false
    return errore_gradiente(t, P_OPERAZIONI) < 1e-7
end

function prova_scambio_intero_in_due_passi()
    # I→H sullo spigolo apicale, poi H→I sulla parete fra A e B, diventata un triangolo: deve
    # dare chi-tocca-chi dello scambio fra prismi di 08_scambi.jl. Sul pezzo di foglio a bordo
    # bloccato: le operazioni generali non lavorano ancora sui fogli periodici.
    prismi = foglio_operazioni()
    i, j = spigolo_apicale_12(prismi)                  # vertici apicali 2p - 1 e 2q - 1
    scambia!(prismi, (i + 1) ÷ 2, (j + 1) ÷ 2; lunghezza_nuova = 1.0) === nothing && return false
    t = foglio_operazioni()
    r = spigolo_a_triangolo!(t, i, j; lunghezza_nuova = 1.0)
    r === nothing && return false
    A, B = r.perdono
    fAB = findfirst(f -> f.etichetta == laterale && A in f.cellule && B in f.cellule, t.facce)
    fAB === nothing && return false
    length(t.facce[fAB].vertici) == 3 || return false
    s = triangolo_a_spigolo!(t, fAB; lunghezza_nuova = 1.0)
    return s !== nothing && vicinato(t) == vicinato(prismi) && isempty(problemi_topologici(t))
end

function prova_scutoide_forzato_torna_indietro()
    regolare = foglio_operazioni()
    minimizza!(regolare, P_OPERAZIONI; tolleranza = 1e-9, iterazioni_massime = 300)
    t = foglio_operazioni()
    spigolo_a_triangolo!(t, spigolo_apicale_12(t)...; lunghezza_nuova = 1.0) === nothing && return false
    r = rilassa_con_operazioni!(t, P_OPERAZIONI; soglia = 0.3, area_soglia = 0.3,
                                lunghezza_nuova = 1.0, tolleranza = 1e-9)
    return r.convergito && vicinato(t) == vicinato(regolare) &&
           isapprox(energia(t, P_OPERAZIONI), energia(regolare, P_OPERAZIONI); rtol = 1e-9)
end

function prova_bordo_libero_trova_equilibrio()
    # risultato negativo #13: con le operazioni generali deve esistere un equilibrio
    t = foglio_esagonale(fill(PRISMA, 19), 2; lato = 5.0, altezza = 15.0)
    r = rilassa_con_operazioni!(t, P_OPERAZIONI; soglia = 0.3, area_soglia = 0.3,
                                lunghezza_nuova = 1.0, tolleranza = 1e-9)
    return r.convergito && isempty(problemi_topologici(t))
end

function prova_cellula_isolata_si_arrotonda()
    t = foglio_esagonale([Intenzioni3D(1.0, 1.0, 0.9, 0.0, 1000.0)], 0; lato = 5.0, altezza = 15.0)
    r = rilassa_con_operazioni!(t, P_3D; soglia = 0.3, area_soglia = 0.3,
                                lunghezza_nuova = 1.0, tolleranza = 1e-9)
    return r.convergito && isempty(problemi_topologici(t)) && energia(t, P_3D) < 456.03
end

@testset "R1 operazioni elementari" begin
    @test prova_invarianti_dei_fogli()
    @test prova_invarianti_violati_si_vedono()
    @test prova_spigolo_apicale_come_scambio_apicale()
    @test prova_triangolo_a_spigolo_inverte()
    @test prova_continuita()
    @test prova_forze_esatte_dopo_le_operazioni()
    @test prova_scambio_intero_in_due_passi()
    @test prova_scutoide_forzato_torna_indietro()
    @test prova_bordo_libero_trova_equilibrio()
    @test prova_cellula_isolata_si_arrotonda()
end
