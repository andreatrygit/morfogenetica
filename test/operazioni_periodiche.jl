# R1, seconda parte (docs/roadmap.md): le operazioni elementari di motore/11_operazioni.jl sui
# fogli periodici. Un vertice compare nelle facce con un'immagine (quanti periodi aggiungere):
# uno spigolo che attraversa il bordo della scatola ha i due estremi con immagini diverse.
#
# Ogni prova è una funzione che restituisce un Bool. Commit rosso: marcate @test_broken finché
# il codice non c'è (oggi le operazioni rispondono `nothing` sui fogli periodici).

foglio_periodico_op() = foglio_scambi()            # 4 × 4 prismi, periodico

# Uno spigolo apicale fra due cellule: con `di_bordo` uno che attraversa la scatola.
function spigolo_apicale_periodico(t; di_bordo)
    for f in t.facce
        f.etichetta == apicale || continue
        n = length(f.vertici)
        for k in 1:n
            k2 = mod1(k + 1, n)
            (f.immagini[k] != f.immagini[k2]) == di_bordo && return (f.vertici[k], f.vertici[k2])
        end
    end
    error("nessuno spigolo apicale con di_bordo = $di_bordo")
end

# Avvicina i due vertici apicali a distanza ε sul loro asse, tenendo conto delle immagini.
function accorcia_apicale_periodico!(t, a1, a2, ε)
    f = first(f for f in t.facce if f.etichetta == apicale && a1 in f.vertici && a2 in f.vertici)
    ia = f.immagini[findfirst(==(a1), f.vertici)]
    ib = f.immagini[findfirst(==(a2), f.vertici)]
    M = Morfogenetica
    Lx, Ly = M._periodo(t, t.x)
    pa, pb = M._pos(t, t.x, a1, ia), M._pos(t, t.x, a2, ib)
    m = (pa .+ pb) ./ 2
    e = M._meno(pb, pa) ./ M._norma(M._meno(pb, pa))
    na = m .- e .* (ε / 2); nb = m .+ e .* (ε / 2)
    t.x[3a1-2:3a1] .= (na[1] - ia[1] * Lx, na[2] - ia[2] * Ly, na[3])
    t.x[3a2-2:3a2] .= (nb[1] - ib[1] * Lx, nb[2] - ib[2] * Ly, nb[3])
    return t
end

function prova_periodico_andata_e_ritorno()
    for di_bordo in (false, true)
        t = foglio_periodico_op()
        prima = firma(t)
        r = spigolo_a_triangolo!(t, spigolo_apicale_periodico(t; di_bordo)...; lunghezza_nuova = 1.0)
        r === nothing && return false
        (isempty(problemi_topologici(t)) && firma(t) != prima) || return false
        # in un foglio periodico ogni spigolo apicale sta fra due facce apicali
        all(all(==(2), f.condivisione) for f in t.facce if f.etichetta == apicale) || return false
        s = triangolo_a_spigolo!(t, r.triangolo; lunghezza_nuova = 1.0)
        (s !== nothing && firma(t) == prima && isempty(problemi_topologici(t))) || return false
    end
    return true
end

function prova_periodico_forze_esatte()
    t = foglio_periodico_op()
    r = spigolo_a_triangolo!(t, spigolo_apicale_periodico(t; di_bordo = true)...; lunghezza_nuova = 1.0)
    r === nothing && return false
    errore_gradiente(t, P_OPERAZIONI) < 1e-7 || return false
    triangolo_a_spigolo!(t, r.triangolo; lunghezza_nuova = 1.0) === nothing && return false
    return errore_gradiente(t, P_OPERAZIONI) < 1e-7
end

function prova_periodico_continuita()
    salti = Float64[]
    for ε in (1e-2, 1e-3)
        t = foglio_periodico_op()
        i, j = spigolo_apicale_periodico(t; di_bordo = true)
        accorcia_apicale_periodico!(t, i, j, ε)
        E = energia(t, P_OPERAZIONI)
        V = sum(volume_cellula(t, c) for c in eachindex(t.cellule))
        r = spigolo_a_triangolo!(t, i, j; lunghezza_nuova = ε)
        r === nothing && return false
        isapprox(sum(volume_cellula(t, c) for c in eachindex(t.cellule)), V; rtol = 1e-12) || return false
        E1 = energia(t, P_OPERAZIONI)
        s = triangolo_a_spigolo!(t, r.triangolo; lunghezza_nuova = ε)
        s === nothing && return false
        push!(salti, abs(E1 - E) + abs(energia(t, P_OPERAZIONI) - E1))
    end
    return 9 < salti[1] / salti[2] < 11
end

function prova_periodico_scambio_intero_in_due_passi()
    # come in test/operazioni.jl, ma nel foglio periodico, con uno spigolo interno alla scatola
    # e uno che ne attraversa il bordo
    for di_bordo in (false, true)
        prismi = foglio_periodico_op()
        i, j = spigolo_apicale_periodico(prismi; di_bordo)
        scambia!(prismi, (i + 1) ÷ 2, (j + 1) ÷ 2; lunghezza_nuova = 1.0) === nothing && return false
        t = foglio_periodico_op()
        r = spigolo_a_triangolo!(t, i, j; lunghezza_nuova = 1.0)
        r === nothing && return false
        A, B = r.perdono
        fAB = findfirst(f -> f.etichetta == laterale && A in f.cellule && B in f.cellule, t.facce)
        fAB === nothing && return false
        triangolo_a_spigolo!(t, fAB; lunghezza_nuova = 1.0) === nothing && return false
        (vicinato(t) == vicinato(prismi) && isempty(problemi_topologici(t))) || return false
    end
    return true
end

function prova_periodico_scutoide_torna_indietro()
    regolare = foglio_periodico_op()
    minimizza!(regolare, P_OPERAZIONI; tolleranza = 1e-9, iterazioni_massime = 300)
    t = foglio_periodico_op()
    spigolo_a_triangolo!(t, spigolo_apicale_periodico(t; di_bordo = true)...;
                         lunghezza_nuova = 1.0) === nothing && return false
    r = rilassa_con_operazioni!(t, P_OPERAZIONI; soglia = 0.3, area_soglia = 0.3,
                                lunghezza_nuova = 1.0, tolleranza = 1e-9)
    return r.convergito && vicinato(t) == vicinato(regolare) && isempty(problemi_topologici(t)) &&
           isapprox(energia(t, P_OPERAZIONI), energia(regolare, P_OPERAZIONI); rtol = 1e-9)
end

@testset "R1 operazioni nei fogli periodici" begin
    @test_broken prova_periodico_andata_e_ritorno()
    @test_broken prova_periodico_forze_esatte()
    @test_broken prova_periodico_continuita()
    @test_broken prova_periodico_scambio_intero_in_due_passi()
    @test_broken prova_periodico_scutoide_torna_indietro()
end
