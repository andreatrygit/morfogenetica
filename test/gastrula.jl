# R5: la catena a scala vera, organizzatore poi invaginazione (programmi/gastrula.jl).
#
# Foglio periodico di 1672 cellule (circa 1 × 1 mm, cellule di 26 µm), Nodal e Lefty come in
# test/organizzatore.jl (Müller et al. 2012, κ = 0.35), il gene bersaglio con tempo di
# dimezzamento di 60 minuti (β = γ, K = 1.4 c0, n = 4: liberi), impegno a metà del massimo,
# Λ = 8 nN. Passata di 5 minuti (la chimica la chiede, risultato negativo #17), meccanica ogni 6
# passate (30 minuti; specifica, decisione 6), solutore sparso con la matrice esatta (16).
#
# Le prove sono comportamenti: nei primi quattro giorni nessuna cellula si impegna (l'organizzatore
# nasce in giorni); poi una sola regione si impegna, e si invagina dal lato basale.

# Con `lato` diverso da 15 µm (test/scala.jl) il foglio copre lo stesso millimetro con cellule
# più piccole, alto sempre 15 µm; le grandezze in unità di cellula scalano con il lato: il volume
# a riposo, σl (perché il foglio stia in piedi alla stessa altezza: altezza/lato = √3·(σa + σb)/σl)
# e Λ (la stessa tensione apicale per unità di lunghezza del tessuto, Λ/lato).
const _GASTRULA = Dict{Float64,Any}()
function foglio_a_taglia(lato)
    lato == 15.0 && return foglio_grande_rilassato()
    σl = 2sqrt(3) * lato / 15
    nx = round(Int, 38 * 15 / lato)
    ny = 2round(Int, 22 * 15 / lato)
    V0 = 3sqrt(3) / 2 * lato^2 * 15
    t = foglio_periodico(fill(Intenzioni3D(1.0, 1.0, σl, 0.0, V0), nx * ny), nx, ny; lato, altezza = 15.0)
    minimizza_sparso!(t, P_PROGRAMMA; tolleranza = 1e-5, iterazioni_massime = 100, hessiana = :esatta)
    return t
end

function corsa_gastrula(lato = 15.0)
    haskey(_GASTRULA, lato) && return _GASTRULA[lato]
    t = foglio_a_taglia(lato)
    Lx, Ly = t.x[end-1], t.x[end]
    griglia(nome, D, k) = Campo(Sostanza(nome, D, k); origine = (0.0, 0.0, 0.0),
                                h = (Lx / 82, Ly / 82, 15.0), n = (82, 82, 1),
                                bordi = (bordo_periodico, bordo_periodico, bordo_chiuso))
    campi = Dict(:nodal => griglia(:nodal, NODAL_LEFTY.Da, NODAL_LEFTY.ka),
                 :lefty => griglia(:lefty, NODAL_LEFTY.Dh, NODAL_LEFTY.kh))
    for c in values(campi)
        c.c .= 1.0
    end
    γ = log(2) / 60
    programma = gastrula(; ka = NODAL_LEFTY.ka, kh = NODAL_LEFTY.kh, κ = NODAL_LEFTY.κ, β = γ, γ,
                         K = 1.4, n = 4, impegno = 0.5, Λ = 8.0 * lato / 15, dt = 5.0)
    e = Embrione(t, programma; parametri = P_PROGRAMMA, dt = 5.0, campi, seme = 1,
                 solutore = :esatto, tolleranza = 1e-5, meccanica_ogni = 6,
                 interruttori = Interruttori(; topologia = false, popolazione = false))
    primo_impegno = 0
    for k in 1:round(Int, 7 * 24 * 60 / 5)              # 7 giorni
        passata!(e)
        primo_impegno == 0 && any(c -> c.fase == 2, e.cellule) && (primo_impegno = k)
    end
    _GASTRULA[lato] = (e, primo_impegno)
    return _GASTRULA[lato]
end

impegnate(e) = [c for c in eachindex(e.cellule) if get(e.cellule[c].stato, :gene, 0.0) > 0.5]

function componenti(t, insieme_)
    resto = Set(insieme_); n = 0
    while !isempty(resto)
        n += 1; coda = [pop!(resto)]
        while !isempty(coda)
            c = pop!(coda)
            for v in vicine_di(t, c)
                v in resto && (delete!(resto, v); push!(coda, v))
            end
        end
    end
    return n
end

function _quota_apicale_cella(t, c)
    cel = t.cellule[c]
    k = findfirst(fi -> t.facce[fi].etichetta == apicale, cel.facce)
    return Morfogenetica._baricentro(Morfogenetica._vertici(t, t.x, t.facce[cel.facce[k]],
                                                            cel.traslazioni[k]))[3]
end

function prova_l_organizzatore_nasce_in_giorni()
    e, primo = corsa_gastrula()
    return primo * 5 > 4 * 24 * 60
end

function prova_una_sola_regione_si_impegna()
    e, _ = corsa_gastrula()
    I = impegnate(e)
    return 20 <= length(I) <= 600 && componenti(e.tessuto, I) == 1
end

function prova_la_regione_si_invagina()
    e, _ = corsa_gastrula()
    I = impegnate(e)
    fuori = setdiff(eachindex(e.cellule), I)
    q(cs) = sum(_quota_apicale_cella(e.tessuto, c) for c in cs) / length(cs)
    return all(e.cellule[c].intenzioni.Λa == 8.0 for c in I) && q(I) - q(fuori) < -3.0
end

@testset "R5 la catena a scala vera" begin
    @lenta @test prova_l_organizzatore_nasce_in_giorni()
    @lenta @test prova_una_sola_regione_si_impegna()
    @lenta @test prova_la_regione_si_invagina()
end
