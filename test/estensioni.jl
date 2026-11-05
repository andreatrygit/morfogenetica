# Estensioni compatibili dei file congelati (CLAUDE.md, regola 1, seconda eccezione).
#
# Un file congelato può guadagnare un'opzione nuova se, con l'opzione al valore predefinito, i
# risultati restano identici bit per bit: qui l'impronta SHA-256 dei risultati di un caso fisso,
# presa prima della modifica. Le estensioni di oggi: la scelta del solutore nella passata
# (`Embrione(...; solutore = :pieno o :sparso)`, 13_programma.jl).

using SHA

# Tre passate di costrizione con chimica e meccanica sul pezzo di 19 cellule a bordo bloccato.
function caso_passata()
    t = foglio_esagonale(fill(PRISMA, 19), 2; lato = 5.0, altezza = 15.0, bordo_fisso = true)
    campo = Campo(Sostanza(:segnale, 6.0, 0.1); origine = (-25.0, -25.0, -5.0), h = 2.5,
                  n = (20, 20, 10))
    stati = [Dict{Symbol,Float64}(:sorgente => (c == 1 ? 1.0 : 0.0)) for c in eachindex(t.cellule)]
    programma = invaginazione(; Λ = 3.0, soglia = 0.001, tasso = 20.0)
    return Embrione(t, programma; parametri = P_PROGRAMMA, dt = 1.0, campi = Dict(:segnale => campo),
                    stati, interruttori = Interruttori(; topologia = false))
end

impronta_passata(e) = bytes2hex(sha256(vcat(reinterpret(UInt8, e.tessuto.x),
                                            reinterpret(UInt8, vec(e.campi[:segnale].c)))))

# presa il 3 novembre 2026, prima dell'estensione
const IMPRONTA_PASSATA = "a78fc939e954775a0cbd977a7ae8c18257963a70c74428807693a1ffc344f340"

function prova_passata_predefinita_identica()
    e = caso_passata()
    for _ in 1:3
        passata!(e)
    end
    return impronta_passata(e) == IMPRONTA_PASSATA
end

# meccanica_ogni (4 novembre; specifica, decisione 6: scale di tempo separate): la meccanica
# avanza solo ogni m passate, con un passo di m·dt. Con un programma che non cambia le
# intenzioni, m = 3 con passata di 1 minuto dà esattamente lo stesso tessuto di m = 1 con
# passata di 3 minuti; fra un passo meccanico e l'altro il tessuto non si muove.
function prova_meccanica_ogni()
    programma() = [Fase(:f, (p, par) -> Uscita(; intenzioni = (Λa = 3.0,)), NamedTuple(),
                        (p, par) -> false; scrive = (:Λa,))]
    tessuto() = foglio_esagonale(fill(PRISMA, 19), 2; lato = 5.0, altezza = 15.0, bordo_fisso = true)
    a = Embrione(tessuto(), programma(); parametri = P_PROGRAMMA, dt = 1.0, meccanica_ogni = 3,
                 interruttori = Interruttori(; topologia = false))
    b = Embrione(tessuto(), programma(); parametri = P_PROGRAMMA, dt = 3.0,
                 interruttori = Interruttori(; topologia = false))
    x0 = copy(a.tessuto.x)
    passata!(a); passata!(a)
    a.tessuto.x == x0 || return false
    passata!(a); passata!(b)
    a.tessuto.x == b.tessuto.x || return false
    passata!(a); passata!(a); passata!(a); passata!(b)
    return a.tessuto.x == b.tessuto.x
end

# La matrice delle derivate seconde esatta (5 novembre): `minimizza_sparso!(...; hessiana =
# :differenze o :esatta)` in 15_sparso.jl e `solutore = :esatto` dell'Embrione in 13_programma.jl.
# Con i valori predefiniti, lo stesso minimo bit per bit: l'impronta di minimizza_sparso! su un
# foglio periodico disturbato, presa prima.
function impronta_sparso()
    P = ParametriMeccanici(K = 0.1, η = 1.0, lunghezza_minima = 0.05, k_contatto = 0.0, d_contatto = 0.5)
    t = foglio_periodico(fill(Intenzioni3D(1.0, 1.0, 1.0, 0.0, 1000.0), 16), 4, 4; lato = 5.0, altezza = 15.0)
    rng = MersenneTwister(31)
    n = length(t.x) - 2
    t.x[1:n] .+= 0.3 .* (rand(rng, n) .- 0.5)
    for c in 1:3
        t.intenzioni[c] = Intenzioni3D(1.0, 1.0, 1.0, 2.0, 900.0)
    end
    minimizza_sparso!(t, P; tolleranza = 1e-9, iterazioni_massime = 300, ancora = copy(t.x), peso = 0.5)
    return bytes2hex(sha256(reinterpret(UInt8, t.x)))
end

const IMPRONTA_SPARSO = "e79d32d0018a2c2569ef466c618b105994cf7d625f505c150902364a11575ef2"

prova_sparso_predefinito_identico() = impronta_sparso() == IMPRONTA_SPARSO

@testset "estensioni compatibili" begin
    @test prova_passata_predefinita_identica()
    @test prova_meccanica_ogni()
    @test prova_sparso_predefinito_identico()
end
