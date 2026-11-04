# Prestazioni del motore congelato (CLAUDE.md, regola 1, eccezione per le prestazioni).
#
# Una modifica a un file congelato fatta solo per la velocità è ammessa se (a) una prova di tempo
# fallisce prima e passa dopo, e (b) i numeri restano identici bit per bit: qui l'impronta SHA-256
# delle concentrazioni prodotte da un caso fisso, presa prima della modifica.
#
# Il caso: la secrezione di molte cellule in `avanza!` (12_chimica.jl). Con 1672 sorgenti e 937
# sotto-passi (l'inibitore dell'organizzatore, test/organizzatore.jl) durava 3 secondi per
# passata: i pesi della griglia avevano un tipo che il compilatore non riesce a dedurre, e ogni
# operazione sulla griglia passava per una ricerca del metodo a tempo di esecuzione. Corretto il
# 3 novembre scrivendo il tipo per esteso: da 2.4 a 0.06 s, stessa impronta.

using SHA

function caso_secrezione()
    campo = Campo(Sostanza(:lefty, 1134.0, 0.0041); origine = (0.0, 0.0, 0.0),
                  h = (12.0, 12.0, 15.0), n = (82, 82, 1),
                  bordi = (bordo_periodico, bordo_periodico, bordo_chiuso))
    campo.c .= 1.0
    rng = MersenneTwister(11)
    sorgenti = [Scambio((984.0 * rand(rng), 984.0 * rand(rng), 7.5), 1.0 + rand(rng), 0.0)
                for _ in 1:1672]
    return campo, sorgenti
end

impronta(campo) = bytes2hex(sha256(reinterpret(UInt8, vec(campo.c))))

# presa il 3 novembre 2026, prima della modifica per le prestazioni
const IMPRONTA_SECREZIONE = "cc851011bc73907462dd266c800faf2990c88cd60712b06e0c31e74d038bed9c"

function prova_secrezione_identica()
    campo, sorgenti = caso_secrezione()
    avanza!(campo, 30.0, sorgenti)
    return impronta(campo) == IMPRONTA_SECREZIONE
end

function prova_secrezione_veloce()
    campo, sorgenti = caso_secrezione()
    avanza!(campo, 1.0, sorgenti)                    # compilazione
    campo, sorgenti = caso_secrezione()
    t = @elapsed avanza!(campo, 30.0, sorgenti)
    return t < 0.5
end

# Il gradiente e l'energia del guscio (07_guscio.jl), 4 novembre. Il gradiente creava un vettore
# nuovo per ogni faccia a ogni chiamata (e uno di immagini per ogni faccia di ogni cellula): sul
# foglio di 1672 cellule 7.6 MiB e circa 7 ms per gradiente, e ogni matrice del solutore sparso
# allocava 4 GiB, con il garbage collector che ferma tutti i fili. Le impronte sono prese prima
# della modifica, su un foglio a bordo bloccato, uno periodico e uno con uno scutoide (facce di
# lunghezze diverse).
function tessuti_per_impronta()
    rng = MersenneTwister(21)
    a = foglio_esagonale([Intenzioni3D(1.0, 1.0, 1.0, c <= 3 ? 2.0 : 0.0, 1000.0) for c in 1:19], 2;
                         lato = 5.0, altezza = 15.0, bordo_fisso = true)
    b = foglio_periodico(fill(Intenzioni3D(1.0, 0.8, 1.2, 0.5, 900.0), 16), 4, 4; lato = 5.0, altezza = 15.0)
    c = foglio_esagonale(fill(Intenzioni3D(1.0, 1.0, 1.0, 1.0, 1000.0), 37), 3; lato = 5.0, altezza = 15.0,
                         bordo_fisso = true)
    ap(t, k) = t.facce[Morfogenetica._faccia_apicale_indice(t, k)].vertici
    spigolo_a_triangolo!(c, Tuple(intersect(ap(c, 1), ap(c, 2)))...; lunghezza_nuova = 1.0)
    for t in (a, b, c)
        n = 3 * Morfogenetica._nvert(t)
        t.x[1:n] .+= 0.2 .* (rand(rng, n) .- 0.5)
    end
    return (a, b, c)
end

function impronta_guscio()
    P = ParametriMeccanici(K = 0.1, η = 1.0, lunghezza_minima = 0.05, k_contatto = 0.0, d_contatto = 0.5)
    dati = UInt8[]
    for t in tessuti_per_impronta()
        g = gradiente!(similar(t.x), t, P)
        append!(dati, reinterpret(UInt8, g))
        append!(dati, reinterpret(UInt8, [energia(t, P)]))
        append!(dati, reinterpret(UInt8, [volume_cellula(t, c) for c in eachindex(t.cellule)]))
    end
    return bytes2hex(sha256(dati))
end

# presa il 4 novembre 2026, prima della modifica
const IMPRONTA_GUSCIO = "f934c86719faeea8727bc11e427913512efb1d38ad6e72f275e9fa15157755b1"

prova_guscio_identico() = impronta_guscio() == IMPRONTA_GUSCIO

function prova_gradiente_senza_allocazioni()
    P = ParametriMeccanici(K = 0.1, η = 1.0, lunghezza_minima = 0.05, k_contatto = 0.0, d_contatto = 0.5)
    t = foglio_periodico(fill(Intenzioni3D(1.0, 1.0, 2sqrt(3), 0.0, 8768.0), 38 * 44), 38, 44;
                         lato = 15.0, altezza = 15.0)
    g = similar(t.x)
    gradiente!(g, t, P)                              # compilazione
    byte = @allocated gradiente!(g, t, P)
    tempo = minimum(@elapsed(gradiente!(g, t, P)) for _ in 1:5)
    return byte < 512 * 1024 && tempo < 0.004
end

@testset "prestazioni del motore congelato" begin
    @test prova_secrezione_identica()
    @test prova_secrezione_veloce()
    @test prova_guscio_identico()
    @test_broken prova_gradiente_senza_allocazioni()
end
