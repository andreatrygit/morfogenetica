# Prestazioni del motore congelato (CLAUDE.md, regola 1, eccezione per le prestazioni).
#
# Una modifica a un file congelato fatta solo per la velocità è ammessa se (a) una prova di tempo
# fallisce prima e passa dopo, e (b) i numeri restano identici bit per bit: qui l'impronta SHA-256
# delle concentrazioni prodotte da un caso fisso, presa prima della modifica.
#
# Il caso: la secrezione di molte cellule in `avanza!` (12_chimica.jl). Con 1672 sorgenti e 937
# sotto-passi (l'inibitore dell'organizzatore, test/organizzatore.jl) durava 3 secondi per
# passata: i pesi della griglia avevano un tipo che il compilatore non riesce a dedurre, e ogni
# operazione sulla griglia passava per una ricerca del metodo a tempo di esecuzione.

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

@testset "prestazioni del motore congelato" begin
    @test prova_secrezione_identica()
    @test_broken prova_secrezione_veloce()
end
