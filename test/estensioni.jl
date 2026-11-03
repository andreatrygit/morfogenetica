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

@testset "estensioni compatibili" begin
    @test prova_passata_predefinita_identica()
end
