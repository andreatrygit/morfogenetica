# R5: invarianza di scala (specifica, sezione 9; quaderno, "Invarianza di scala";
# prove/attesi.md). La catena a scala vera (test/gastrula.jl) con cellule più piccole: lato 12 µm
# invece di 15, sullo stesso millimetro (2688 cellule invece di 1672), il foglio sempre alto 15
# µm.
#
# Le tre unità dei parametri:
#   lunghezze della forma (fisse)    la scatola, l'altezza del foglio, la griglia della chimica, le
#                                    diffusività e le costanti di Nodal e Lefty
#   lunghezze della cellula          il lato, il volume a riposo, σl (il foglio sta in piedi alla
#   (scalano con il lato)            stessa altezza) e Λ (la stessa tensione apicale per unità di
#                                    lunghezza del tessuto); la secrezione è già per unità di
#                                    volume (attiva_inibisci: ∝ V)
#   tempi (fissi)                    la passata, la meccanica ogni 6 passate, il gene bersaglio
#
# Le prove chiedono la stessa forma in unità della forma: la stessa area impegnata (non lo
# stesso numero di cellule), lo stesso tempo, la stessa profondità. Tolleranze scelte prima
# della corsa: il rumore di espressione è per cellula, e con cellule più piccole la sua media su
# un'area scende di √(15²/12²) = 1.25 volte; la nascita della macchia da un rumore più piccolo
# tarda di ln(1.25)/λ, pochi punti percentuali.
#
# Ogni prova è una funzione che restituisce un Bool.

area_cellula(lato) = 3sqrt(3) / 2 * lato^2

function prova_stessa_area_a_due_taglie()
    e15, _ = corsa_gastrula(15.0)
    e12, _ = corsa_gastrula(12.0)
    I15, I12 = impegnate(e15), impegnate(e12)
    (componenti(e15.tessuto, I15) == 1 && componenti(e12.tessuto, I12) == 1) || return false
    r = length(I12) * area_cellula(12.0) / (length(I15) * area_cellula(15.0))
    return 0.75 <= r <= 1 / 0.75
end

function prova_stesso_tempo_a_due_taglie()
    _, p15 = corsa_gastrula(15.0)
    _, p12 = corsa_gastrula(12.0)
    p12 > 0 || return false
    return 0.85 <= p12 / p15 <= 1 / 0.85
end

function profondità(e)
    I = impegnate(e)
    fuori = setdiff(eachindex(e.cellule), I)
    q(cs) = sum(_quota_apicale_cella(e.tessuto, c) for c in cs) / length(cs)
    return q(I) - q(fuori)
end

function prova_stessa_profondità_a_due_taglie()
    r = profondità(corsa_gastrula(12.0)[1]) / profondità(corsa_gastrula(15.0)[1])
    return 0.7 <= r <= 1 / 0.7
end

@testset "R5 invarianza di scala" begin
    @lenta @test prova_stessa_area_a_due_taglie()
    @lenta @test prova_stesso_tempo_a_due_taglie()
    @lenta @test prova_stessa_profondità_a_due_taglie()
end
