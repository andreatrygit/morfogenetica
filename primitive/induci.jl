# induci: un gene bersaglio che si accumula con il segnale. È il modo plausibile di leggere un
# morfogeno: la cellula non decide da una soglia sulla concentrazione istantanea, ma dal prodotto
# di un gene che il segnale induce, che si accumula mentre il segnale c'è e decade quando manca.
# Così filtra il rumore e i picchi brevi, e la risposta dipende dalla durata dell'esposizione.
#
# Regola, sui canali "dentro di sé" (il prodotto del gene, stato[voce]) e "misura assoluta" (la
# concentrazione s alla propria posizione):
#
#   dg/dt = β·H(s) − γ·g,     H(s) = sⁿ / (Kⁿ + sⁿ)
#
# integrata in modo esatto nella passata, con s costante (letto all'inizio della passata):
# g ← g·e^(−γ dt) + (β H(s)/γ)·(1 − e^(−γ dt)). Esatta, quindi stabile per ogni passata: il
# problema del passo esplicito (risultato negativo #17) qui non c'è. Predicato di fine: nessuno;
# una fase successiva può leggere stato[voce] (per esempio `costringi(; voce)`).
#
# Meccanismo: Nodal attiva Smad2, che induce la trascrizione dei geni bersaglio; la risposta di
#   ciascun gene dipende dalla velocità e dal momento della sua induzione, non da una soglia sulla
#   concentrazione del ligando (geni a lungo e a corto raggio nel pesce zebra).
# Fonte: Dubrulle et al., "Response to Nodal morphogen gradient is determined by the kinetics of
#   target gene induction", eLife 4, e05042 (2015), doi:10.7554/eLife.05042.

export induci

function _regola_induci(p::Percezione, par)
    s = max(get(p.concentrazioni, par.sostanza, 0.0), 0.0)
    H = s^par.n / (par.K^par.n + s^par.n)
    g = get(p.stato, par.voce, 0.0)
    f = exp(-par.γ * par.dt)
    return Uscita(; stato = Dict{Symbol,Float64}(par.voce => g * f + par.β * H / par.γ * (1 - f)))
end

_fine_induci(p::Percezione, par) = false

"""
    induci(; sostanza, voce, β, γ, K, n, dt = 5.0)

La fase in cui il prodotto di un gene bersaglio, `stato[voce]`, si accumula con la
concentrazione di `sostanza`: velocità massima `β` (unità/min), decadimento `γ` (1/min),
soglia `K` e cooperatività `n` della curva di Hill. `dt` (min) deve essere la passata
dell'embrione: la regola non conosce l'orologio della simulazione, solo il proprio passo.
"""
induci(; sostanza, voce, β, γ, K, n, dt = 5.0) =
    Fase(:induci, _regola_induci,
         (sostanza = sostanza, voce = voce, β = Float64(β), γ = Float64(γ), K = Float64(K),
          n = Float64(n), dt = Float64(dt)),
         _fine_induci; scrive = (), stato = (voce,))
