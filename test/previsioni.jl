# Le previsioni esatte del modello a vertici laterale (docs/piano.md, sezione 1).
# Sono l'oracolo delle prove: stanno qui e non nel motore.

"Semiangolo del cuneo all'equilibrio per una striscia uniforme: sin α = Δ/Γl."
cuneo_semiangolo(Γa, Γb, Γl) = asin((Γa - Γb) / Γl)

"Altezza della cellula con semiangolo α imposto e area A."
cuneo_altezza(Γa, Γb, Γl, A, α) = sqrt((Γa + Γb) * A / (cos(α) * (Γl - (Γa - Γb) * sin(α))))

"Rapporto altezza/larghezza nel foglio piatto (Γa = Γb)."
piatto_rapporto(Γa, Γb, Γl) = (Γa + Γb) / Γl

"Condizione di validità: il lato apicale resta positivo."
cuneo_valido(Γa, Γb, Γl) = Γl^2 - (Γa - Γb)^2 > (Γa - Γb) * (Γa + Γb)
