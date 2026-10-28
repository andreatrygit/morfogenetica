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

# Guscio 3D: prisma esagonale regolare a volume V fissato (docs/piano.md, sezione 8).
# c = 3√3/2 è il rapporto fra l'area dell'esagono e il quadrato del lato.
const ESAGONO = 3sqrt(3) / 2

"Rapporto altezza/lato in un foglio di prismi, facce laterali condivise: √3·(σa+σb)/σl."
prisma_foglio_rapporto(σa, σb, σl) = sqrt(3) * (σa + σb) / σl

"Prisma isolato, sei facce laterali libere: lato s con s³ = 3σl·V/((σa+σb)·c²). Vale solo con σa = σb."
prisma_isolato_lato(σa, σb, σl, V) = cbrt(3σl * V / ((σa + σb) * ESAGONO^2))
