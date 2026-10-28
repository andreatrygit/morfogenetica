# Numeri attesi (motore di riferimento)

Misure del motore JavaScript congelato in `riferimento/`, al 27 ottobre 2026. Servono come
**riferimento**, non come bersaglio: nel motore nuovo le unità cambiano (micrometri e minuti al
posto di "unità di forma" e "tick"), quindi i valori assoluti non si trasferiscono. Ciò che si
trasferisce sono i **rapporti** e i **comportamenti qualitativi**: che la forma si chiuda, che si
fermi, che ripari una ferita, e che tutto questo non dipenda dal raggio cellulare.

Come si eseguono (dalla cartella `riferimento/`):

```
node -e "globalThis.window=globalThis; for(const f of ['macchine.js','motore3d.js','suite.js'])
  eval(require('fs').readFileSync(f,'utf8'));"
```

Gli script sono `suite.js` (le quattro forme con ferita), `ablaz.js` (ablazione dei meccanismi),
`p2_gastrula.js` e `p3_gastrula.js` (gastrulazione 2D e 3D), `p_simmetria.js` (la legge di
simmetria), `diag2.js` (profili radiali), `vedi.py` (vista ASCII di una sezione).

## Le quattro forme di base

Un solo motore per 2D e 3D (interruttore `piano`). Ferita = asportato un quarto del corpo.
Tre ripetizioni ciascuna, 12 su 12 riuscite.

| forma | n | bordo esterno | bordo interno | errore ρ | quiete | ferita |
| --- | --- | --- | --- | --- | --- | --- |
| sfera 3D (r=6, R=60) | 841–853 | 59–66 | — | 1.3–1.8 | 3/3 | 3/3 riparate |
| disco 2D (r=4, R=60) | 216–221 | 55–65 | — | 0.8–1.2 | 3/3 | 3/3 riparate |
| guscio 3D (r=6, R=60, Rin=35) | 678–697 | 53–67 | 28–45 | 1.5–2.0 | 3/3 | 3/3 riparate |
| anello 2D (r=4, R=60, Rin=35) | 150–152 | 52–64 | 30–44 | 1.3–2.0 | 3/3 | 3/3 riparate |

Le forme cave sono più rumorose delle piene: il bordo interno è tenuto da una condizione sulla
*storia* della cellula (la giovinezza), non da una misura istantanea.

## Invarianza di scala

La stessa forma a raggio cellulare 6, 4.5 e 3 (3D). Tre unità di parametro distinte: lunghezze
della forma (fisse), lunghezze cellulari (scalano con r), tempi e contatti (fissi).
Una violazione documentata: **la larghezza del foro segue il raggio cellulare, non la lunghezza
della forma.**

## Gastrulazione 2D

Guscio sottile (Rin=45, parete di due cellule), arco apicale contiguo che si allarga, KC=0.5.

| fase | fronte (ρ) | canale aperto fino a (ρ) | pezzi |
| --- | --- | --- | --- |
| prima | — | 60–67 | 1 |
| fine costrizione | 36–42 | 45–50 | 1 |
| rilasciata (200–800 passi) | — | 59–66 (**richiusa**, 3/3) | 1 |
| `sempre mantieni` (4000 passi) | 36–42 invariato | 43–50 (**ferma**, 3/3) | 1 |
| dominio differenziato, costrizione spenta (4000 passi) | — | 43–54 (**ferma**, 3/3) | 1 |
| dominio ferito a metà | — | 59–64 (**perduta**, 2/2) | 1 |

Strain delle coppie costrette, stabile: −0.20 … −0.26. A KC=0.8 l'arco si autodistrugge; a KC=0.4
metà delle prove non parte.

## Gastrulazione 3D

Fronte da ρ≈57 a 38 e poi 31; canale aperto fino a ρ 43; tessuto in un pezzo; sotto `mantieni`
ferma 4000 passi. **Difetto dichiarato:** le cellule scendono da 643 a 418 e il guscio si
restringe (bordo 39–60 invece di 51–67) — la calotta, costringendosi, tira dentro tutto il guscio.
In 2D non succedeva.

## Ablazione dei meccanismi (anello 2D, 3 ripetizioni)

Con la polarità funzionante, sull'anello **statico** paga solo l'adesione apicale debole.

| variante | n | est | int | errore ρ | quiete |
| --- | --- | --- | --- | --- | --- |
| tutte le regole | 138–150 | 52–64 | 30–44 | 0.9–1.9 | 3/3 |
| senza polarità | 152–155 | 52–68 | 29–45 | 0.7–2.0 | 3/3, più rapida |
| senza divisione planare | 144–148 | 52–64 | 31–43 | 1.1–1.5 | 3/3, più rapida |
| senza deposizione basale | 146–152 | 51–64 | 32–43 | 0.6–1.0 | 3/3 |
| senza adesione apicale debole | 143–148 | 51–63 | 30–44 | 0.9–2.2 | **1/3** |

Polarità attiva: 60–67 cellule polarizzate su ~145 (cioè tutte quelle di superficie), orientamento
corretto — verso fuori sul bordo, verso il lume all'interno.

## La legge di simmetria

`elegge(biglietto = f)` rompe la simmetria se e solo se f è costante sulle orbite del gruppo
attuale. Misura: l'angolo fra i due poli, 8 ripetizioni per modo (2D).

| modo | biglietto di P2 | angolo fra i poli |
| --- | --- | --- |
| sequenza | d₁ (varia sull'orbita) | **170° ± 6°** |
| parallelo, biglietto identico | ρ, lo stesso per entrambe | **0° in 8 prove su 8** (stessa cellula) |
| parallelo, errore indipendente 0.3 | ρ + errore ≪ struttura | 38° (5/8 stessa cellula) |
| parallelo, errore indipendente 8 | ρ + errore ≫ struttura | 80° ± 48° |
| parallelo, errore indipendente 25 | ρ + errore ≫ struttura | 72° ± 58° |

Variazione residua di ρ sulle cellule di superficie, per confronto: deviazione 3.2, estremi 50–61.
Attesa per due punti casuali su un cerchio: 90° ± 52°.

**Previsione ancora da verificare:** la sequenza dà 170° invece di 180°, un difetto di ~10°
significativo rispetto all'errore della media. L'oggetto realizzato non ha il piano equatoriale: è
C∞v, non D∞h. Il difetto dovrebbe ritrovarsi come **asimmetria sistematica fra le due bocche del
foro del toro** — da misurare quando il toro sarà rifatto.
