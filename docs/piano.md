# Piano: dalla specifica al motore

28 ottobre 2026. Fase 1 del piano di lavoro: colmare i buchi della specifica prima di scrivere il
motore. Tre buchi: **la previsione esatta** per il primo test fisico, **le unità e i parametri**,
**il motore esterno di confronto**. In fondo l'ordine dei lavori della Fase 2 e le decisioni che
servono dall'utente.

## 1. Il primo test fisico: cosa deve uscire

La specifica chiede di confrontare la curvatura con "la previsione analitica del modello a
vertici", ma non la scrive. Eccola, per il caso più semplice in cui è esatta.

### Il modello

Sezione di un monostrato (modello a vertici *laterale*, come in Krajnc e Ziherl 2015 e Polyakov et
al. 2014): ogni cellula è un quadrilatero con un lato apicale, uno basale e due laterali condivisi
con le vicine. Energia:

```
E = Σ_cellule (Γa · apicale + Γb · basale) + Γl · Σ_spigoli laterali lunghezza
```

con l'area di ogni cellula fissa ad A (incomprimibilità). Γa, Γb, Γl sono tensioni di linea
(forza; in 3D diventano tensioni di superficie).

### Il risultato

Si cerca il minimo fra le configurazioni in cui tutte le cellule sono lo stesso trapezio isoscele,
con l'apice dal lato concavo. Con Σ = Γa + Γb e Δ = Γa − Γb (l'eccesso di tensione apicale) il
conto è chiuso:

- **ogni cellula diventa un cuneo** con semiangolo **sin α = Δ / Γl**;
- **altezza** h = Γl · √(Σ·A) / (Γl² − Δ²)^¾, **larghezza a metà altezza** w = A / (h · cos α);
- **curvatura** del foglio κ = 2Δ / (Γl · w), cioè raggio R = Γl · w / (2Δ);
- **senza differenza** (Δ = 0) il foglio resta piatto, con rapporto altezza/larghezza
  **h/w = Σ / Γl**.

Due cose da notare. La prima: l'angolo α **non dipende dall'area**, che fissa soltanto la scala:
è l'invarianza di scala della specifica, sezione 9, ottenuta per costruzione. La seconda: il
conto vale solo finché il lato apicale resta positivo, cioè se

```
Γl² − Δ² > Δ · (Γa + Γb)
```

Oltre quel limite la cellula collassa a triangolo (apice zero) e il modello va regolarizzato con
una lunghezza minima di spigolo: il motore deve **rifiutare** quella configurazione, non
simularla in silenzio.

Con l'angolo **imposto** (un anello chiuso di N cellule, α = π/N) la stessa richiesta, lato
apicale positivo, diventa **sin(π/N) < Γl / (2Γa)**. È la soglia che la batteria verifica: con
Γa = 1.3 e Γl = 1 l'anello collassa per N ≤ 7 e regge per N ≥ 9.

**Verificato numericamente** (Julia, minimizzazione a forza bruta su 6 milioni di configurazioni,
A = 1): con Γ = (1.3, 1.0, 1.0) l'angolo misurato è 0.305 contro 0.305 previsto, l'altezza 1.628
contro 1.628; con Γ = (1.1, 1.0, 0.8) l'angolo è 0.125 contro 0.125, l'altezza 1.638 contro 1.639.
Con Γ = (1.5, 1.0, 1.2) i valori non tornano (0.416 contro 0.430): la condizione di validità è
violata (1.19 < 1.25) e la cellula è diventata un triangolo, esattamente come previsto. Il conto
va riscritto in `prove/` come test eseguibile.

### La scala delle prove

Il traguardo del README (una *fila* di cellule che costringe dentro un foglio piatto) non ha una
formula chiusa: le vicine non costrette resistono. Lo si raggiunge per gradini, ciascuno con un
numero atteso:

| prova | cosa | numero atteso |
| --- | --- | --- |
| F1 | foglio piatto, Δ = 0 | h/w = Σ/Γl, esatto |
| F2 | anello chiuso di N cellule (α = π/N imposto) | h = √(Σ·A / (cos α · (Γl − Δ·sin α))), esatto |
| F3 | striscia aperta, tutte le cellule costrette | sin α = Δ/Γl al centro; l'errore deve scendere al crescere di N (effetto dei bordi) |
| F4 | invarianza di scala | stesso α per A diverse, h e w che scalano come √A |
| F5 | **fila di k cellule costrette in un foglio piatto** (il traguardo) | per k grande l'angolo totale tende a 2kα; per k piccolo il confronto è con il motore esterno (sezione 3) |
| F6 | 3D: foglio di prismi esagonali, Δ = 0 | rapporto altezza/larghezza dal minimo del singolo prisma (Hannezo, Prost e Joanny 2014) |

Avvertenza dai risultati negativi (#4): **la sezione 2D non predice il 3D per la flessione.** Il 2D
qui serve a *validare la fisica* contro un conto esatto, non a tarare forme. F6 e il confronto
esterno in 3D vengono subito dopo.

Oltre a queste servono le prove di **meccanica pura**, che non dipendono dal modello biologico:

- le forze coincidono con le differenze finite dell'energia;
- l'energia scende sempre durante il rilassamento;
- la somma delle forze e dei momenti è zero (le forze sono appaiate);
- l'area o il volume di ogni cellula si conserva entro la tolleranza;
- il risultato non cambia dimezzando il passo di tempo.

## 2. Unità e parametri

Il risultato della sezione 1 dice una cosa importante: **la forma all'equilibrio dipende solo dai
rapporti delle tensioni** (Δ/Γl, Σ/Γl) e dall'area. La scala assoluta delle tensioni fissa solo la
*velocità* con cui ci si arriva. È una fortuna, perché le tensioni assolute sono il dato peggio
misurato: con l'ablazione laser si misura la velocità di ritrazione dopo il taglio, che è la
forza divisa per un attrito sconosciuto, e senza misurare quell'attrito la forza non si ricava.

Proposta:

| grandezza | unità | stato |
| --- | --- | --- |
| lunghezze | µm | la dimensione della cellula (A in 2D, volume in 3D) va presa da una misura per il tessuto scelto: **da cercare** |
| rapporti di tensione Δ/Γl, Σ/Γl | adimensionali | sono i parametri fisici veri; **liberi**, scelti entro la validità |
| tensione di riferimento Γl | nN | **libera**: fissa solo la scala dei tempi |
| attrito η (dinamica sovrasmorzata) | nN·min/µm | **libero**; insieme a Γl dà il tempo caratteristico τ = η·ℓ/Γl |
| tempo | min | si tara τ quando servirà, con le velocità di ritrazione misurate dopo ablazione: in cellule MDCK vanno da 0.01 a 1 µm/s (arXiv 1604.03733) |
| incomprimibilità | — | vincolo, oppure penalità K con prova che il risultato non dipende da K |

Restano da trovare con una fonte, prima di passare alle forme (non servono per F1–F5): la taglia
delle cellule nel tessuto di riferimento, i rapporti di tensione misurati in un'invaginazione vera
(il solco ventrale di *Drosophila*, Polyakov et al. 2014, è il candidato naturale), i coefficienti
di diffusione e degradazione per la chimica.

## 3. Il motore esterno di confronto

| motore | linguaggio | licenza | adatto? |
| --- | --- | --- | --- |
| **tyssue** (DamCB) | Python | GPL-3.0 | **sì**: è una libreria proprio per modelli a vertici di epiteli, in 2D e 3D, installabile da PyPI |
| Chaste (Oxford) | C++ | BSD-3 (dalla 3.1) | sì per il 2D a vertici, ma è pesante da installare e da usare |
| EmbryoMaker (Marin-Riera et al., *Bioinformatics* 2016) | — | dichiarato open source | **per ora no**: non ho trovato il repository né la licenza |

Proposta: **tyssue**. Lo si usa come programma separato dentro `confronto/`, senza mescolare il
suo codice al nostro, quindi la GPL non ci vincola. Le licenze di tyssue e Chaste vengono da
indici e pagine di progetto: vanno ricontrollate sul file LICENSE prima dell'uso.

## 4. Ordine dei lavori per la Fase 2

1. **Tipi e geometria 2D laterale**: vertici, cellule a quattro lati etichettati apicale, basale,
   laterale. Scritti in modo che il guscio 3D a due superfici sia la stessa astrazione con più
   facce, non un secondo motore.
2. **Energia e forze**: una sola funzione energia; le forze come suo gradiente, calcolato con
   differenziazione automatica (ForwardDiff.jl) oppure a mano e verificato con essa.
3. **Rilassamento sovrasmorzato** e le prove di meccanica pura.
4. **Prove F1–F4**, poi **F5** contro tyssue.
5. **Tag `motore-meccanica-v1`**: da qui la meccanica si tocca solo con un test fisico che fallisce.
6. Il guscio 3D (F6, poi l'invaginazione contro tyssue in 3D), gli scambi di vicine registrati, la
   chimica su griglia contro un profilo noto. Ognuno con il proprio tag quando è validato:
   aggiungere un modulo non è modificarne uno congelato.

## 5. Decisioni prese

1. Si comincia dal **2D laterale**.
2. **tyssue** come motore di confronto.
3. Incomprimibilità come **penalità**. In pratica K va tenuta morbida (0.01 nN/µm³ con aree di
   100 µm²): con K = 1 il passo di Eulero stabile diventa così piccolo che il rilassamento non
   finisce. L'area all'equilibrio scende allora del 10–20%, ma le previsioni valgono per l'area
   misurata, e la prova F4 verifica che cambiare K non sposta la forma.

## 6. Stato della Fase 2 (28 ottobre 2026)

Fatto e verificato (`scripts/prove.sh`, 368 prove su 368):

- geometria laterale, energia, forze esatte scritte a mano, rilassamento sovrasmorzato, misure;
- meccanica pura: forze uguali alle differenze finite, risultante e momento nulli, energia che
  scende a ogni passo, errore di Eulero del primo ordine;
- F1 (foglio piatto), F2 (anello, esatto fino a 1e-8), F3 (striscia costretta: al centro l'angolo
  coincide con sin α = Δ/Γl a meno di 1e-6), F4 (invarianza di scala e indipendenza da K);
- il rifiuto delle configurazioni degeneri, sulla soglia dell'anello;
- **il contatto fra parti lontane del tessuto** (`motore/05_contatto.jl`): un vertice e uno
  spigolo di cellule non vicine si respingono sotto d_contatto, come energia. Spento, una
  striscia che si arrotola oltre un giro si attraversa (56 incroci); acceso, si ferma quando
  tocca sé stessa. Se il tessuto si attraversa comunque, il motore lo rifiuta;
- **F5, il confronto con tyssue** (`confronto/`): stessa energia, due motori indipendenti, su
  una striscia piatta con una e con quattro cellule costrette e su un anello. Coincidono a meno
  di 1e-6 su angoli, altezze e aree, e l'energia fino a 1e-13 relativo.

Da tenere:

- **Piegarsi costa poco, quindi è lento.** Con Δ piccolo una striscia di 40 cellule non converge
  in 5 milioni di passi di Eulero. Risolto dall'integratore (sezione 7).
- tyssue si ferma con forze residue di 0.01 nN: per il confronto lo si fa ripartire fino a 1e-7.

La Fase 2 è chiusa: la meccanica è congelata come `motore-meccanica-v1`. Il tag git non si può
pushare da questo ambiente: al suo posto c'è il commit vuoto "CONGELATO: MOTORE-MECCANICA-V1" e,
soprattutto, la prova "file congelati", che confronta l'impronta SHA-256 dei cinque file con
`prove/congelati.txt` e fallisce se uno cambia.

## 7. L'integratore (fase 3, primo lavoro)

`motore/06_integratore.jl`, file nuovo: i cinque file congelati non sono toccati.

- **`minimizza!`: Newton smorzato** (Levenberg–Marquardt), con la matrice delle derivate
  seconde per differenze del gradiente esatto e un passo accettato solo se l'energia scende.
  Per gli equilibri. Una striscia di 40 cellule converge in 151 iterazioni (Eulero: non in 5
  milioni di passi), una di 80 in 598, l'anello in 16; l'angolo coincide con la teoria fino a
  1e-11. Limite dichiarato: la matrice è piena, va bene per il 2D con qualche centinaio di
  vertici, non per il 3D, dove servirà una versione sparsa.
- **`passo_implicito!` ed `evolvi!`: Eulero implicito** per la dinamica sovrasmorzata, scritto
  come spostamento minimizzante (il minimo di E(x) + η/(2dt)·|x − x_n|²). Stabile con qualunque
  passo: con dt = 50 min, oltre cento volte il limite di Eulero esplicito, l'energia non sale
  mai. Per dt → 0 ritrova la traiettoria di Eulero esplicito, con errore del primo ordine.
- FIRE, provato prima, è stato scartato: 13 volte più veloce di Eulero ma ancora lontano da
  Newton (una striscia di 80 cellule non converge in 2 milioni di iterazioni).

Le prove (`test/integratore.jl`) ritrovano con Newton F2, F3 e F5 (contro tyssue) e lo stesso
equilibrio di Eulero, e verificano il contatto. 459 prove su 459.

Congelato come `motore-integratore-v1` (commit "CONGELATO: MOTORE-INTEGRATORE-V1"), dopo averlo
reso indipendente dal tipo di tessuto, così il guscio 3D lo usa senza toccarlo.

## 8. Il guscio 3D (fase 3, secondo lavoro)

`motore/07_guscio.jl`, file nuovo. Ogni cellula è un prisma fra un poligono apicale e uno
basale, con facce laterali condivise; tensione di superficie per faccia (apicale, basale,
laterale), tensione di linea sull'anello apicale, penalità sul volume. Usa l'integratore
congelato senza toccarlo.

Le previsioni esatte, per prismi esagonali regolari a volume fissato:

- **in un foglio** (facce laterali condivise): **altezza/lato = √3·(σa + σb)/σl**;
- **isolato** (σa = σb): lato s con s³ = 3σl·V/((σa + σb)·c²), c = 3√3/2, cioè metà del rapporto
  del foglio.

Verificate (`test/guscio.jl`, 643 prove su 643 in tutto):

- geometria esatta su prismi regolari, forze uguali alle differenze finite (a meno di 1e-7) e
  appaiate;
- F6a, prisma isolato, e F6b, **foglio periodico**: entrambe coincidono con la formula a meno
  di 1e-9. Newton converge in 10 iterazioni;
- F6c invarianza di scala e indipendenza da K;
- una cellula che alza l'anello apicale (Λ = 2 nN) si restringe all'apice da 50 a 15 µm².

Il foglio periodico non ha bordi: le cellule a destra toccano quelle a sinistra, e le due
dimensioni del periodo sono incognite come le posizioni. È la prova pulita, perché:

- **un pezzo di foglio con bordo libero non ha equilibrio** senza scambi di vicine: le facce
  laterali del bordo, di una cellula sola, accorciano gli spigoli fino alla lunghezza minima.
  Registrato come risultato negativo #13;
- **con σl < σa una cellula isolata si arrotonda**: lascia la forma di prisma, tre spigoli
  apicali su sei collassano e l'energia scende sotto quella del miglior prisma regolare. È
  fisica (una cellula libera si arrotonda), ma richiede di nuovo uno scambio topologico.

Con un anello apicale forte (Λ = 5 nN) l'apice si chiude del tutto: in una cellula vera è
l'estrusione dal foglio, che il modello non ha ancora, e il motore si ferma prima.

### F7, il confronto 3D con tyssue

tyssue non sa fare fogli periodici in 3D, quindi il confronto usa un pezzo di 37 prismi con il
**bordo bloccato** (vertici fissi, una novità di `07_guscio.jl`), in due scenari: la cellula
centrale che stringe l'anello apicale e una placca di sette cellule. Esito
(`confronto/LEGGIMI.md`):

- le due energie coincidono fino a 1e-12 relativo anche su posizioni deformate;
- il gradiente di tyssue non è la derivata esatta della sua energia (1% nel termine di volume),
  e per questo il suo minimizzatore non converge; il nostro sì, in 7–8 iterazioni, e trova un
  minimo più basso del punto in cui tyssue si ferma.

Trovato strada facendo, nel nostro motore: con i vertici bloccati, Newton li sentiva ancora
nella matrice delle derivate seconde e il passo li spostava, così il rilassamento si
trascinava. Ora la posizione di un vertice bloccato si legge da una copia fissa, e per
l'integratore non è un'incognita.

Congelato come `motore-guscio-v1` (commit "CONGELATO: MOTORE-GUSCIO-V1"), dopo aver separato
la costruzione dalla topologia (`tessuto_da_poligoni`, `poligoni`), così i cambi di topologia
ricostruiscono facce e cellule senza toccarlo.

## 9. Scambi di vicine (fase 3, terzo lavoro)

`motore/08_scambi.jl`, file nuovo (specifica, decisione 4). Quando lo spigolo apicale fra due
cellule A e B scende sotto una soglia, diventa uno spigolo perpendicolare fra le altre due
cellule C e D che lo toccavano; lo stesso sulla superficie basale, e la faccia laterale fra A e
B diventa una faccia fra C e D. Ogni scambio è registrato. Al bordo bloccato non è ammesso.

Verificato (`test/scambi.jl`, 693 prove su 693 in tutto):

- la topologia: A e B perdono un lato, C e D ne guadagnano uno, ogni vertice resta di tre
  cellule, i volumi positivi, le forze esatte;
- lo scambio è reversibile: rifatto sullo stesso spigolo riporta alla stessa rete di cellule;
- l'energia è continua: con lo spigolo lungo ε il salto vale circa 49·ε nN·µm, e il volume
  totale non cambia;
- `scambi_di_vicine!` scambia gli spigoli sotto soglia e solo quelli;
- **un difetto da scambio forzato è metastabile**: dopo uno scambio forzato nel reticolo
  regolare il tessuto si ferma con due cellule a 5 lati e due a 7, 24 nN·µm sopra il reticolo
  regolare, e non torna indietro da solo; con lo scambio inverso esplicito sì.

- lo scambio parte **solo se è corta tutta la parete**, apicale e basale: nel guscio a prismi
  una cellula ha gli stessi vicini sopra e sotto. Due cellule che stringono l'anello apicale
  chiudono la giunzione apicale ma non quella basale, e lo scambio fatto lo stesso fa saltare
  l'energia di 8000 nN·µm: risultato negativo #14 (le cellule vere qui diventano scutoidi);
- **scambi spontanei**: un pezzo di tessuto fluido (σl = 0.4) ancorato al bordo e stirato
  (estensione convergente) risponde in modo elastico fino a s = 2.1; a s = 2.15 parte una
  valanga di scambi e l'energia scende di 170 nN·µm sotto lo stato elastico. È una soglia di
  snervamento: sotto il tessuto è solido, sopra scorre. La valanga arriva al bordo bloccato, dove
  lo scambio non è ammesso, e lì si ferma.

Congelato come `motore-scambi-v1` (commit "CONGELATO: MOTORE-SCAMBI-V1"). Le operazioni che
mancano (estrusione, scambio al bordo) entrano come file nuovi.

## 10. Estrusione (fase 3, quarto lavoro)

`motore/09_estrusione.jl`, file nuovo (specifica, decisione 4, e il repertorio: "morire:
apoptosi, oppure estrusione dal foglio"). Una cellula triangolare con apice e base tutti e due
sotto una soglia d'area esce dal foglio: i suoi tre vertici apicali si fondono in uno, così i
tre basali, e le tre vicine si chiudono su quel vertice. È l'operazione T2 dei modelli a
vertici, sulle due superfici del guscio. Ogni estrusione è registrata.
`rilassa_con_topologia!` alterna Newton, scambi ed estrusioni.

Verificato (`test/estrusione.jl`, 739 prove su 739 in tutto):

- **una cellula apoptotica esce dal foglio da sola**: con il volume bersaglio portato da 1000 a
  1 µm³ si restringe, perde tre lati con tre scambi di vicine, diventa un triangolo ed esce; il
  tessuto si richiude (da 36 a 35 cellule, due vertici in meno) e converge;
- interruttore: senza operazioni topologiche la stessa cellula si blocca, esagonale, e il
  rilassamento non converge;
- l'estrusione è un evento discreto che libera energia: le pareti del triangolo, alte ancora
  ~15 µm, spariscono con lui (caduta dello 0.3% dell'energia totale), e dopo il tessuto si
  richiude abbassandola ancora;
- **risultato negativo #15**: una cellula che chiude solo l'apice (0.015 µm²) con la base larga
  (175 µm²), la cellula a bottiglia, non si può estrudere nel guscio a prismi. È lo stesso limite
  del #14: apice e base cambiano insieme.

Congelato come `motore-estrusione-v1` (commit "CONGELATO: MOTORE-ESTRUSIONE-V1").

## 11. Il bordo libero, e il limite comune

Il terzo lavoro previsto, lo scambio al bordo libero, non si può fare nel guscio a prismi. La
diagnosi (risultato negativo #13, aggiornato): le pareti libere delle cellule di bordo si
torcono, e collassano a turno all'apice o alla base, mai insieme. Ogni operazione sui prismi
cambia insieme le due superfici e schiaccerebbe uno spigolo di 8 µm.

I risultati #13, #14 e #15, e la cellula isolata che si arrotonda, hanno **la stessa radice**:
nel guscio a prismi una cellula ha lo stesso poligono all'apice e alla base. Le cellule vere no
(gli scutoidi). L'estensione che li risolve tutti è un guscio con **superfici apicale e basale a
topologia propria**: ogni superficie fa i suoi scambi, e le facce laterali non sono più
quadrilateri ma poligoni che collegano lo spigolo apicale e quello basale di due cellule
vicine, anche quando sono vicine solo da una parte. È un cambio di rappresentazione, da
decidere prima di scriverlo.

## 12. Scutoidi: la cellula come poliedro qualunque (fase 3, quinto lavoro)

Decisione presa il 29 ottobre: risolvere subito il limite comune dei risultati #13, #14 e #15.
Una cellula resta un **monostrato** (va dalla superficie apicale a quella basale), ma è un
**poliedro qualunque**, non più un prisma. È il compromesso fra il guscio a prismi e il modello
3D pieno che la specifica esclude (cellule che si rimescolano in volume): le operazioni
topologiche avvengono sulle superfici e lungo le pareti, non nel volume.

Energia, aree, volumi e forze del guscio (`07_guscio.jl`) lavorano già su facce poligonali
qualsiasi: a supporre i prismi erano solo la costruzione e le operazioni topologiche. Gli
scutoidi entrano quindi come operazioni nuove, in `motore/10_scutoidi.jl`.

**Fatto: lo scambio apicale con scutoide.** Quando si chiude solo la giunzione apicale fra A e
B: in alto lo spigolo passa a C e D; la parete fra A e B resta e in alto finisce in un vertice
nuovo, interno al foglio, dove si incontrano A, B, C e D; sopra quel vertice nasce una parete
triangolare fra C e D. Verificato (`test/scutoidi.jl`): ogni cellula resta un poliedro chiuso
(V − E + F = 2), volumi positivi, forze esatte, volume totale conservato, energia continua (il
salto va a zero come ε).

**Trovato per strada, e corretto in un file congelato** (regola 1, con il test fisico che
fallisce): aggiungere un punto su uno spigolo cambiava aree e volumi di facce non piane, perché
il centro di triangolazione era il baricentro dei vertici. Ora è il baricentro del contorno,
che non cambia suddividendo uno spigolo. Prima della correzione lo scambio apicale faceva
saltare l'energia di 183 nN·µm. Conseguenza su F7: con tyssue, che usa il baricentro dei
vertici, le energie coincidono solo dove le facce sono piane (sezione 8).

**Che cosa manca.** In un tessuto uniforme uno scutoide forzato torna indietro: il vertice
interno risale all'apice e lo spigolo nuovo fra C e D si riaccorcia. Servono le operazioni che
completano il quadro, cioè le due operazioni elementari dei modelli a vertici 3D (Okuda et al.,
"reversible network reconnection"):

1. **da spigolo a triangolo**: uno spigolo corto, condiviso da tre facce, diventa una faccia
   triangolare (lo scambio apicale ne è il caso sulla superficie);
2. **da triangolo a spigolo**: una faccia triangolare piccola diventa uno spigolo. È l'inversa,
   e serve a far tornare indietro lo scutoide, o a completarlo in uno scambio intero quando il
   vertice interno arriva alla base.

Con queste due, sulle superfici e nel volume del foglio, si risolvono anche il bordo libero
(#13) e la cellula che si arrotonda; l'estrusione apicale della cellula a bottiglia (#15)
chiede in più che una cellula possa lasciare una delle due superfici.

## 13. Le operazioni elementari generali (fase 3, sesto lavoro; roadmap R1)

`motore/11_operazioni.jl`, fatto in red-green: prima le dieci prove in `test/operazioni.jl`
(commit rosso), poi il codice (commit verde, 30 ottobre).

Il mezzo esterno conta come una pseudo-cellula per ogni etichetta di faccia libera (apicale −1,
basale −2, laterale −3). Così attorno a ogni spigolo ci sono tre regioni e attorno a ogni
vertice quattro, e le operazioni sono le stesse in superficie, al bordo e dentro il foglio:

- `spigolo_a_triangolo!` (I→H): lo spigolo con attorno A, B, E e agli estremi in più C e D
  diventa un triangolo fra C e D, nel piano perpendicolare allo spigolo, di lato
  `lunghezza_nuova`; ogni vertice nuovo va verso la faccia da cui nasce. Sullo spigolo apicale
  dà la stessa topologia dello scambio apicale di `10_scutoidi.jl`;
- `triangolo_a_spigolo!` (H→I): l'inversa, con lo spigolo lungo la normale del triangolo;
- `contrai_spigolo!`: sul contorno libero di una cellula di bordo (la cellula e due mezzi) I→H
  non è ammessa, perché creerebbe una faccia fra due mezzi; lì lo spigolo corto si contrae e
  l'angolo si fonde nell'altro estremo. Le prove l'hanno chiesta: senza, il foglio a bordo libero
  restava fermo con uno spigolo del contorno alla lunghezza minima, fra un angolo e un vertice
  dove arriva una vicina;
- `problemi_topologici`: gli invarianti. Al bordo libero valgono due eccezioni: il vertice
  d'angolo ha tre regioni, e lo spigolo dove la superficie libera di una cellula fa angolo
  (fra due sue facce verso il mezzo laterale) ne ha due.

Risultati (prove in `test/operazioni.jl`):

| prova | esito |
| --- | --- |
| invarianti dei fogli (a bordo bloccato, periodico, dopo lo scambio apicale) | nessuna violazione |
| invarianti violati (una faccia tolta) | visti |
| I→H sullo spigolo apicale = `scambio_apicale!` | stessa topologia |
| H→I dopo I→H | topologia di partenza |
| continuità, ε = 10⁻² e 10⁻³ | salto di energia ∝ ε (rapporto fra 9 e 11), volume conservato a 10⁻¹² |
| forze esatte dopo I→H e dopo H→I | errore sotto 10⁻⁷ |
| scambio intero in due passi | chi-tocca-chi dello scambio fra prismi di `08` |
| scutoide forzato | torna indietro, stesso vicinato ed energia del reticolo regolare |
| foglio di 19 prismi a bordo libero (#13) | equilibrio in 10 cicli, 24 contrazioni, 7747.9 nN·µm |
| cellula isolata, σl = 0.9 | equilibrio in 2 cicli, 6 contrazioni (esagoni → triangoli), 455.92 nN·µm, sotto i 456.03 dove si fermava senza operazioni |

La contrazione lascia fermo l'estremo che resta, quindi il salto di energia è proporzionale
alla lunghezza dello spigolo, come per le altre due.

**Fogli periodici** (31 ottobre, ancora red-green: cinque prove in
`test/operazioni_periodiche.jl`, poi il codice). In un foglio periodico un vertice compare in
ogni faccia con un'immagine. Prima di un'operazione i vertici coinvolti si portano in una
cornice comune con due ricalibrazioni che non cambiano la geometria: un vertice si sposta di un
periodo e le sue immagini scendono di uno; una faccia abbassa le sue immagini e alza della
stessa quantità la sua traslazione in ogni cellula. Dopo, nelle facce toccate i vertici
coinvolti hanno immagine (0, 0), e l'operazione è quella di prima. Gli spigoli si distinguono
anche per lo scarto fra le immagini dei due estremi (anche nel conto della tensione di linea);
un'operazione rifiutata lascia il tessuto com'era. Provato nel foglio 4 × 4, su uno spigolo
interno alla scatola e su uno che ne attraversa il bordo: andata e ritorno, forze esatte anche
rispetto alla scatola, continuità, scambio intero in due passi uguale a quello di `08`, scutoide
forzato che torna indietro all'energia del reticolo regolare. Batteria: 781 prove passate.

`10` e `11` sono congelati insieme (`motore-scutoidi-v1`, `motore-operazioni-v1`).
`scambio_apicale!` resta solo non periodico: nei fogli periodici fa lo stesso
`spigolo_a_triangolo!` sullo spigolo apicale.

## Fonti

- Krajnc e Ziherl, *Theory of epithelial elasticity*, Phys. Rev. E 92, 052713 (2015):
  https://link.aps.org/doi/10.1103/PhysRevE.92.052713
- Hannezo, Prost e Joanny, *Theory of epithelial sheet morphology in three dimensions*, PNAS 111,
  27–32 (2014): https://pmc.ncbi.nlm.nih.gov/articles/PMC3890844
- Polyakov et al., *Passive mechanical forces control cell-shape change during Drosophila ventral
  furrow formation*, Biophys. J. 107, 998–1010 (2014): https://pmc.ncbi.nlm.nih.gov/articles/PMC4142243
- Velocità di ritrazione dopo ablazione in MDCK: https://arxiv.org/pdf/1604.03733
- tyssue: https://pypi.org/project/tyssue e https://tyssue.readthedocs.io/
- Chaste: https://www.cs.ox.ac.uk/chaste
- EmbryoMaker, *Bioinformatics* 32(2), 2016, doi:10.1093/bioinformatics/btv527

Il conto della sezione 1 è mio, non preso da un articolo: la verifica numerica è quella riportata
sopra, e va rifatta come prova nel motore.
