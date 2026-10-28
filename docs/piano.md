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

## 5. Decisioni che servono

1. **Si comincia dal 2D laterale** per validare la fisica contro il conto esatto, sapendo che non
   tara il 3D?
2. **tyssue** come motore di confronto?
3. **Incomprimibilità**: vincolo esatto o penalità? Proposta: penalità, più semplice e standard,
   con una prova che il risultato non cambia alzando K.

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
