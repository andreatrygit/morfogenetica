# Roadmap

29 ottobre 2026. Questo è il file da cui si riparte quando si apre una sessione nuova. Dice che
cosa stiamo facendo, dove siamo, come si lavora e che cosa viene dopo, passo per passo. La storia
dettagliata di ogni lavoro fatto e dei numeri misurati è in `docs/piano.md`; i fallimenti, che
non vanno ripetuti, in `docs/risultati-negativi.md`.

**Quando un passo è fatto, aggiornare questo file nello stesso commit** (la sezione "Dove siamo"
e il passo stesso).

## 1. Che cosa stiamo facendo

Non un simulatore generico di cellule. Il progetto è una disciplina morfogenetica: N cellule
portano ciascuna una copia dello stesso programma (una tabella di fasi), leggono solo quattro
canali locali (dentro di sé, dalle vicine, dalla matrice, misure assolute), usano un repertorio
chiuso di comportamenti chiedibili a una cellula vera, e da questo devono nascere forme, senza
seme e senza piano centrale (`README.md`, `docs/specifica-motore-julia.md`,
`docs/algebra-delle-primitive-morfogenetiche.md`).

Il motore in Julia è il **substrato fisico** di quella disciplina: meccanica a energia, geometria
che possa essere un cuneo, topologia controllata, chimica. La specifica chiede di validarlo e
congelarlo *prima* di costruirci sopra, per non ripetere l'errore del motore JavaScript, ritoccato
a ogni forma. Finora si è lavorato solo sul substrato meccanico e topologico; **lo strato che dà
senso al progetto (percezione, fasi, primitive, programmi) non c'è ancora**. Tenerlo presente:
ogni passo in più sulla meccanica rimanda la prima forma.

## 2. Dove siamo

Moduli in `motore/`, inclusi a mano in `src/Morfogenetica.jl` (la cache di Julia non vede i file
nuovi di una cartella). I congelati hanno l'impronta SHA-256 in `prove/congelati.txt`: la prova
"file congelati" fallisce se cambiano.

| file | cosa | stato |
| --- | --- | --- |
| `01`–`05` | modello a vertici laterale 2D (sezione del foglio): geometria, energia, Eulero, misure, contatto | congelati, `motore-meccanica-v1` |
| `06_integratore.jl` | Newton smorzato per gli equilibri, Eulero implicito per la dinamica; generico sul tipo di tessuto | congelato, `motore-integratore-v1` |
| `07_guscio.jl` | guscio 3D: cellule poliedriche fra superficie apicale e basale, foglio periodico, bordo bloccato | congelato, `motore-guscio-v1`; modificato il 29/10 con test fisico (centro delle facce = baricentro del contorno) |
| `08_scambi.jl` | scambi di vicine fra prismi (T1 su tutta la parete) | congelato, `motore-scambi-v1` |
| `09_estrusione.jl` | estrusione di una cellula triangolare piccola (T2), `rilassa_con_topologia!` | congelato, `motore-estrusione-v1` |
| `10_scutoidi.jl` | scambio apicale con scutoide (vertice interno al foglio), solo non periodico | congelato, `motore-scutoidi-v1` |
| `11_operazioni.jl` | operazioni elementari generali: spigolo↔triangolo, contrazione al bordo libero, invarianti; anche nei fogli periodici | congelato, `motore-operazioni-v1` |
| `12_chimica.jl` | una sostanza che diffonde su griglia 3D, con degradazione, secrezione e assorbimento delle cellule | **in corso, non congelato** (R3 verde) |

Prove: `scripts/prove.sh`, 795 passate e nessuna rotta al 31 ottobre, circa 6 minuti (il grosso è Newton a
matrice piena sui tessuti 3D). Confronti esterni con tyssue in `confronto/` (`LEGGIMI.md`).

Il tag git non si può pushare da questo ambiente: un modulo congelato è segnato da un **commit
vuoto** "CONGELATO: MOTORE-…-V1" e dall'impronta in `prove/congelati.txt`.

## 3. Come si lavora

- **Red-green.** Prima le prove, in un commit "rosso" dove sono marcate `@test_broken` (falliscono
  senza bloccare la batteria e risultano come rotte); poi il codice, nel commit "verde" che le
  fa passare e toglie `@test_broken`. Dove un numero si trova solo con l'esperimento (una soglia,
  un tempo), prima si scrivono le prove delle proprietà che devono valere comunque (continuità,
  conservazione, invarianti, forze esatte, reversibilità) e il numero diventa poi una prova di
  regressione, dichiarata come tale.
- **Congelamento.** Un modulo validato si congela: impronta in `prove/congelati.txt` più commit
  vuoto. Dopo, si tocca solo con un test fisico che fallisce, nello stesso commit dell'impronta
  nuova, con la ragione scritta accanto all'impronta.
- **Regole della specifica** (`CLAUDE.md`): quattro canali, nessun burattinaio, comportamenti
  chiedibili a una cellula vera, ogni meccanismo con un interruttore, casualità con seme, unità
  vere, risultati negativi come regressioni.
- **Lingua**: tutto in italiano, codice compreso.
- **Git**: si lavora sul branch indicato dalla sessione; ogni lavoro finisce in una PR verso
  `main`, che l'utente unisce (o chiede di unire). Unire con merge, non squash, così i commit
  vuoti di congelamento restano nella storia.
- **Ambiente**: Julia via juliaup (`CLAUDE.md`, sezione "Installare Julia"); per rigenerare i
  confronti con tyssue serve Python 3.12 con NumPy < 2 (`confronto/LEGGIMI.md`).

## 4. I passi

### R1 — Le due operazioni elementari, in forma generale (fatto)

Stato: **fatto e congelato** (31 ottobre). Il codice è in `motore/11_operazioni.jl`; le prove in
`test/operazioni.jl` (dieci, 30 ottobre) e `test/operazioni_periodiche.jl` (cinque, fogli
periodici, 31 ottobre), tutte scritte prima del codice. `10` e `11` sono congelati insieme.

Che cosa c'è:

- `problemi_topologici(t)`: gli invarianti (V − E + F = 2 per cellula, tre regioni per spigolo,
  quattro per vertice, facce di almeno tre vertici, nessuna faccia doppia, una faccia apicale e
  una basale per cellula). Due eccezioni al bordo libero, dove una cellula sola tocca due mezzi:
  il vertice d'angolo ha tre regioni, e lo spigolo dove la superficie libera di una cellula fa
  angolo ne ha due (la cellula e il mezzo laterale);
- `spigolo_a_triangolo!` (I→H) e `triangolo_a_spigolo!` (H→I);
- `contrai_spigolo!`: una terza operazione, che le prove hanno chiesto. Sul contorno libero
  (lo spigolo fra la cellula, il mezzo laterale e quello apicale o basale) I→H non si può fare,
  perché creerebbe una faccia fra due mezzi; lo spigolo che si accorcia lì si contrae, l'angolo
  si fonde nell'altro estremo. È ciò che risolve il risultato negativo #13;
- `rilassa_con_operazioni!`: Newton e operazioni alternati, fino all'equilibrio;
- i fogli periodici: prima di un'operazione i vertici coinvolti si portano in una "cornice"
  comune, ricalibrando vertici e facce (immagini e traslazioni compensate, geometria identica);
  poi l'operazione lavora come nel caso non periodico. Se la scatola è così piccola che la
  cornice non esiste, l'operazione non si fa. Un'operazione rifiutata lascia il tessuto com'era.

Numeri (prove di regressione, non bersagli): il foglio di 19 prismi a bordo libero arriva
all'equilibrio in 10 cicli con 24 contrazioni, energia 7747.9 nN·µm; la cellula isolata in 2
cicli con 6 contrazioni (l'esagono apicale e quello basale diventano triangoli), energia 455.92
nN·µm, sotto i 456.03 dove si fermava senza operazioni, con tre spigoli schiacciati alla
lunghezza minima (`test/guscio.jl`).

Perché: risolvono insieme i risultati negativi #13 (bordo libero), #14 (giunzione che si chiude
solo all'apice), la cellula isolata che si arrotonda, e lo scutoide che vuole tornare indietro;
e sono la topologia di un modello 3D a vertici completo, quindi preparano il multistrato (R7)
senza costruirlo.

Le due operazioni dei modelli 3D a vertici (Okuda et al., "reversible network reconnection"):

1. **da spigolo a triangolo** (I→H): uno spigolo corto, attorno a cui stanno tre "cellule" A, B,
   E, i cui estremi hanno in più C e D, diventa una faccia triangolare fra C e D; i tre vertici
   del triangolo stanno fra C, D e due di A, B, E;
2. **da triangolo a spigolo** (H→I): una faccia triangolare piccola fra C e D diventa uno
   spigolo attorno a cui stanno le altre tre. È l'inversa.

Scelta di progetto: il **mezzo esterno come pseudo-cellula** (una per ogni mezzo: apicale,
basale, laterale libero). Così ogni spigolo ha esattamente tre "cellule" attorno e ogni vertice
quattro, e le stesse due operazioni valgono sulla superficie, al bordo e nel volume. Lo scambio
apicale di `10_scutoidi.jl` è I→H su uno spigolo apicale (E = mezzo apicale). Vincolo: un'operazione
non deve mai creare una faccia fra due mezzi (aprirebbe un buco nel foglio).

Le prove, scritte prima del codice (commit rosso) e ora verdi:

- invarianti dopo ogni operazione, in `problemi_topologici`: ogni cellula ha
  V − E + F = 2; ogni spigolo ha tre "cellule" attorno; ogni vertice quattro; facce di almeno tre
  vertici; nessuna faccia doppia; monostrato (ogni cellula tocca le due superfici) quando il
  tessuto lo dichiara;
- I→H su uno spigolo apicale dà la stessa topologia di `scambio_apicale!`;
- H→I è l'inversa di I→H: stessa topologia di partenza;
- continuità dell'energia (salto che va a zero come ε) e volume totale conservato, nei due
  versi; forze esatte dopo ogni operazione;
- **lo scambio intero come due passi**: I→H sullo spigolo apicale, poi H→I sul triangolo A|B
  quando lo spigolo basale si accorcia, dà la stessa topologia dello scambio fra prismi di
  `08_scambi.jl` (sul pezzo di foglio a bordo bloccato: i fogli periodici non sono ancora
  supportati);
- lo scutoide forzato in un tessuto uniforme torna indietro e il tessuto ritrova il reticolo
  regolare (con H→I);
- regressioni riviste: il foglio a bordo libero (#13) e la cellula isolata arrivano a un
  equilibrio. Il #14 (la giunzione che si chiude solo all'apice) è coperto dallo scambio in due
  passi e dallo scutoide che torna indietro; una prova dinamica, con l'anello apicale che
  stringe, verrà con R4.

Nei fogli periodici (`test/operazioni_periodiche.jl`), su uno spigolo interno alla scatola e
su uno che ne attraversa il bordo: andata e ritorno con gli invarianti, forze esatte (anche
rispetto alla scatola), continuità, scambio intero in due passi uguale a quello di `08`, scutoide
forzato che torna indietro.

Lasciato com'è: `scambio_apicale!` di `10` resta solo non periodico; nei fogli periodici si usa
`spigolo_a_triangolo!`, che fa la stessa cosa.

### R2 — Aggiornare la specifica sul 3D pieno (fatto)

La specifica (sezione 3) escludeva il modello 3D pieno. Decisione del 29 ottobre: non escluso,
**rimandato**, con la strada preparata da R1. Scritto nella specifica; resta da riportarlo nel
quaderno dell'algebra quando lo si riprende.

### R3 — Chimica (specifica, decisione 5) (in corso)

Stato: **commit verde fatto** (31 ottobre). Le dieci prove di `test/chimica.jl`, scritte prima,
passano; il codice è in `motore/12_chimica.jl`: griglia cartesiana 3D a voxel (passo per asse, bordi chiusi o
periodici), schema esplicito con sotto-passi scelti per la stabilità, degradazione e secrezione
integrate esattamente nel passo (così il bilancio stazionario non ha errore di passo), secrezione e assorbimento alla posizione con pesi trilineari (assorbimento
esponenziale, così la concentrazione non diventa negativa), `centro_cellula` per il guscio.
Parametri di riferimento dalla letteratura: Dpp nel disco dell'ala, D = 0.1 µm²/s = 6 µm²/min e
λ ≈ 20 µm (Kicheva et al., Science 315, 521, 2007), quindi k = D/λ² = 0.015/min.

Numeri (`docs/piano.md`, sezione 14): λ misurata con errore 2·10⁻⁴ a h = 1 µm e 8·10⁻⁴ a
h = 2 µm (secondo ordine); varianza 6Dt esatta; sorgente puntiforme entro il 3% a 5–8 µm;
errore nel tempo del primo ordine (rapporto 2.07 dimezzando il passo).

Resta, prima di congelare `12`: la griglia deve seguire la scatola quando la scatola cambia
(oggi il periodo è fissato alla costruzione; nei fogli periodici la scatola è un'incognita di
Newton). Si decide con R4, quando chimica e meccanica girano nella stessa passata. La
regressione del #1 chiede la percezione e le fasi (R4).


Almeno una sostanza che diffonde, con sorgenti, consumo e degradazione, risolta come PDE su
griglia con passo proprio, accoppiata alle posizioni delle cellule (secrezione e lettura alla
propria posizione). Prove analitiche prima: profilo di una sorgente puntiforme, stato stazionario
esponenziale con lunghezza √(D/k), conservazione della massa senza degradazione, convergenza con
la griglia e con il passo. Coefficienti di diffusione e degradazione da letteratura, con la fonte.
Regressione: il risultato negativo #1 (la necrosi da limite di diffusione non stabilizza un lume)
va riprovato con la chimica vera.

### R4 — Percezione e programma (specifica, sezioni 5–8)

La struttura di percezione con i soli quattro canali, costruita senza accesso allo stato globale
(il tipo lo impedisce); la tabella di fasi copiata in ogni cellula; l'avanzamento monotono
(evento di fine rilevato da sé, oppure una vicina già nella fase dopo); le intenzioni meccaniche
come uscita di una fase (tensioni, volume bersaglio, richieste di dividersi, morire,
differenziare). Prima primitiva: `costringi` come anello apicale. I tre burattinai del codice
vecchio (sezione 8 della specifica) diventano prove che non devono ripresentarsi.

### R5 — Prime forme

L'invaginazione 3D di una placca che stringe l'anello apicale (oggi possibile: l'apice si
stringe e il foglio si piega), poi la catena di gastrulazione compilata (specifica, sezione 5).
Confronto con `prove/attesi.md` sui rapporti e sui comportamenti, non sui valori.

### R6 — Infrastruttura (quando una forma la chiede)

- **Solutore sparso**: Newton con matrice piena non regge migliaia di cellule; servono la matrice
  delle derivate seconde sparsa e un solutore adatto. Già ora la batteria dura 5 minuti:
  separare prove veloci e lente.
- **Contatto in 3D**: oggi solo nel modello laterale; senza, un foglio che si piega può
  attraversarsi.
- Le operazioni topologiche nei fogli periodici.

### R7 — Multistrato e mesenchima

Con R1 fatto, la struttura dati e la topologia sono già quelle di un modello 3D completo: manca
togliere il vincolo di monostrato, generalizzare le etichette (la tensione di una faccia dipende
da che cosa c'è dall'altra parte, cioè dalla coppia di tipi: è già la decisione sull'adesione
della specifica), costruire tessuti pieni, e decidere se le cellule mesenchimali (nella specifica,
centri con vicinato di Delaunay) diventano anch'esse poliedri o restano un secondo tipo. Serve un
confronto esterno per il 3D in volume: tyssue ha un modello bulk, ma il suo gradiente non è
esatto (`confronto/LEGGIMI.md`).

## 5. Decisioni aperte

- R7: cellule mesenchimali come poliedri o come centri.
- Quale motore esterno per il confronto del 3D in volume.
- I valori dei parametri dalla letteratura (dimensioni delle cellule, rapporti di tensione,
  coefficienti di diffusione): oggi sono quasi tutti dichiarati liberi (`docs/piano.md`, sezione 2).

## 6. Per riprendere

1. Leggere questo file, poi `CLAUDE.md`.
2. Installare Julia se manca e lanciare `scripts/prove.sh`: deve passare tutto, salvo le prove
   marcate `@test_broken` del passo in corso.
3. Riprendere dal primo passo non finito della sezione 4, con il metodo della sezione 3.
