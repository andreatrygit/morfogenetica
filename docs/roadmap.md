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
| `11_operazioni.jl` | operazioni elementari generali: spigolo↔triangolo, contrazione al bordo libero, invarianti; anche nei fogli periodici | congelato, `motore-operazioni-v1`; corretto il 2/11 con test che fallivano (facce che restavano con due vertici) |
| `12_chimica.jl` | una sostanza che diffonde su griglia 3D, con degradazione, secrezione e assorbimento delle cellule; la griglia segue la scatola | congelato, `motore-chimica-v1` |
| `13_programma.jl` | percezione a quattro canali, tabella di fasi, onda, passata completa (popolazione, matrice, casualità con seme) | congelato, `motore-programma-v1`; modificato il 2/11 con test che fallivano (le fasi dichiarano che cosa scrivono) |
| `14_popolazione.jl` | divisione nel guscio a prismi (angolo dato o regola di Hertwig) | congelato, `motore-popolazione-v1` |
| `15_sparso.jl` | il solutore sparso per i tessuti di migliaia di cellule | congelato, `motore-sparso-v1`; esteso il 5/11 (`hessiana = :esatta`) |
| `16_hessiana.jl` | la matrice delle derivate seconde esatta, per il solutore sparso | congelato, `motore-hessiana-v1` |
| `primitive/costringi.jl` | la prima primitiva: anello apicale teso dove la sostanza supera la soglia | fatta (le primitive non si congelano: si ridisegnano) |
| `primitive/secerni.jl` | secernere un segnale se marcata | fatta |
| `primitive/attiva_inibisci.jl` | attivatore e inibitore (Nodal/Lefty): l'organizzatore senza seme | fatta |
| `primitive/induci.jl` | un gene bersaglio che si accumula con il segnale | fatta |
| `programmi/composizione.jl` | `insieme` (‖, rifiuta i conflitti) e `con_fine` | fatta |
| `programmi/invaginazione.jl` | la prima forma: l'invaginazione di una placca | fatta (R5) |
| `programmi/gastrula.jl` | la catena a scala vera: l'organizzatore senza seme, poi la sua invaginazione | fatta (R5) |

Prove: `scripts/prove.sh`, la batteria intera: 861 passate e una rotta dichiarata
(`prova_gradiente_veloce`) al 4 novembre, circa 37 minuti (il grosso è Newton sui tessuti 3D, a matrice piena o sparsa); `scripts/prove.sh veloci` salta i 21 gruppi
marcati `@lenta` e dura circa 2 minuti (780 prove). Quando lanciare quale: `CLAUDE.md`. Confronti esterni con tyssue in `confronto/` (`LEGGIMI.md`).

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

### R3 — Chimica (specifica, decisione 5) (fatto)

Stato: **fatto e congelato** (2 novembre, con `adatta_alla_scatola!` venuta da R4). Commit verde
del 31 ottobre. Le dieci prove di `test/chimica.jl`, scritte prima,
passano; il codice è in `motore/12_chimica.jl`: griglia cartesiana 3D a voxel (passo per asse, bordi chiusi o
periodici), schema esplicito con sotto-passi scelti per la stabilità, degradazione e secrezione
integrate esattamente nel passo (così il bilancio stazionario non ha errore di passo), secrezione e assorbimento alla posizione con pesi trilineari (assorbimento
esponenziale, così la concentrazione non diventa negativa), `centro_cellula` per il guscio.
Parametri di riferimento dalla letteratura: Dpp nel disco dell'ala, D = 0.1 µm²/s = 6 µm²/min e
λ ≈ 20 µm (Kicheva et al., Science 315, 521, 2007), quindi k = D/λ² = 0.015/min.

Numeri (`docs/piano.md`, sezione 14): λ misurata con errore 2·10⁻⁴ a h = 1 µm e 8·10⁻⁴ a
h = 2 µm (secondo ordine); varianza 6Dt esatta; sorgente puntiforme entro il 3% a 5–8 µm;
errore nel tempo del primo ordine (rapporto 2.07 dimezzando il passo).

La griglia segue la scatola dei fogli periodici (`adatta_alla_scatola!`, prova in
`test/programma.jl`): i passi x e y diventano Lx/nx e Ly/ny e le concentrazioni si scalano con
il volume del voxel, così la massa non cambia. La passata lo fa prima della chimica. Resta la
regressione del #1 (necrosi da limite di diffusione), che chiede la morte delle cellule.


Almeno una sostanza che diffonde, con sorgenti, consumo e degradazione, risolta come PDE su
griglia con passo proprio, accoppiata alle posizioni delle cellule (secrezione e lettura alla
propria posizione). Prove analitiche prima: profilo di una sorgente puntiforme, stato stazionario
esponenziale con lunghezza √(D/k), conservazione della massa senza degradazione, convergenza con
la griglia e con il passo. Coefficienti di diffusione e degradazione da letteratura, con la fonte.
Regressione: il risultato negativo #1 (la necrosi da limite di diffusione non stabilizza un lume)
va riprovato con la chimica vera.

### R4 — Percezione e programma (specifica, sezioni 5–8) (fatto)

Stato: **commit verde fatto** (2 novembre). Le undici prove di `test/programma.jl`, scritte
prima, passano. Il codice è in `motore/13_programma.jl` (percezione, fasi, passata) e
`primitive/costringi.jl` (la prima primitiva); `adatta_alla_scatola!` in `motore/12_chimica.jl`,
ora congelato. Le prove: la percezione non può portare lo stato globale (si percorrono i tipi); la regola riceve
solo la percezione; l'onda delle fasi è monotona e avanza di un contatto per passata, e spenta
lascia avanzare solo chi ha visto la fine da sé; l'orologio è della cellula (burattinaio 3);
il cono di luce, cioè nessuna decisione globale (burattinai 1 e 2); l'ordine delle cellule non
conta; gli interruttori; `costringi` fa i cunei dove la sostanza secreta dalla cellula centrale
supera la soglia; la marca è un'ipotesi sul presente (spenta la sorgente, le cellule smettono
da sé); la griglia segue la scatola; le corse sono riproducibili.

Il numero di `costringi` (regressione, `docs/piano.md` sezione 15): sul pezzo di 37 cellule a
bordo bloccato, la cellula centrale secerne, le 7 dentro la soglia tendono l'anello (Λ = 3 nN),
e dopo 30 passate da 1 minuto il rapporto fra area apicale e basale è 0.31–0.45 dentro e 1.43
nel secondo anello (le vicine si allargano all'apice).

**Chiusura di R4: fatta** (2 novembre), `13` e `14` congelati. Commit rosso: Otto prove in `test/popolazione.jl`
(divisione che conserva il volume e gli invarianti, divisione secondo Hertwig, la figlia che
eredita fase, stato e programma, popolazione ogni m passate e il suo interruttore, morte come
estrusione, matrice che persiste, casualità con seme) e una in `test/operazioni.jl`: un difetto
trovato in `11` congelato (sotto). Il codice va in `13` e in un file nuovo,
`motore/14_popolazione.jl`. Commit verde: tutte passano.

Nella passata ora ci sono tutti i passi della specifica: vicinato, chimica, percezione e
programma, meccanica, topologia, popolazione (ogni `ogni` passate), matrice, registro. Dettagli
che le prove hanno chiesto (`docs/piano.md`, sezione 16):

- la cellula in apoptosi ha il programma fermo e il volume bersaglio che scende con un tempo
  caratteristico (`tempo_apoptosi`, 5 min, libero), non di colpo: un crollo da 1000 a 1 µm³ in
  un passo non lascia convergere il passo implicito;
- nel guscio a prismi la topologia fa uno scambio di vicine per passata, il più corto: più
  scambi insieme in una cellula minuscola, con pareti nuove più lunghe della cellula, la
  rovesciano;
- se il passo implicito non converge (la topologia è in ritardo) la passata va avanti e lo
  registra fra gli eventi, come fanno le funzioni di rilassamento validate;
- se il piano di divisione passa per un vertice (la regola di Hertwig su un esagono regolare
  allungato lo fa) lo si ruota di 0.01 rad.

Scelta: divisione e morte solo nel guscio a prismi (scambi di `08`, estrusione di `09`), che le
sa già fare. Nel modello generale una cellula che si restringe tutta insieme arriva al prisma
triangolare, ma passando per un volume negativo, e lì si ferma perché nel modello generale
non c'è un'estrusione: risultato negativo #16. La morte nel modello generale è un lavoro di R7.

I difetti di `11`, corretti con le prove che fallivano e l'impronta nuova:
`spigolo_a_triangolo!` non controllava che le facce dello spigolo, che perdono un vertice,
restassero almeno triangoli (`prova_spigolo_di_una_faccia_triangolare_si_rifiuta`), e lo stesso
`triangolo_a_spigolo!` per le facce con due vertici del triangolo (trovato dall'esperimento del
#16, `prova_risultato_negativo_16`). Rompevano gli invarianti e sollevavano un errore invece di
rifiutare. Il #14 come prova dinamica (l'anello che stringe
fino a uno scambio) si fa con R5.


La struttura di percezione con i soli quattro canali, costruita senza accesso allo stato globale
(il tipo lo impedisce); la tabella di fasi copiata in ogni cellula; l'avanzamento monotono
(evento di fine rilevato da sé, oppure una vicina già nella fase dopo); le intenzioni meccaniche
come uscita di una fase (tensioni, volume bersaglio, richieste di dividersi, morire,
differenziare). Prima primitiva: `costringi` come anello apicale. I tre burattinai del codice
vecchio (sezione 8 della specifica) diventano prove che non devono ripresentarsi.

### R5 — Prime forme (in corso)

Stato: **prima forma fatta** (2 novembre): l'invaginazione di una placca, in red-green. Sei
prove in `test/invaginazione.jl`, tutte verdi. Il codice è in `programmi/composizione.jl`
(`insieme`, il ‖ dell'algebra, e `con_fine`), `programmi/invaginazione.jl` (il programma) e
`primitive/secerni.jl`; il motore congelato non è stato toccato. Le prove sono comportamenti
(`prove/attesi.md`, quaderno): la placca che costringe scende dal lato basale, senza
costrizione resta piatta, più forte la costrizione più profonda, mantenuta resta, rilasciata si
riapre (la piega è elastica); `insieme` rifiuta due regole che toccano la stessa intenzione
(la legge del conflitto manca ancora, specifica sezione 5). Un esperimento prima delle prove:
nel foglio periodico 8 × 8, sette cellule con Λ = 3 nN abbassano l'apice di 4.8 µm.

Numeri (regressione, `docs/piano.md` sezione 17): foglio periodico 8 × 8 all'equilibrio,
placca di 7 cellule, passo meccanico 10 min. Con Λ = 3 nN la placca scende di 4.26 µm in 12
passate e va all'equilibrio a 4.47 µm; con Λ = 1.5 nN di 2.34 µm; rilasciata, torna a 0.27 µm
in 12 passate (la piega si muove con un tempo di circa 35 minuti). Tre scelte delle prove
corrette nel verde, dichiarate in testa al file.

L'organizzatore è per ora una condizione iniziale; l'elezione che lo sceglie senza seme viene
con la catena di gastrulazione.

**Il conflitto alla composizione: fatto** (2 novembre). `insieme` riconosceva un
conflitto confrontando i valori delle intenzioni con quelli a riposo: una regola che scrive il
valore a riposo (`costringi(Λ = 0)`, che spegne) non si distingueva da una che non tocca niente,
e `costringi(Λ = 0) ‖ costringi(Λ = 3)` dava Λ = 3 in silenzio. Correzione: ogni fase dichiara
che cosa scrive, `insieme` rifiuta i conflitti quando compone e la passata rifiuta una regola
che scrive fuori dalla dichiarazione. Quattro prove in `test/composizione.jl`, scritte prima;
la correzione tocca `13` congelato (il campo `scrive` di `Fase` e il controllo nella passata),
con l'impronta nuova e la ragione in `prove/congelati.txt`. `costringi` dichiara `(:Λa,)`,
`secerni` `(:secrezioni,)`; una fase senza dichiarazione può scrivere tutto ma non si compone.
Lo stato interno e le dichiarazioni alle vicine restano controllati sui valori, alla corsa.

**Uscite a delta e dichiarazioni di stato e vicine: fatto** (2 novembre). Lo
stato interno e le dichiarazioni alle vicine si univano ancora per valori, e ogni regola
restituiva lo stato intero: un contatore composto con `costringi` dava un conflitto falso, e
una scrittura del valore a riposo restava invisibile. Analisi delle alternative (lock, MVCC,
proxy) in `docs/piano.md`, sezione 18. Scelta: dichiarazioni statiche anche per le voci di
stato e dichiarazione (il conflitto alla composizione, per ogni cellula) più uscite a delta,
che contengono solo ciò che la regola scrive (il registro esatto delle scritture, che la
passata confronta con la dichiarazione). Cinque prove in `test/composizione.jl`, scritte prima;
tocca `13` congelato, con l'impronta nuova. Una regola ora restituisce solo ciò che scrive:
`Uscita(; intenzioni = (Λa = 3.0,), stato = Dict(:n => 1.0))`; lo stato non scritto resta, le
intenzioni non scritte tornano a riposo, la dichiarazione è il messaggio della passata.

**L'organizzatore senza seme: fatto** (3 novembre), `primitive/attiva_inibisci.jl`. Il criterio di ammissione è
ora più stretto (`CLAUDE.md`, regola 3): ogni primitiva ha un meccanismo cellulare con una fonte,
e `test/primitive.jl` lo controlla. L'elezione a biglietti del quaderno è un protocollo di
calcolo senza meccanismo noto: la sostituisce l'attivazione locale e inibizione a lungo raggio
(Gierer e Meinhardt; Nodal/Lefty). Con i parametri misurati (Müller et al. 2012: Cyclops e
Lefty2) c'è un regime di Turing solo con una saturazione dell'autoattivazione (κ ≈ 0.29), e
l'instabilità ha lunghezza d'onda di circa 1.1 mm e tempo di crescita di circa 13 ore: una sola
macchia nasce senza seme in un tessuto di circa 1 mm, in giorni (come negli aggregati di Hydra e
nei gastruloidi; nel pesce zebra la posizione la fissa invece un pre-pattern materno). Scelta
dell'utente (3 novembre): la primitiva provata a scala vera, sola chimica (`test/organizzatore.jl`,
tre prove `@lenta`, scritte prima, verdi); per la catena sul foglio piccolo, un pre-pattern
dichiarato. La catena completa a scala vera chiede il solutore sparso (R6). Numeri (foglio
periodico di 1672 cellule, circa 1 × 1 mm, passata di 30 min): con il seme 1 il contrasto
dell'attivatore è 1.18 dopo un giorno, 8 dopo tre, 57 dopo quattro, una sola macchia; con il
seme 2 una sola macchia altrove, più lenta (2.1 a quattro giorni); con le diffusività uguali
nessuna macchia (1.003). Costo: 23 minuti per le tre corse, quasi tutto nella secrezione dentro
`avanza!` di `12` (pesi della griglia con tipo instabile). **Rifatto il 3 novembre** con κ = 0.35
e passate di 5 minuti (risultato negativo #17: con 30 minuti il circuito oscillava per un difetto
numerico; la passata deve stare sotto 2δ/ω²): una sola macchia verso il quinto o sesto giorno,
contrasto 17 (seme 1) e 14 (seme 2) a 7 giorni, in posti diversi; a diffusività uguali niente
(1.04). Decisione dell'utente (3 novembre):
eccezione dichiarata per le prestazioni (`CLAUDE.md`, regola 1). Fatto in red-green:
`test/prestazioni.jl`, una prova d'identità (l'impronta delle concentrazioni di un caso fisso con
1672 sorgenti, presa prima) e una di tempo (meno di 0.5 s). La correzione è una riga in `12`
(il tipo dei pesi scritto per esteso): da 2.4 a 0.06 s, stessa impronta; la batteria intera,
con l'organizzatore, scende a 8 minuti e mezzo.

**Il gene bersaglio: fatto** (3 novembre), `primitive/induci.jl`, integrato in modo esatto nella
passata (quindi senza il problema del passo esplicito del #17). A regime la macchia di Nodal arriva a
circa 1.75 volte il livello di base, e una soglia sulla concentrazione istantanea accende e spegne
le cellule al bordo. La forma plausibile è un gene bersaglio che si accumula con il segnale e
decade senza (Dubrulle et al., eLife 2015: contano la cinetica e la durata dell'induzione, non
una soglia sulla concentrazione). Tre prove in `test/induci.jl`, scritte prima, verdi; `costringi`
può leggere il prodotto del gene (`voce`) invece della concentrazione.

**La catena a scala vera, primo tentativo** (3–4 novembre): `attiva_inibisci ‖ induci`, con
l'impegno quando il prodotto del gene supera metà del massimo, poi l'onda e `costringi` sul
prodotto del gene. Alla 132ª ora simulata l'organizzatore si è impegnato da solo, l'onda ha
portato tutte le 1672 cellule nella seconda fase, 78 cellule costringono e la regione è scesa di
3.6 µm. Fermata lì: la meccanica costava circa 15 s per passata (5–10 iterazioni di Newton, ognuna
con una matrice nuova da 2.7 s; la fattorizzazione ne costa 0.3).

I rimedi (4 novembre), con le prove scritte prima:
- `meccanica_ogni` in `13` (specifica, decisione 6: la meccanica ogni m passate, con passo m·dt;
  estensione compatibile): **fatto**, `prova_meccanica_ogni` verde, impronta predefinita
  identica. La chimica di Nodal/Lefty chiede passate di 5 minuti, la meccanica no (il suo tempo
  è circa 35 minuti): con m = 6 i passi meccanici sono sei volte meno.
- Newton modificato (riusare la matrice fra le iterazioni): **provato e scartato**. In questo
  regime la forza non si dimezza con la matrice vecchia, che va ricalcolata quasi sempre, e i
  tentativi falliti costano in più (passo da 12 a 17 s).
- Matrice dai gradienti parziali (`15`): per ogni gruppo di colore il gradiente si calcola solo
  sulle cellule che contengono i suoi vertici, un sotto-tessuto. **Fatto**: la matrice scende da
  2.7 a 1.75 s, il passo da 12 a circa 9 s.
- Matrice in parallelo sui fili (`15`, `scripts/prove.sh` con `--threads=auto`): **fatto**,
  identica bit per bit con 1 e con 4 fili, ma quasi senza guadagno. Ogni matrice alloca 4 GiB in
  8 milioni di pezzi (il gradiente di `07` crea un vettore nuovo per ogni faccia), e il garbage
  collector ferma tutti i fili insieme.
- Obiettivo di `prova_passo_grande_veloce` (meno di 5 s per passo): non raggiunto con i soli
  rimedi di `15` (circa 9 s); raggiunto dopo, con il gradiente senza allocazioni (sotto).
- Perché la batteria intera è passata da 13 a 26 minuti (misurato il 4 novembre): non per i fili.
  Con 4 fili i file pesanti vanno uguali o più veloci che con uno (`invaginazione` 248 s contro
  305); ma con un filo solo file il cui codice non è cambiato (`scambi`, `estrusione`) sono 1.4–1.5
  volte più lenti della misura del 3 novembre: la macchina, oggi, va più piano. Il resto sono le
  prove lente nuove (organizzatore, passo grande, scala vera).
- **Gradiente di `07` senza allocazioni: fatto**, con l'eccezione per le prestazioni: un vettore
  di lavoro riusato al posto di uno nuovo per faccia, in gradiente, energia, volumi e controllo.
  L'impronta di gradiente, energia e volumi su tre tessuti (bloccato, periodico, con uno scutoide)
  è identica bit per bit. Le allocazioni del gradiente sul foglio grande passano da 7.9 MB a 208
  byte; il suo tempo solo da 9.95 a 8.4 ms (il costo è calcolo), quindi l'obiettivo dei 4 ms resta
  rosso (`prova_gradiente_veloce`, `@test_broken`).
- Trovata e tolta anche una perdita mia in `15`: i pezzi della matrice si univano con un
  `reduce` che ricopiava tutto a ogni aggiunta (3.4 GiB per matrice).
- Con le due correzioni i fili lavorano: la matrice passa da 2.7 s (all'inizio) a 1.2 s con un
  filo e 0.65 s con quattro, identica bit per bit; **il passo meccanico grande scende a 4.2 s**, e
  `prova_passo_grande_veloce` (meno di 5 s) è verde. Riusata anche l'analisi simbolica della
  fattorizzazione (l'ordinamento non cambia fra le iterazioni): dopo le prove pesanti il passo
  costa 3.8–3.9 s invece di 4.4–4.75, e il foglio grande con la regione larga 300 µm va
  all'equilibrio in 84 iterazioni e 111 s (erano 237). Un errore trovato per strada: dopo una
  fattorizzazione fallita l'oggetto riusato restava nello stato di errore e μ cresceva senza fine
  (`prova_scala_vera` non convergeva più); ora dopo un fallimento si ricomincia da capo.
- Corretta, nel frattempo, una prova mal posta: `prova_stesso_minimo` confrontava le posizioni
  assolute di un foglio periodico senza ancora, che si può traslare rigidamente senza cambiare
  energia; ora le confronta a meno della traslazione (la forma coincide a 10⁻¹²).

**La catena a scala vera: fatta** (4 novembre), `programmi/gastrula.jl`. Tre prove `@lenta` in
`test/gastrula.jl`, scritte prima, verdi: nei primi quattro giorni nessuna cellula si impegna
(l'organizzatore nasce in giorni); poi una sola regione si impegna (fra 20 e 600 cellule,
connessa); le cellule impegnate tendono l'anello e la regione scende di oltre 3 µm dal lato
basale. Il programma è organizza ; invagina, con `meccanica_ogni` = 6 e il solutore sparso.
Numeri (seme 1): la prima cellula si impegna alla 128ª ora (5.3 giorni), si impegnano 331 cellule
in una sola regione, la regione scende di 17.9 µm, nessun passo meccanico irrisolto. La prova
dura circa 19 minuti (14.7 con la matrice esatta, sotto).

**Invarianza di scala: prove scritte** (5 novembre), `test/scala.jl`, prima della corsa: la
catena a scala vera con cellule di lato 12 µm invece di 15 sullo stesso millimetro (2688 cellule),
il foglio sempre alto 15 µm; in unità di cellula scalano il volume a riposo, σl e Λ. Chiedono la
stessa area impegnata (entro il 25%), lo stesso tempo (entro il 15%), la stessa profondità (entro
il 30%).

Prossimi passi di R5: la stessa forma a due taglie di cellula (invarianza di scala,
`prove/attesi.md`); il #14 come prova dinamica (l'anello che stringe fino a uno scambio);
la placca più larga che si chiude in una coppa (servirà il contatto in 3D, R6); poi la catena di
gastrulazione con elezione, campo, costringi e differenzia (specifica, sezione 5). La batteria
dura ormai circa 10 minuti: separare prove veloci e lente (R6) diventa urgente.


L'invaginazione 3D di una placca che stringe l'anello apicale (oggi possibile: l'apice si
stringe e il foglio si piega), poi la catena di gastrulazione compilata (specifica, sezione 5).
Confronto con `prove/attesi.md` sui rapporti e sui comportamenti, non sui valori.

### R6 — Infrastruttura (quando una forma la chiede)

**Il solutore sparso: fatto** (3 novembre), `motore/15_sparso.jl`. Lo chiede la catena di gastrulazione a
scala vera: la catena sul foglio piccolo con un pre-pattern non è plausibile nei numeri (Nodal
ha una portata di circa 76 µm e Lefty di circa 450, il foglio piccolo misura 60 µm: il segnale lo
riempirebbe tutto), quindi si va a circa 1 mm e 1700 cellule, dove Newton a matrice piena non
regge (circa 20 000 incognite). Il codice va in `motore/15_sparso.jl` (`minimizza_sparso!`:
lo stesso metodo di `minimizza!`, con la matrice delle derivate seconde sparsa, calcolata per
differenze del gradiente esatto con una colorazione a distanza due, e una fattorizzazione
sparsa); la passata guadagna l'opzione `solutore = :pieno o :sparso` con la seconda eccezione
della regola 1, l'estensione compatibile (decisione dell'utente, 3 novembre). Prove:
`test/sparso.jl` (quattro, marcate `@test_broken`: matrice uguale alla piena, stesso minimo,
passata sparsa come la piena, un foglio di 1672 cellule all'equilibrio) e `test/estensioni.jl`
(l'impronta della passata predefinita, presa prima). Tutte verdi. Numeri: sul foglio di 1672
cellule (6688 vertici, circa 20 000 incognite) bastano 44 colori, cioè 264 valutazioni del
gradiente per matrice invece di circa 40 000; con una regione larga 300 µm che costringe
(Λ = 8 nN) l'equilibrio arriva in circa 4 minuti. Corretto nel verde, nella prova: il foglio
deve partire all'equilibrio (σl = 2√3 per prismi di lato e altezza 15 µm, tensioni libere).
`minimizza_sparso!` resta fuori dal congelamento finché la catena a scala vera non l'ha usato.

- **Solutore sparso**: Newton con matrice piena non regge migliaia di cellule; servono la matrice
  delle derivate seconde sparsa e un solutore adatto.
- Prove veloci e lente separate: fatto il 3 novembre (`test/lente.jl`, `scripts/prove.sh veloci`).
- **Contatto in 3D**: oggi solo nel modello laterale; senza, un foglio che si piega può
  attraversarsi.
- Le operazioni topologiche nei fogli periodici: fatto con R1.
- Nel modello generale: l'estrusione (la T2 dei modelli a vertici 3D) e pareti nuove
  proporzionate alla cellula (risultato negativo #16), poi divisione e morte senza il guscio a
  prismi.

**La matrice delle derivate seconde esatta: fatto** (5 novembre), `motore/16_hessiana.jl`. Il
collo di bottiglia della meccanica a scala vera era la matrice, calcolata per differenze del
gradiente (264 gradienti). Il gradiente di ogni termine di faccia (σA con l'anello apicale, e il
contributo di volume) è scritto in forma chiusa, sulle sole posizioni della faccia; la sua
matrice jacobiana, cioè le derivate seconde, viene dalla differenziazione automatica in avanti
(ForwardDiff), esatta a precisione di macchina, senza scrivere a mano le derivate seconde del
baricentro del contorno. Ogni cellula raccoglie in una matrice piena il suo volume
(K·(∇V ∇Vᵀ + (V − V0) ∇²V)) e le facce di cui è la prima cellula. `15` e `13` hanno guadagnato
un'opzione (estensioni compatibili: `minimizza_sparso!(...; hessiana = :esatta)`,
`Embrione(...; solutore = :esatto)`), con l'impronta del predefinito presa prima
(`test/estensioni.jl`). Prove in `test/hessiana.jl`, tutte verdi: uguale alla matrice per
differenze e simmetrica, stesso minimo, la passata come quella sparsa, almeno due volte più
veloce sul foglio grande. Numeri (foglio di 1672 cellule, 4 fili): la matrice da circa 0.7 a 0.3
s, il passo meccanico grande da 4.4 a 2.9 s, con le stesse iterazioni. Scartato: differenziare
due volte l'energia (ForwardDiff.hessian), giusto ma 6 s per matrice, dieci volte la versione
per differenze. Il resto del passo è la fattorizzazione e le energie della ricerca del passo.
La catena a scala vera (`test/gastrula.jl`) ora usa `solutore = :esatto`: stessi numeri (prima
cellula impegnata alla 128ª ora, 331 cellule, −17.9 µm, nessun passo irrisolto) in 14.7 minuti
invece di circa 19. Congelato il 5 novembre (`CONGELATO: MOTORE-HESSIANA-V1`).

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
