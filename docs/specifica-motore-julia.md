# Motore morfogenetico: specifica per la riscrittura in Julia

Oct 27, 2026 · @Andrea

## 1. Perché si riscrive

Il motore attuale (JavaScript, cellule sferiche, `motore3d.js`) non si riusa. Non perché sia scritto male, ma perché è stato costruito al contrario: ogni forma nuova faceva emergere un difetto, e il difetto veniva corretto nel motore. Dopo una decina di forme il motore è un **sedimento di correzioni**, non un modello. Alcune correzioni erano giuste, altre erano aggiramenti di difetti di rappresentazione, e dal codice non si distinguono.

Tre esempi concreti, tutti misurati:

- **La cellula non può essere un cuneo.** Una sfera con una sola lunghezza di riposo per contatto non può diventare stretta all'apice e larga alla base. La costrizione apicale fra sfere accorcia il foglio senza curvarlo. L'aggiramento è stato una parete di due file che si comporta come una striscia bimetallica: funziona, ma è un trucco imposto da una rappresentazione sbagliata, e in 3D una calotta sferica (curva in due direzioni) non si piega con gli stessi parametri che bastavano a un arco piano.
- **Non c'è un'energia.** Le forze sono inventate e tarate a mano, quindi vanno appaiate a mano. `costringi` è stato implementato come accorciamento della *coppia* proprio per garantire l'appaiamento. Gli equilibri non sono definiti e i parametri non hanno significato fisico.
- **Non c'è chimica.** Il limite di diffusione per il lume è stato finto con una distanza su grafo, e si è rivelato autodistruttivo: la regola cancella la propria premessa (chi muore rende il vicino meno profondo). Non era sfortuna, era chimica finta.

A questo si aggiunge un difetto di metodo che va corretto per primo: **il motore non deve essere modificabile da una forma.** Dopo la validazione, una modifica al motore richiede un test *fisico* che falla, non una forma che si comporta male.

Il lavoro fatto non si butta: diventa **dato di riferimento** — le misure, le forme che funzionano, e soprattutto i risultati negativi. Non fondamento.

## 2. Le sei decisioni di fondo

Queste sono prese in prestito dai motori maturi (PhysiCell, Chaste, CompuCell3D, Morpheus, tyssue, EmbryoMaker). Non ci servono le loro librerie: ci servono le loro decisioni. Ognuna è legata al fallimento misurato che la motiva.

| decisione | cosa significa | il fallimento che la motiva |
| --- | --- | --- |
| **1. Un funzionale di energia** | Ogni interazione meccanica discende da un'energia; la dinamica la scende in regime sovrasmorzato. Le forze sono gradienti, quindi appaiate per costruzione. | `costringi` appaiato a mano; costanti di molla ritarate a ogni forma; il "corpo plastico" che irrigidì tutto e fu rimosso |
| **2. Geometria che possa essere un cuneo** | Vertici e facce; tensione di superficie distinta per faccia apicale, basale e laterale; tensione di linea sugli spigoli. | nessuna piega in un monostrato; la striscia bimetallica come aggiramento; la calotta 3D che non si invagina |
| **3. Unità vere** | Micrometri, minuti, nanonewton, pascal. Ogni parametro o viene da una misura in letteratura, o è dichiarato libero. | "unità di forma" e "tick": un modello che non può essere smentito da un esperimento non può nemmeno essere confermato |
| **4. Cambio topologico controllato** | Scambio di vicine e rimozione di cellula come operazioni discrete, esplicite, con le loro condizioni e registrate. | la lacerazione quando le due labbra del solco si costrinsero l'una contro l'altra (pezzi: 2); il tessuto deve *scorrere*, non strapparsi |
| **5. Chimica vera, anche minima** | Almeno una sostanza che diffonde, con sorgenti, consumo e degradazione, risolta come PDE su griglia, con passo proprio. | il limite di diffusione finto con una distanza su grafo: 6 prove a tre soglie, nessuna quiete, nessuna cavità |
| **6. Scale di tempo separate** | Chimica veloce, rilassamento meccanico intermedio, divisione e morte lente; passi diversi e rapporto dichiarato. | un tick per tutto; parte dei problemi di stallo e dei cicli futili vengono da lì |

Una conseguenza da accettare subito: **le decisioni 2, 3 e 5 rendono impossibile riusare il codice attuale.** È una riscrittura, non un refactoring.

## 3. Geometria cellulare

**Vertici e facce.** I gradi di libertà sono i vertici; ogni cellula è un poliedro con una lista di facce, ciascuna etichettata apicale, basale o laterale. L'energia ha tre termini: elasticità di volume (incomprimibilità), tensione di superficie per faccia con valori distinti per etichetta, tensione di linea sugli spigoli.

In questa forma `costringi` non è più una forza inventata né un trucco a due file: è **alzare la tensione dell'anello apicale di una cellula**. Minimizzare quell'energia a volume costante *è* diventare un cuneo. Il difetto che ci ha bloccato sparisce per rappresentazione, non per taratura.

Sparisce anche un bug che era un effetto collaterale di avere solo distanze fra centri: la tensione apicale è proprietà della **faccia apicale**, non della coppia di cellule, quindi due cellule che si toccano fuori dalle rispettive facce apicali non hanno nessun motivo di accorciarsi. È la lacerazione del solco, risolta alla radice.

**Il 3D come guscio a due superfici.** Restrizione di scopo deliberata, per arrivare in fondo: una maglia triangolata apicale e la sua copia basale; ogni cellula è il prisma fra un poligono apicale e il corrispondente basale, con le facce laterali che li collegano. Si ottengono cunei naturali, l'anello apicale naturale, un guscio chiuso per la blastula, e — decisivo — i cambiamenti topologici restano confinati alla superficie, quindi sono scambi di vicine bidimensionali, trattabili. Un modello a vertici 3D *pieno*, con cellule che si riarrangiano in volume, è il punto dove questi progetti si impantanano, e non ci serve: tutto ciò che vogliamo piegare è un foglio.

**Due rappresentazioni, una algebra tipata.** Il prezzo della restrizione è che il mesenchima non c'è: una palla piena di cellule, il lume che si riempie, una cellula che migra da sola non sono prismi di una maglia. Quindi:

- cellula **epiteliale** = prisma in una maglia a due superfici;
- cellula **mesenchimale** = centro con il proprio volume, vicinato per triangolazione.

Le primitive si scrivono una volta e scelgono l'implementazione sul tipo — ed è qui che la scelta su *tutti* gli argomenti si paga da sé: `aderisci` fra epiteliale ed epiteliale, fra epiteliale e mesenchimale, fra mesenchimale e matrice sono tre casi distinti; e `costringi` su un centro **non esiste**, il tipo lo rifiuta.

Questo non è solo ordine: l'algebra diventa **tipata**, e il tipo dice quali forme sono raggiungibili con quali cellule. È un enunciato nuovo da aggiungere al quaderno.

## 4. Il vicinato

Per le cellule a centro, "chi tocca chi" si definisce con la **triangolazione di Delaunay** (il duale della costruzione di Wigner–Seitz), non con un raggio di taglio.

Motivo misurato: nel motore attuale il vicinato è `1.25 × la lunghezza di riposo`, scelto a mano. Allargandolo a 1.35 per far funzionare una costrizione anisotropa, ho degradato silenziosamente tutto il resto e ho dovuto tornare indietro. Un parametro arbitrario nel posto peggiore — quello da cui dipendono tutti i canali di informazione.

Per la stessa ragione **Voronoi non si usa per la forma**: le celle di Voronoi sono convesse, e una cellula che si invagina davvero è a fiasco, cioè non convessa. Quella forma è il cuore del fenomeno e non si può perdere. Wigner–Seitz serve per il vicinato, i vertici per la forma.

## 5. Come un programma si compila nelle primitive della cellula

Un programma non viene eseguito *sopra* le cellule: viene **copiato dentro ognuna**, identico, come un genoma. Nessuna cellula è speciale, nessuna sa di far parte di una sfera o di un toro.

La compilazione di una catena `p₁ ; p₂ ; … ; pₙ` produce una **tabella di fasi**. Ogni cellula porta una copia della tabella più l'indice della fase che sta eseguendo. Per ogni fase la tabella contiene tre cose, e nient'altro:

1. la **regola**, una funzione pura della sola percezione;
2. i **parametri**;
3. il **predicato di fine**, anch'esso locale.

**L'avanzamento di fase non ha orologio centrale.** Una cellula passa da *k* a *k+1* quando una delle due è vera: ha rilevato da sé l'evento di fine di *k*, oppure una vicina è già in *k+1*. La seconda clausola è l'onda: la fase si propaga per contatto. L'avanzamento è **monotono**, non torna indietro mai — è questo che rende una sequenza una sequenza senza un direttore d'orchestra.

Nell'algebra: `;` è concatenazione di fasi; `sempre mantieni(p)` è una fase il cui predicato di fine è sempre falso; `‖` è due fasi attive insieme, ed è esattamente il caso di cui ci manca la legge quantitativa.

### La gastrulazione compilata

| fase | regola (solo percezione) | evento di fine |
| --- | --- | --- |
| 0 cresci | fascia \[Rin, R\]: divido dove la figlia resterebbe dentro, muoio se ne sono fuori | sono quieta e nessuna vicina ha avuto eventi per T passi |
| 1 eleggi | biglietto congelato = il mio ρ; propago il massimo con un'età che scade; mi astengo se una vicina è già chiusa o se non ho vicine | il mio biglietto è imbattuto da 45 passi → **sono io il polo** |
| 2 campo | d = minimo fra le vicine di (d + distanza); il polo tiene 0 | il mio d cambia meno di ε per 60 passi |
| 3 costringi | se d < a e la mia faccia apicale guarda via dal centro: alzo la tensione dell'anello apicale; arresto mitotico in me e entro due contatti; **a cresce sul mio orologio personale** da quando sono entrata in questa fase | a ha raggiunto a\_max |
| 4 differenzia | se stavo costringendo divento tipo B: post-mitotica, immune alla fascia madre, tensione apicale a riposo | nessuno (terminale) |

Due vincoli emersi dalle prove, da portare nella riscrittura:

- **Le fasi in sequenza condividono parametri.** `a_max` è contemporaneamente l'evento di fine della fase 3 e la taglia che servirebbe al dominio della fase 4. La composizione non è libera: un dominio **eredita la propria taglia dal processo che l'ha fatto**, non la riceve come dato.
- **Una fascia non si applica prima che il suo campo sia converguto.** Alla prima iterazione tutte le distanze sono ignote, e "ignota" non è un motivo per morire. Nella prima versione il dominio differenziato si estingueva in 200 passi; lo stesso errore era costato il campo di profondità. Vale come regola generale: *una distanza ancora in discesa non autorizza nessuna decisione irreversibile.*

## 6. Una passata del motore

Una passata è **un passo meccanico** (ordine della decina di secondi di tempo simulato). Dentro ci stanno più passi di chimica, perché la diffusione è veloce; la popolazione si aggiorna solo ogni *m* passate, perché dividersi e morire sono lenti.

1. **Vicinato.** Si ricostruiscono i contatti (triangolazione per i centri, adiacenza di maglia per gli epiteliali). I contatti definiscono "le vicine" per tutti i canali.
2. **Chimica**, *n* sotto-passi: diffusione e degradazione di ogni sostanza sulla griglia; secrezione o consumo di ogni cellula alla propria posizione. È l'unico punto in cui le concentrazioni cambiano.
3. **Percezione.** Ogni cellula assembla i propri ingressi, e **solo dai quattro canali**: il proprio stato interno; ciò che le vicine dichiarano, letto nel proprio riferimento; la matrice alla propria posizione; le misure assolute. La struttura di percezione si costruisce qui, e **non ha accesso allo stato globale**. Se una regola ha bisogno di qualcosa che non c'è, non si aggiunge una scorciatoia: si dichiara che quella regola non è ammissibile.
4. **Passo di programma.** Ogni cellula esegue la regola della propria fase sulla percezione e produce: il nuovo stato interno; l'avanzamento di fase o no; le **intenzioni meccaniche** (tensioni per faccia, tensione di linea apicale, coefficienti di adesione per tipo di partner, asse di polarità, volume bersaglio); e le **richieste** di dividersi con un piano di divisione, morire, differenziare. Tutto calcolato sullo stato della passata precedente, quindi l'ordine delle cellule non conta.
5. **Meccanica.** Si costruisce l'energia dalla geometria attuale e dalle intenzioni: elasticità di volume per cellula, tensione di superficie per faccia, tensione di linea sugli spigoli apicali, adesione come energia interfacciale negativa dipendente dalla coppia di tipi, ancoraggio alla matrice. Poi si scende l'energia di un passo, in regime sovrasmorzato. Le forze sono gradienti: appaiate per costruzione.
6. **Topologia.** Si applicano i cambiamenti discreti e dichiarati: scambio di vicine dove una giunzione si è accorciata sotto soglia, rimozione di una cellula collassata. Registrati, non impliciti.
7. **Popolazione**, solo ogni *m* passate: si eseguono le richieste accolte al passo 4. La figlia nasce con la propria geometria, la madre si dimezza, la figlia **erede fase, identità e programma**; le morte liberano il loro posto e la loro matrice.
8. **Matrice.** Ogni cellula depone sul proprio lato basale; il deposto persiste; chi muore lo libera.
9. **Registro.** Si annotano gli eventi — serviranno alle vicine come percezione nella passata successiva — e si calcolano le misure per l'osservatore. Queste ultime possono essere globali quanto si vuole: **nessuna cellula le legge.**

Due idealizzazioni dichiarate, da non nascondere:

- l'aggiornamento è **sincrono**, mentre le cellule vere sono asincrone. Lo facciamo per riproducibilità.
- l'orologio interno di una cellula è un contatore di passate: assumiamo che sappia contare il tempo, che è vero (i cicli cellulari sono orologi) ma con una precisione che una cellula vera non ha.

## 7. Repertorio chiuso dei comportamenti

Il motore offre questi comportamenti e **nessun altro**. Una primitiva che ne richiede uno fuori elenco non è ammissibile, e la risposta non è aggiungerlo al motore ma ridisegnare la primitiva.

- **dividersi**, con un piano di divisione (planare nel foglio, o orientato)
- **morire**: apoptosi, oppure estrusione dal foglio
- **aderire**, con intensità dipendente dalla coppia di tipi e dalla faccia
- **polarizzarsi**: apico-basale, e planare (un asse nel piano del foglio)
- **cambiare forma**: accorciare una faccia, cambiare volume bersaglio
- **muoversi e intercalarsi** (per i mesenchimali; per gli epiteliali è lo scambio di vicine)
- **deporre e leggere matrice**
- **emettere e leggere un segnale** (una sostanza che diffonde)
- **cambiare identità**, in modo irreversibile

**Criterio di ammissione, unico:** si deve poter dire cosa si chiederebbe a una cellula vera o a un organoide per ottenerlo. Se la risposta è "nulla, è un dispositivo del simulatore", il comportamento non entra.

Due precedenti utili, da ricordare come esempi di ciò che è stato giustamente espulso: lo **stampo** (una griglia di voxel che ricordava la forma finale: nessuna cellula la possedeva) e l'**impalcatura** (uno stato cellulare dedicato, con privilegi di crescita propri, che teneva aperto un lume). Entrambi funzionavano; entrambi erano dispositivi del simulatore.

## 8. Nessun burattinaio

La disciplina si fa rispettare **per costruzione**, non per buona volontà: la regola di una fase riceve una *struttura di percezione* e nient'altro, e quella struttura viene costruita al passo 3 della passata, che non ha accesso allo stato globale. Il motore può calcolare quantità globali quanto vuole — ma solo dopo, al passo 9, per l'osservatore.

Ogni grandezza che una cellula usa va dichiarata come uno dei quattro canali:

1. **dentro di sé**: stato interno, vettori nel proprio riferimento, orologio proprio
2. **dalle vicine**: ciò che dichiarano, letto nel proprio riferimento (quindi con l'errore di rotazione relativa)
3. **dalla matrice**: alla propria posizione — densità, tensore delle fibre, etichette
4. **misura assoluta**: curvatura delle fibre, concentrazione di una sostanza

### I burattinai nel codice attuale, da non riportare

Tre regole che sembrano locali e non lo sono. Sono già state riscritte in forma locale nella tabella della sezione 5, ed è uno dei guadagni della riscrittura:

| burattinaio | cosa fa davvero | forma locale |
| --- | --- | --- |
| la fase 0 avanza su `S.calm` | un contatore **globale** di quiete del simulatore | sono quieta io, e nessuna vicina ha avuto eventi per T passi |
| la fase 2 avanza su un controllo di stabilità | confronto fatto su **tutte** le cellule insieme | il mio d cambia meno di ε per 60 passi |
| il raggio dell'arco cresce sul tick globale | orologio del simulatore, non della cellula | cresce sul mio orologio da quando sono entrata nella fase |

Un quarto caso è più subdolo e va capito bene, perché è *legittimo*: l'elezione viene risolta nel codice con una scansione di tutte le cellule alla ricerca di chi ha vinto. Il **predicato** però è locale ("il mio biglietto è imbattuto da 45 passi"), quindi la scansione è solo l'osservatore che guarda chi si è dichiarato. La distinzione da tenere: *una scansione globale per leggere è ammessa, una per decidere no.*

E un principio che vale più di tutti, trovato pagandolo: **`mantieni` tiene accesa la regola, non la marca.** Congelando la marca, le due pareti del solco — entrambe marcate e ora in contatto — si costringevano l'una contro l'altra, le coppie passavano da 30 a 199 e il tessuto si spaccava. Rivalutando il predicato a ogni passo, chi non è più sulla superficie libera smette da sé. In generale: **una marca è un'ipotesi sul presente, non un ricordo.**

## 9. Metodo

**Il motore si congela quando è validato.** Dopo, una modifica richiede un test *fisico* che falla — non una forma che si comporta male. Questa è la regola che corregge l'errore di fondo del lavoro fatto finora.

**Ogni meccanismo ha un interruttore**, e la batteria riporta cosa accade spegnendolo, una prova per meccanismo, con ripetizioni. Motivo: in questo lavoro cinque meccanismi epiteliali sono sembrati inerti per una sessione intera, e in realtà **non erano mai attivi** — il test di polarità usava una soglia assoluta su una densità che la cellula stessa produceva, e una cellula su dieci risultava polarizzata. Senza ablazione sistematica non se ne sarebbe saputo nulla.

**Casualità con seme, corse riproducibili.** La varianza fra prove ha confuso più di un'ablazione; senza seme non si distingue un effetto da un caso.

**Batteria con numeri attesi**, eseguita da sola a ogni modifica: crescita, quiete, riparazione di una ferita, invarianza di scala, e la stessa forma a due o tre raggi cellulari diversi.

**I risultati negativi entrano come test di regressione.** Sono la parte più preziosa e la più facile da perdere. I principali, da far valere come fatti:

- la necrosi da limite di diffusione **non** stabilizza un lume (la regola cancella la propria premessa)
- un campo di distanza min-plus inizializzato in basso **non** cresce oltre un salto
- la costrizione apicale fra sfere **non** piega un foglio
- la sezione 2D **non** predice il 3D per la flessione (una calotta è curva in due direzioni)
- il rimodellamento giunzionale **non** tiene una coppa invaginata
- il dominio differenziato tiene la coppa ma **non** la ripara

**Un esperimento di confronto esterno.** La stessa invaginazione epiteliale, nel nostro motore e in uno maturo (un modello a vertici, o EmbryoMaker). Se coincidono, la fisica tiene. Da verificare prima di scegliere: licenze e linguaggi effettivi — le indicazioni date a voce vengono dalla memoria e non sono state controllate.

## 10. Cosa si porta, cosa si butta

**Si porta** — e è la parte di valore:

- i **documenti**: l'algebra delle primitive, le regole di buona formazione, gli enunciati, e soprattutto i risultati negativi;
- le **primitive come idee**: informazione posizionale senza seme (il centro virtuale come ipotesi condivisa, eletta e non imposta), i quattro canali, la matrice come mappa persistente scritta dalle cellule e che sopravvive a loro, le elezioni con biglietto congelato, chiusura e successione, l'onda monotona delle fasi, `sempre mantieni` come legge di conflitto;
- le **misure** come riferimento: le forme che funzionano con i loro numeri (bordi, errore ρ, tempi di quiete, riparazione delle ferite);
- il **banco nel browser** come strumento di presentazione, con il motore JavaScript tenuto in vita *solo* per quello e non più sviluppato.

**Si butta**: il motore JavaScript come fondamento, e con lui `spingi`, `piatto`, `mangia`, lo `stampo`, l'`impalcatura`, `impegna` e i tempi tarati a mano. Quasi tutti hanno già un sostituto trovato e misurato (`costringi`, la fascia a due bordi, il dominio differenziato), e quelli che non ce l'hanno sono annotati come aperti.

### Struttura del repository

```
docs/      algebra, regole di buona formazione, enunciati, risultati negativi,
           questa specifica
motore/    Julia: geometria, energia, topologia, chimica, passata
primitive/ Julia: una regola per file, pura, con il suo predicato di fine
programmi/ le catene, come tabelle di fasi
prove/     batteria con i numeri attesi; i risultati negativi come regressioni
riferimento/ il motore JavaScript congelato, il banco, le misure storiche
confronto/ l'esperimento nel motore esterno
```

### Primo traguardo proposto

Non una forma: un **test fisico**. Un foglio epiteliale piatto, una fila di cellule che alza la tensione dell'anello apicale, e la misura della curvatura che ne risulta confrontata con la previsione analitica del modello a vertici. Se quello torna, la fisica tiene e si può ricominciare a costruire forme — e questa volta senza toccare il motore.
