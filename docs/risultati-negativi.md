# Risultati negativi

Cose provate che **non** funzionano, con il motivo. Sono la parte più preziosa del lavoro e la
più facile da perdere: nel motore nuovo vanno fatte valere come **test di regressione**, perché
la tentazione di riprovarle tornerà.

Ognuna è stata misurata, non supposta.

---

### 1. La necrosi da limite di diffusione non stabilizza un lume

**Cosa:** una cellula troppo lontana da una superficie bagnata dal mezzo muore, come il centro
necrotico di uno sferoide. Doveva aprire la cavità di un guscio senza dichiararla.

**Esito:** 6 prove a tre soglie diverse (20, 30, 45 unità): nessuna apre un buco, nessuna
raggiunge la quiete. Solo logoramento perpetuo.

**Perché, strutturalmente:** *la regola cancella la propria premessa.* Appena il centro muore,
il tessuto che lo circondava non è più profondo, ricresce, muore di nuovo. Il limite di diffusione
spiega come una cavità *si apre*, non come *resta aperta*.

**Regressione:** con una regola di morte per profondità e nessun altro meccanismo, un corpo pieno
non deve raggiungere un lume stabile.

---

### 2. Un campo di distanza min-plus inizializzato in basso non cresce

**Cosa:** il campo di profondità era inizializzato a 0 su tutte le cellule e relassato verso l'alto.

**Esito:** misurava 8 dove doveva misurare 54. Il minimo resta incollato a zero, e la distanza non
supera mai un salto.

**Correzione:** inizializzare le non-esposte a un valore grande e far *scendere* il campo. Corretto
così, in un disco di raggio 60 misura 0 al bordo e 54 al centro.

**Regressione:** ogni campo di distanza deve essere verificato sul suo profilo radiale, non sul fatto
che "sembra funzionare".

---

### 3. La costrizione apicale fra sfere non piega un foglio

**Cosa:** `costringi` implementato come accorciamento della lunghezza di riposo fra due cellule che
costringono la stessa faccia.

**Esito:** in un monostrato non si forma nessuna piega. Il foglio si accorcia, non si curva.

**Perché:** la costrizione apicale è **sub-cellulare** — la cellula diventa a cuneo, stretta
all'apice e larga alla base. Con una sola lunghezza di riposo per contatto la coppia si accorcia
simmetricamente e non produce nessuna coppia di forze.

**Aggiramento usato (da non riportare):** una parete di due file, dove la fila esterna costringe e
quella interna no — una striscia bimetallica. Funziona ma è un trucco imposto da una
rappresentazione sbagliata.

**Nel motore nuovo:** sparisce per rappresentazione. Con vertici e facce, `costringi` è alzare la
tensione dell'anello apicale, e il cuneo è il minimo dell'energia a volume costante.

---

### 4. La sezione 2D non predice il 3D per la flessione

**Cosa:** sviluppare i meccanismi in 2D (economico) e portarli in 3D.

**Esito:** con i parametri che in 2D davano una coppa profonda (fronte da ρ≈57 a 36), in 3D il
fronte non si muove (resta a 52–54, calotta di 12–28 cellule).

**Perché:** un arco piano è curvo in una direzione, una calotta sferica in due, e resiste alla
flessione molto di più.

**Conseguenza di metodo:** il 2D resta utile per **trovare** i meccanismi — è stato decisivo per
capire la faccenda del cuneo — ma non per **tararli**.

---

### 5. Il rimodellamento giunzionale non tiene una coppa invaginata

**Cosa:** un contatto tenuto accorciato abbastanza a lungo adotta la lunghezza accorciata come
nuova lunghezza di riposo (per contatto, non per cellula). Doveva rendere la piega permanente
senza tenere accesa la costrizione.

**Esito:** la coppa si richiude, 2 prove su 2. Il tessuto si affolla (n da 100 a 113, strain delle
coppie a −0.48) invece di tenere la forma.

---

### 6. Il dominio differenziato tiene la coppa ma non la ripara

**Cosa:** le cellule che hanno costretto diventano una sottopopolazione differenziata — post-mitotica
e non più giudicata dalla banda della forma madre. Regola: "resta come sei".

**Esito positivo:** con la costrizione **spenta**, la coppa resta ferma 4000 passi, 3 prove su 3
(canale a 43, 51, 54; dominio di 40, 42, 32 cellule stabili). Questo **sostituisce `impegna`**.

**Esito negativo:** asportata metà del dominio, **non ricresce** (36 → 13 cellule, ferme) e la coppa
è perduta, 2 prove su 2: il tessuto circostante riempie il solco perché per la sua banda quello
spazio è suo.

**Il compromesso, esplicito:**

| meccanismo | memoria richiesta | ripara? |
| --- | --- | --- |
| `stampo` (griglia di voxel) | globale, nessuna cellula la possiede | sì |
| `sempre mantieni(costringi)` | nessuna, ma il motore resta acceso | sì, finché acceso |
| dominio differenziato | nessuna, motore spento | **no** |

---

### 7. Un dominio con campo proprio: la fascia non si applica prima della convergenza

**Cosa:** dare alla sottopopolazione un polo eletto fra i suoi, una geodetica attraverso il solo
dominio, e una fascia propria — per comprare *sia* la stabilità *sia* la riparazione.

**Esito:** non funziona ancora (la crescita del dominio non riparte dopo la ferita), ma due cose
imparate:

- **alla prima iterazione tutte le distanze sono ignote, e "ignota" non è un motivo per morire.**
  Nella prima versione il dominio si estingueva in 200 passi. Vale come regola generale: *una
  distanza ancora in discesa non autorizza nessuna decisione irreversibile.*
- **il raggio della fascia del dominio non è un parametro libero:** dev'essere il raggio a cui
  l'arco si è fermato, lo stesso numero in un altro riferimento. Un dominio **eredita la propria
  taglia dal processo che l'ha fatto**.

---

### 8. Cinque meccanismi epiteliali sembravano inerti: non erano mai attivi

**Cosa:** divisione planare, deposizione basale della lamina, adesione apicale debole, estrusione
per compressione e anoikis, distinzione fra superficie bagnata e interstizio. Spenti uno per uno
sull'anello, nessuno cambiava il risultato oltre il rumore.

**La vera causa:** il test di polarità apicale confrontava la densità di matrice oltre la faccia
libera con una **soglia assoluta** — ma una cellula depone matrice attorno a sé, e la soglia non
veniva mai superata. Una cellula su dieci risultava polarizzata.

**Correzione:** la polarità è un **confronto fra i due lati**, non una soglia; e in un monostrato
nemmeno il baricentro dello spazio libero basta, perché la cellula ha facce libere sui due lati e
la media si annulla. La versione che funziona legge il **gradiente della matrice** attorno a sé e
punta l'apicale dove la lamina manca: 60 cellule polarizzate su 79.

**Conseguenza di metodo:** ogni meccanismo deve avere un interruttore e la batteria deve riportare
cosa accade spegnendolo. Senza ablazione sistematica un meccanismo morto sembra un meccanismo inutile.

---

### 9. La crescita dipendente dall'ancoraggio era biologicamente falsa

**Cosa:** una figlia può nascere solo dove c'è già matrice a cui attaccarsi.

**Perché è falsa:** un epitelio che si espande **depone lamina nuova man mano**. Pretendere la
matrice prima della figlia impedisce al tessuto di allargarsi: il bordo rientrava di una dozzina
di unità (48–63 invece di 52–64).

**Correzione:** l'ancoraggio vincola la **sopravvivenza**, non l'espansione.

---

### 10. `mantieni` tiene accesa la regola, non la marca

**Cosa:** per mantenere una deformazione, congelare l'insieme delle cellule marcate.

**Esito:** le due pareti del solco — entrambe marcate e ora in contatto — si costringono l'una
contro l'altra: le coppie passano da 30 a 199 e **il tessuto si spacca in due**.

**Correzione:** rivalutare il predicato a ogni passo, così chi non è più sulla superficie libera
smette da sé. In generale: **una marca è un'ipotesi sul presente, non un ricordo.**

---

### 11. Il secondo polo non si trova con un predicato locale di massimo

**Cosa:** dalla legge di simmetria (un'elezione rompe la simmetria solo se il biglietto è costante
sulle orbite) segue che la selezione di P2, il cui biglietto d₁ ha massimo determinato nell'antipodo,
è *equivariante* — quindi l'apparato dell'elezione sembrava inutile e bastava "sono un massimo
locale di d₁".

**Esito:** conclude più in fretta (122 passi contro 166) ma **si dichiarano da 3 a 8 candidati per
prova**, e l'angolo fra i poli peggiora da 170° ± 6° a 135° ± 34°.

**Perché:** un'elezione fa **due lavori distinti**, e la legge di simmetria parla solo del primo —
(1) rompere la simmetria quando il biglietto è costante sulle orbite, (2) trovare il **massimo
globale** di un biglietto che varia ma è rumoroso e ha massimi locali spuri. Per P2 il lavoro 1 non
serve, il lavoro 2 sì.

**Da ricordare come distinzione permanente:** **rompere ≠ massimizzare.**

---

### 12. Una regolina scartata con la misura: il margine di nascita

**Cosa:** per evitare il ciclo futile (una cellula che rifà la stessa figlia nello stesso posto
all'infinito) si era provato a restringere il margine di nascita a R − r/2.

**Esito:** risolveva gli stalli ma **restringeva la sfera del 12%** e la rendeva meno rotonda.
Scartata in favore della quiescenza contata su eventi reali e della pausa della figliata non vitale.

---

### 13. Un pezzo di foglio 3D con bordo libero non ha equilibrio

**Cosa:** un pezzo di foglio di prismi esagonali (7 cellule) con il bordo libero, portato al
minimo dell'energia con le tensioni uguali su tutte le facce.

**Esito:** non converge. Gli spigoli del bordo si accorciano fino alla lunghezza minima e lì il
minimizzatore si ferma, con una forza residua di 0.17 nN.

**Perché:** al bordo le facce laterali appartengono a una cellula sola, quindi niente bilancia la
loro tensione, e ogni parete libera si torce. Misurato (29 ottobre, 19 cellule, K = 0.1): gli
spigoli che collassano sono quelli di bordo delle cellule d'angolo, e collassano **a turno in
alto o in basso**: 43–44 apicale 0.24 µm e basale 8.3 µm, 45–46 apicale 8.3 µm e basale
0.36 µm, e così via. La cellula di bordo vuole un lato in meno all'apice e un altro in meno alla
base, come la cellula isolata che si arrotonda.

**Conseguenza:** nel guscio a prismi non esiste un'operazione al bordo che risolva il caso. Una
fusione di vertici fatta insieme sulle due superfici schiaccerebbe uno spigolo di 8 µm (come lo
scambio nel #14). È lo stesso limite dei #14 e #15: servirebbero superfici apicale e basale con
topologie proprie.

**Regressione:** `test/guscio.jl`. Le prove di equilibrio 3D usano il foglio periodico, che non
ha bordi.

---

### 14. Nel guscio a prismi non c'è uno scambio di vicine solo apicale

**Cosa:** due cellule vicine stringono l'anello apicale (Λ = 4 nN, foglio periodico 6 × 6), per
far chiudere la giunzione fra loro e provocare uno scambio di vicine spontaneo.

**Esito:** la giunzione apicale si chiude (0.05 µm) ma quella basale resta lunga 9.3 µm. Lo
scambio, che nel guscio a prismi cambia insieme le due superfici, schiaccia e ruota anche lo
spigolo basale: l'energia salta da 13040 a 21424 nN·µm e il tessuto non si riassesta (lo
scambio si ripete avanti e indietro).

**Perché:** nel guscio ogni cellula è un prisma, con gli stessi vicini sopra e sotto. Le cellule
vere possono avere vicini diversi all'apice e alla base (la forma detta scutoide); il modello
non le rappresenta.

**Correzione:** lo scambio si fa solo quando sono corti tutti e due gli spigoli della parete,
apicale e basale. Una giunzione solo apicale che si chiude resta un caso aperto: chiede cellule a
scutoide o una contrattilità basale.

**Regressione:** `test/scambi.jl`.

---

### 15. Nel guscio a prismi una cellula che chiude solo l'apice non si può estrudere

**Cosa:** una cellula stringe forte l'anello apicale (Λ = 8 nN, foglio periodico 6 × 6), come
una cellula che si prepara a uscire dal foglio.

**Esito:** l'apice si chiude (0.015 µm²) ma la base si allarga (175 µm²): è la cellula a
bottiglia. Non si può estrudere: l'estrusione, che nel guscio fonde insieme i vertici apicali e
quelli basali, farebbe sparire di colpo una base di 175 µm².

**Perché:** come il #14, nel guscio a prismi apice e base cambiano insieme. Una cellula vera che
chiude l'apice esce dal lato basale o si invagina; il modello non rappresenta né l'una né
l'altra cosa.

**Cosa invece funziona:** l'estrusione di una cellula piccola tutta intera, cioè apoptotica, che
riduce il suo volume, perde lati con gli scambi di vicine fino al triangolo ed esce
(`test/estrusione.jl`).

**Regressione:** `test/estrusione.jl`.

### 16. Nel modello generale l'apoptosi arriva al prisma triangolare, ma non esce

**Cosa:** la cellula centrale di un pezzo di foglio a bordo bloccato (19 cellule) entra in
apoptosi nel modello generale (`11_operazioni.jl`, scutoidi ammessi): il volume bersaglio scende
con tempo 8 minuti, passo implicito di 1 minuto, operazioni elementari dopo ogni passo.

**Esito:** la cellula si restringe tutta insieme e i sei lati arrivano alla soglia nello stesso
momento. Alla passata 43 sei operazioni da spigolo a triangolo fanno triangoli la faccia apicale
e quella basale, ma la cellula ha ancora sei pareti alte 15 µm (scutoidi); alla 44 il volume è
transitoriamente negativo (−0.27 µm³); alla 47 la cellula è un prisma triangolare. Lì si ferma:
nel modello generale non c'è un'estrusione, e la meccanica resta bloccata alla lunghezza minima.
L'esperimento ha anche trovato due difetti delle operazioni di `11` (facce che restavano con due
vertici), corretti con le loro prove.

**Perché:** le operazioni scelgono lo spigolo più corto; in una cellula che si restringe in modo
uniforme sono tutti corti insieme, e le pareti nuove, di 1 µm, sono più grandi della cellula.
Manca l'operazione che toglie una cellula piccola (la T2 dei modelli a vertici 3D), e una regola
che tenga le pareti nuove proporzionate alla cellula.

**Cosa invece funziona:** nel guscio a prismi la stessa apoptosi esce dal foglio alla passata 30,
con uno scambio di vicine per passata e l'estrusione di `09` (`test/popolazione.jl`,
`prova_la_morte_e_un_estrusione`). Divisione e morte nella passata chiedono per ora il guscio a
prismi (`topologia = :prismi`).

**Regressione:** `test/popolazione.jl`, `prova_risultato_negativo_16`.

### 17. Una passata troppo lunga fa oscillare un circuito di retroazione

**Cosa:** l'organizzatore senza seme (`attiva_inibisci`, Nodal/Lefty con i numeri di Müller et
al. 2012, κ = 0.29) su un foglio di 1672 cellule, con passate di 30 minuti.

**Esito:** il livello medio di Nodal in tutto il tessuto oscilla con un periodo di circa 20 ore
e un'ampiezza che cresce (fra 0.7 e 1.4 volte lo stato uniforme in 50 ore); una macchia nasce
lo stesso, ma la catena che ne legge la soglia si accende e si spegne. Con passate di 5 minuti
lo stesso livello resta fra 0.95 e 1.04.

**Perché:** le secrezioni di ogni cellula sono decise all'inizio della passata e restano fisse
fino alla successiva, mentre la chimica avanza in sotto-passi. Per un circuito di retroazione è
un passo di Eulero esplicito: lo stato uniforme linearizzato ha autovalori −δ ± iω, e il passo è
stabile solo se la passata è più corta di 2δ/ω². Con κ = 0.29 lo stato uniforme è appena
smorzato (δ = 4·10⁻⁵/min, ω = 0.0066/min, periodo 16 ore) e il limite è 1.9 minuti; con 30
minuti ogni passata amplifica dell'1.8%. La specifica voleva passate di decine di secondi.

**Cosa invece funziona:** κ = 0.35 (lo stato uniforme si smorza in circa 57 ore, limite 12.9
minuti) e passate di 5 minuti. Le prove dell'organizzatore della prima versione (κ = 0.29, 30
minuti) sono state rifatte: le conclusioni reggono, i numeri del contrasto erano gonfiati.

**Regressione:** `test/organizzatore.jl`, `prova_risultato_negativo_17`.
