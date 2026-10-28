/* Libreria delle macchine a stati morfogenetiche — versione a cerchi concentrici, senza semi.

   Ogni cellula porta un vettore v: "il centro è di là, a questa distanza", espresso con una
   bussola comune a tutto il tessuto (come la polarità planare). A ogni tick legge il vettore
   di ogni vicino, lo somma alla posizione del vicino misurata da lei, e fa la media con il suo.
   ρ = |v| è la sua distanza dal centro virtuale; i cerchi concentrici sono i luoghi a ρ costante.
   Il centro è un'ipotesi condivisa: nessuna cellula lo occupa, e se si distrugge la zona
   centrale gli altri continuano a puntare dove era.

   Perché la bussola: con le sole distanze (ρ dei vicini, senza direzioni condivise) il campo
   non ha un'ancora. Localmente vale "ρ cresce di 1 per unità", che è vero per la distanza da
   un punto ma anche da un cerchio o da una retta: il campo scivola verso un cratere e il disco
   cresce senza fine. Prima l'ancora era il seme; senza seme serve un asse comune.

   Matrice extracellulare: ogni cellula stende fibre tangenti al proprio cerchio e le rilegge.
   La fibra dà la direzione del centro (perpendicolare), la sua curvatura dà il verso. La lettura è
   relativa al corpo della cellula, quindi non dipende dalla bussola interna: dove c'è matrice i vicini
   si scambiano solo ρ, e la direzione la dà la fibra. La matrice diventa la bussola del tessuto.
   Le fibre portano anche un'etichetta con il ρ di chi le ha deposte: la matrice è una mappa completa,
   direzione e distanza, che sopravvive alle cellule.
   La bussola interna serve solo all'inizio, prima che ci siano fibre. La matrice resta dove è stata
   deposta: una ferita "solo cellule" trova la mappa intatta.

   Programma unico: una cellula esposta mette una figlia dove ρ cade nella fascia.
     Disco:  fascia [0, R]
     Anello: fascia [Rin, R]. Il lume è il fuori-fascia interno: nessuna impalcatura, nessuna
             cavità dichiarata. Una volta aperto resta aperto perché la parete si polarizza —
             la faccia verso il lume non ha matrice sotto, non aderisce, non vi manda figlie.
   Una ferita si ripara con la stessa regola della crescita.

   Ingressi x (calcolati dal motore, sempre locali):
     x.rho          la mia distanza dal centro virtuale
     x.R, x.Rin, x.r
     x.esposta      varco angolare tra i vicini >= 110°
     x.puoCrescere  esiste una direzione libera in cui la figlia avrebbe ρ nella fascia
                    e dove c'è matrice a cui attaccarsi, e non oltre la propria faccia apicale
     x.bloccata     non posso dividermi da "pazienza" tick
     x.fuori        ρ fuori da [Rin, R] (con margine r/2) da qualche tick
     x.compressa    schiacciata dalle vicine da qualche tick
     x.staccata     senza appigli: poca adesione, o troppa superficie libera attorno
   Stato della cellula c: c.state, c.migra (tick di rotolamento rimasti) */
(function(){

const BASE = [
  { id: "MIGRA",      token: "s-migra",      desc: "figlia appena nata, rotola sui vicini" },
  { id: "PROLIFERA",  token: "s-prolifera",  desc: "esposta, e una figlia cadrebbe nella fascia" },
  { id: "QUIESCENTE", token: "s-quiescente", desc: "circondata, inibizione da contatto" },
  { id: "BORDO",      token: "s-bordo",      desc: "esposta, ma ogni figlia uscirebbe dalla fascia" },
  { id: "APOPTOSI",   token: "s-apoptosi",   desc: "spinta fuori dalla fascia: si restringe e sparisce" }
];
const pos = (stati, map) => stati.map(s => ({ ...s, pos: map[s.id] }));

function nellaFascia(c, x) {
  if (c.state === "APOPTOSI" || c.state === "CAVITA") return c.state;
  if (x.fuori) return "APOPTOSI";
  if (c.migra > 0) return "MIGRA";
  if (x.puoCrescere) return "PROLIFERA";
  if (x.esposta) return "BORDO";
  return "QUIESCENTE";
}

const DISCO = {
  id: "disco",
  nome: "Disco",
  sintesi: "Fascia [0, R] attorno a un centro virtuale. Ogni cellula esposta cresce finché la figlia resta entro R. Nessun seme: la prima cellula è solo la prima.",
  parametri: ["R"],
  dividono: ["PROLIFERA"],
  muoiono: ["APOPTOSI"],
  stati: pos(BASE, { MIGRA: [190, 40], PROLIFERA: [190, 140], QUIESCENTE: [318, 90], BORDO: [318, 200], APOPTOSI: [72, 200] }),
  transizioni: [
    { da: "PROLIFERA",  a: "MIGRA",      cond: "divisione: la figlia nasce in MIGRA" },
    { da: "MIGRA",      a: "PROLIFERA",  cond: "fine rotolamento (30 tick), rivaluta" },
    { da: "PROLIFERA",  a: "QUIESCENTE", cond: "varco < 110°: circondata" },
    { da: "QUIESCENTE", a: "PROLIFERA",  cond: "si apre un varco (ferita): la riempie" },
    { da: "PROLIFERA",  a: "BORDO",      cond: "ogni figlia avrebbe ρ > R" },
    { da: "BORDO",      a: "PROLIFERA",  cond: "si apre un varco verso l'interno" },
    { da: "BORDO",      a: "APOPTOSI",   cond: "spinta oltre R + r/2" }
  ],
  decide(c, x) { return nellaFascia(c, x); }
};

const ANELLO = {
  id: "anello", nome: "Anello",
  sintesi: "La stessa regola del disco, su una fascia con due bordi: [Rin, R]. Una cellula che si trova sotto Rin è "
    +"fuori fascia esattamente come quella che supera R, e muore per la stessa ragione. Il lume non è un oggetto "
    +"dichiarato né tenuto aperto da un'impalcatura: è il fuori-fascia interno. Una volta aperto, resta aperto perché "
    +"la parete si polarizza — la faccia che guarda il lume non ha matrice sotto, non aderisce e non vi manda figlie.",
  parametri: ["R", "Rin"],
  dividono: ["PROLIFERA"],
  muoiono: ["APOPTOSI"],
  decide(c, x) { return nellaFascia(c, x); }
};

window.MACCHINE = { disco: DISCO, anello: ANELLO };
window.MORFO_COSTANTI = {
  r: 6, varco: 1.92, pienoDa: 11, pazienza: 120, fuoriTick: 6, giovinezza: 250,
  /* Figliata non vitale: se una figlia muore prima di 'giovane' tick, la madre smette di dividersi per 'pausaFigliata'.
     Chiude il ciclo futile di chi rifà all'infinito una figlia che la pressione della superficie spinge fuori. */
  giovane: 120, pausaFigliata: 400,
  /* Estrusione da affollamento: compressione = somma delle sovrapposizioni con le vicine, in frazione di raggio.
     Oltre soglia per compressioneTick la cellula viene espulsa e muore. Anoikis: con adesione totale sotto
     ancoraggioMin per anoikisTick, la cellula non ha appigli e muore. Sono misure meccaniche, non contatori di tentativi. */
  adesioneApicale: 0.1, polarita: 0.6, planare: 0.64, basale: 0.5, compressione: 0.9, compressioneTick: 60, ancoraggioMin: 0.6, stacco: 0.6, anoikisTick: 60,
  /* Consenso sul centro: quanto della media dei vicini si accetta a ogni tick. */
  consenso: { alfa: 0.5 },
  /* Matrice extracellulare: griglia di caselle (unità), quanto si deposita per tick in una casella vuota,
     quanto lentamente si rimodella una matrice matura, soglie di densità e coerenza per fidarsi, passo per leggere la curvatura e curvatura minima leggibile (rad).
     Vicino al centro (ρ < vicinoCentro) le fibre sono troppo curve: niente deposito. Si leggono solo da ρ ≥ leggiDa:
     all'inizio la matrice è fatta da poche cellule e un errore letto lì si propaga a tutto il tessuto.
     La curvatura dà anche ρ in assoluto (2·passo / angolo girato): pesoCurv è quanto conta rispetto a un
     vicino. Senza un'àncora ρ, propagato solo per differenze, andrebbe alla deriva.
     L'àncora vera è la memoria: le fibre portano il ρ di chi le ha deposte, e chi arriva dopo lo rilegge
     (pesoMemoria, relativo a un vicino). Resta anche quando il nucleo muore, come nell'anello. */
  matrice: { griglia: 4, deposito: 0.05, rimodella: 0.002, soglia: 0.2, coerenza: 0.5, passo: 10, /* portata della lettura: 10 unità con cellule di raggio 6, cioè meno di un diametro. È la portata di un prolungamento, non una lettura a distanza arbitraria. */ curvMin: 0.03, vicinoCentro: 12, leggiDa: 25, pesoCurv: 0, curvMax: 400, pesoMemoria: 1, nucleo: false,
    /* Controllo di curvatura: la lettura locale (2·passo/angolo) è troppo rumorosa per correggere ρ cellula per cellula,
       ma la sua media su una zona ampia no. curvMedia: quanto pesa la mia lettura nella media con i vicini (più basso = zona
       più ampia); curvBias: lo scarto sistematico della lettura in un tessuto sano; curvGuadagno: quanto in fretta ρ segue la
       curvatura media; fuoriMargine e fuoriTick: oltre R + margine secondo le fibre, per tanti tick, la cellula muore e
       cancella la matrice sotto di sé, così il buco guarisce con il ρ giusto. */
    curvLocale: false, ancoraggio: 0.5, curvRhoMin: 20, curvMedia: 0.02, curvBias: 0.03, curvGuadagno: 0.005, fuoriMargine: 15, fuoriTick: 60,
    /* Ellisse: le fibre portano anche l'asse di polarità di chi le ha deposte; asse = quanto la lettura riallinea
       l'asse della cellula a ogni tick; pcp = quanto la cellula allinea il proprio asse a quello delle vicine che tocca
       (polarità planare). Taglio: larghezza della striscia che muore tra i due frammenti. */
    asse: 0.2, pcp: 0.3, striscia: 15, raggioLettura3d: 9 }
};
})();
