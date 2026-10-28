/* Primitive topologiche sul motore 3D. Ogni primitiva è una regola locale con un evento di fine.
   Il programma (livello esterno) è una fase per cellula; una cellula passa alla fase successiva
   quando rileva localmente l'evento di fine della precedente:
     CRESCI ──quiete──> ELEGGI1 ──eletto P1 + campo d₁ stabile──> ELEGGI2 ──eletto P2──> POLI
   Tutto è locale: ogni cellula legge solo sé stessa e le vicine (c.nb, aggiornato dal motore). */
(function(){
const P={ quiete:120,        // tick senza attività entro la propria portata per dichiarare quiete
          scadenza:40,       // passaggi oltre i quali un biglietto è scartato
          durata:45,         // tick da imbattuta per essere eletta
          stabile:60,        // tick senza variazioni di un campo per considerarlo stabile
          raggioDisco:27,    // distanza dal polo entro cui la superficie diventa disco
          attesaSpinta:30,   // tick in DISCHI prima di spingere (la marca si diffonde)
          spinta:0.6,        // spostamento attivo per tick delle cellule di disco (sotto ~0.5 il tessuto pieno non cede)
          piatto:0.5,        // quanto in fretta una cellula di disco va alla quota media delle compagne
          pesoGradiente:0.5, // peso della lettura locale del gradiente rispetto a una compagna
          tieni:20,          // tick di contatto continuo che chiudono la spinta
          colla:4,           // adesione tra i due dischi dopo il contatto, relativa a quella normale
          durataAttacca:300, // tick di fusione prima di staccarsi dal resto
          coesioneEpi:3,     // adesione tra cellule dell'epitelio unite da giunzioni, relativa a quella normale
          apicale:0.2,       // adesione tra epiteli che si toccano per la prima volta (pareti opposte di un lume)
          durataStacca:200,  // tick tra lo stacco e la perforazione
          fora:true,         // il cuscino muore e apre il foro
          assesto:200,       // tick dopo il foro prima di stampare la forma
          mangia:20,         // tick davanti a un disco dopo i quali una cellula del corpo muore (0 = il disco spinge soltanto)
          fuoriStampo:60,    // tick fuori dallo stampo prima di morire
          programma:"toro",  // "toro": due poli, i dischi si fondono e aprono un foro. "gastrula": un polo, un disco che invagina
          rhoFine:0.3,       // gastrula: frazione di R sotto la quale il disco ha finito di scendere
          reclutaDisco:false,// "sempre marca": l'epitelio dietro al disco entra nel disco e riveste il canale
          reclutaTocchi:2,   // quante cellule di disco deve toccare per essere reclutata
          reclutaDietro:0.5, // di quanto il suo ρ deve superare il loro (in raggi di cellula) per contare come "dietro"
          canale:false,      // il disco marca il canale nella matrice mentre scende, e chi ci entra muore
          raggioCanale:1.0,  // raggio del canale marcato, in raggi di cellula
          fuoriCanale:30,    // tick dentro al canale prima di morire
          };
const dist=(a,b)=>{const x=a[0]-b[0],y=a[1]-b[1],z=a[2]-b[2];return Math.sqrt(x*x+y*y+z*z);};
// direzioni di prova per le figlie nel toro (spirale di Fibonacci) e passo di divisione, come nel motore
const DIRS=[];{const N=48,g=Math.PI*(3-Math.sqrt(5));for(let i=0;i<N;i++){const y=1-2*(i+.5)/N,r=Math.sqrt(1-y*y);DIRS.push([Math.cos(g*i)*r,y,Math.sin(g*i)*r]);}}
const SW=150;
function Primitive(E){
  const S=E.S, R={eventi:[],tick0:{}};
  // risoluzione: le grandezze di cellula scalano con il raggio r0 del motore; i tempi di attraversamento di una
  // cellula (mangia) scalano con lui perché il disco avanza a velocità assoluta
  const r0=E.r||6, sc=r0/6, PASSO=2*r0*0.85, SV=4*sc, SN=Math.ceil(2*SW/SV);
  // lo stampo: occupazione su una griglia di 4 unità; una cellula stampa una palla di raggio r+1 attorno a sé
  const ST=new Uint8Array(SN*SN*SN), si=x=>Math.max(0,Math.min(SN-1,Math.floor((x+SW)/SV)));
  const leggi=p=>ST[(si(p[2])*SN+si(p[1]))*SN+si(p[0])];
  function perVox(p,rad,fn){ for(let k=si(p[2]-rad);k<=si(p[2]+rad);k++) for(let j=si(p[1]-rad);j<=si(p[1]+rad);j++) for(let i=si(p[0]-rad);i<=si(p[0]+rad);i++){
      const x=-SW+(i+.5)*SV-p[0], y=-SW+(j+.5)*SV-p[1], z=-SW+(k+.5)*SV-p[2]; if(x*x+y*y+z*z<=rad*rad) fn((k*SN+j)*SN+i); } }
  const stampa=p=>perVox(p,r0+sc,v=>{if(ST[v]!==2)ST[v]=1;}); // il canale (2) vince sul corpo (1): chi lo invade non ne cancella la memoria
  // il disco, scendendo, marca nella matrice il vuoto che lascia dietro di sé: "qui è canale" (valore 2).
  const stampaVuoto=p=>perVox(p,r0*P.raggioCanale,v=>{if(ST[v]!==1)ST[v]=2;});
  // una ferita che cancella la matrice cancella anche lo stampo
  S.onClear=(p,rad)=>perVox(p,rad,v=>{ST[v]=0;});
  const dying=E.dying;
  const init=c=>{ if(c.fase) return; c.fase="CRESCI"; c.q=0; c.e1={bv:-Infinity,ba:0,win:0}; c.e2={bv:-Infinity,ba:0,win:0}; c.d1=Infinity; c.d1st=0; c.jit=Math.random(); c.polo=0; };
  const log=(nome)=>{ if(!(nome in R.tick0)){ R.tick0[nome]=S.tick; R.eventi.push(nome+" a t="+S.tick); } };
  // un'elezione: biglietto più alto vince; il biglietto viaggia un contatto per tick con un'età che scade
  function elezione(cells,key,ticket,cand){
    for(const c of cells){ if(dying(c)||!cand(c)) continue; const e=c[key], my=ticket(c);
      let bv=-Infinity,ba=0; for(const o of c.nb){ if(dying(o)||!cand(o)) continue; const oe=o[key]; if(oe.ba+1>P.scadenza) continue;
        if(oe.bv>bv||(oe.bv===bv&&oe.ba+1<ba)){bv=oe.bv;ba=oe.ba+1;} }
      e.seen=bv; if(my>=bv){e.nbv=my;e.nba=0;} else {e.nbv=bv;e.nba=ba;} e.win=(my>=bv)?e.win+1:0; }
    for(const c of cells){ if(dying(c)||!cand(c)) continue; const e=c[key]; e.bv=e.nbv; e.ba=e.nba; } }
  function step(){
    const cells=S.cells; cells.forEach(init);
    // successione: un polo che muore passa la marca alla vicina viva più vicina (la marca non deve morire con la cellula)
    // dopo l'impegna i poli hanno finito il loro compito: niente successione (altrimenti la marca sopravvive al foro)
    for(const c of cells){ if(!c.polo||!dying(c)||c.impegnata) continue; let best=null,bd=Infinity;
      for(const o of c.nb){ if(dying(o)||o.polo) continue; const q=dist(o.x,c.x); if(q<bd){bd=q;best=o;} }
      if(best){ best.polo=c.polo; if(c.polo===1){best.d1=0;} R.successioni=(R.successioni||0)+1; } c.polo=0; }
    // primitiva: attendi quiete. q = tick da cui nessuna cellula attiva è a portata (si propaga di un contatto per tick)
    // "attiva" vuol dire che è successo qualcosa, non che qualcuno potrebbe agire: una cellula che sta in PROLIFERA
    // senza mai dividersi non è attività. Eventi reali: rotolamento (MIGRA), divisione, nascita, inizio di una morte.
    for(const c of cells){ if(dying(c)) continue; const attiva=c.state==="MIGRA"||(S.tick-(c.ev==null?-999:c.ev))<=1;
      let m=Math.min(255,c.q+1); for(const o of c.nb){ if(!dying(o)&&o.q+1<m) m=o.q+1; } c.qn=attiva?0:m; }
    for(const c of cells) if(!dying(c)) c.q=c.qn;
    for(const c of cells){ if(dying(c)) continue;
      if(c.fase==="CRESCI"&&c.q>=P.quiete){ c.fase="ELEGGI1"; log("quiete"); }
      // la fase si propaga: chi tocca una cellula già avanti la segue (l'evento è arrivato anche qui)
      if(c.fase==="CRESCI"&&c.nb.some(o=>o.fase&&o.fase!=="CRESCI")) c.fase="ELEGGI1";
      // il biglietto si congela entrando nell'elezione: un valore che cambia a ogni tick non fa vincere nessuno
      if(c.fase!=="CRESCI"&&c.t1==null) c.t1=c.rho+0.5*c.jit; }
    // primitiva: elegge(P1, biglietto = ρ) tra le cellule in ELEGGI1 o oltre
    // l'elezione si chiude: chi ha già ricevuto il campo d₁ (cioè sa che P1 esiste) non compete più.
    // Senza chiusura, le cellule nate dopo per ricambio rieleggevano un secondo polo.
    const in1=c=>c.fase==="ELEGGI1"&&!isFinite(c.d1);
    elezione(cells,"e1",c=>c.t1,in1);
    for(const c of cells){ if(dying(c)||!in1(c)) continue;
      if(!c.polo&&c.e1.win>=P.durata){ c.polo=1; log("eletto P1"); } }
    // primitiva: campo(d₁, da P1): distanza geodetica lungo i contatti. Serve dall'elezione di P1 fino al contatto dei
    // dischi (la spinta legge il gradiente); prima e dopo non si calcola, per risparmiare tempo
    const duePoli=P.programma!=="gastrula";
    const serveD1=R.tick0["eletto P1"]!=null&&R.tick0["fine della spinta"]==null, serveD2=duePoli&&R.tick0["eletto P2"]!=null&&R.tick0["contatto stabile: attacca"]==null;
    if(serveD1) for(const c of cells){ if(dying(c)) continue; let d=c.polo===1?0:Infinity; for(const o of c.nb){ if(dying(o)) continue; const v=o.d1+dist(o.x,c.x); if(v<d) d=v; }
      c.d1n=d; }
    if(serveD1) for(const c of cells){ if(dying(c)) continue; if(Math.abs(c.d1n-c.d1)>0.5||!isFinite(c.d1n)) c.d1st=0; else c.d1st++; c.d1=c.d1n;
      // nella gastrula il campo d₁ stabile porta direttamente alla marca del disco: non c'è un secondo polo
      if(c.fase==="ELEGGI1"&&isFinite(c.d1)&&c.d1st>=P.stabile){ c.fase=duePoli?"ELEGGI2":"POLI"; log("campo d1 stabile"); }
      if(c.fase==="ELEGGI1"&&c.nb.some(o=>o.fase==="ELEGGI2"||o.fase==="POLI")) c.fase=duePoli?"ELEGGI2":"POLI";
      if((c.fase==="ELEGGI2"||c.fase==="POLI")&&c.t2==null) c.t2=(isFinite(c.d1)?c.d1:0)+0.5*c.jit; }
    // primitiva: elegge(P2, biglietto = d₁)
    if(duePoli) for(const c of cells){ if(dying(c)) continue; c.p2kn=c.polo===2||c.nb.some(o=>!dying(o)&&o.p2k); }
    for(const c of cells) if(!dying(c)) c.p2k=c.p2kn;
    const in2=c=>duePoli&&c.fase==="ELEGGI2"&&!c.p2k;
    if(duePoli) elezione(cells,"e2",c=>c.t2,in2);
    for(const c of cells){ if(dying(c)||!in2(c)) continue;
      if(!c.polo&&c.e2.win>=P.durata){ c.polo=2; c.fase="POLI"; log("eletto P2"); } }
    for(const c of cells){ if(!dying(c)&&c.fase==="ELEGGI2"&&(c.p2k||c.nb.some(o=>o.fase==="POLI"||o.fase==="DISCHI"))) c.fase="POLI"; }
    // primitiva: campo(d₂, da P2)
    if(serveD2) for(const c of cells){ if(dying(c)) continue; if(c.d2==null){c.d2=Infinity;c.d2st=0;} let d=c.polo===2?0:Infinity;
      for(const o of c.nb){ if(dying(o)||o.d2==null) continue; const v=o.d2+dist(o.x,c.x); if(v<d) d=v; } c.d2n=d; }
    if(serveD2) for(const c of cells){ if(dying(c)) continue; if(Math.abs(c.d2n-c.d2)>0.5||!isFinite(c.d2n)) c.d2st=0; else c.d2st++; c.d2=c.d2n; }
    // primitiva: marca(D1, d₁ ≤ r_d ∧ superficie) ‖ marca(D2, d₂ ≤ r_d ∧ superficie)
    // superficie = cellula esposta (BORDO). La marca si decide una volta, entrando nella fase.
    for(const c of cells){ if(dying(c)) continue;
      const pronta=c.fase==="POLI"&&(duePoli?(isFinite(c.d2)&&c.d2st>=P.stabile):true), propaga=c.fase==="POLI"&&c.nb.some(o=>o.fase==="DISCHI");
      if(pronta||propaga){ c.fase="DISCHI"; log(duePoli?"dischi marcati":"disco marcato");
        const sup=c.state==="BORDO"; c.disco=sup&&c.d1<=P.raggioDisco?1:(duePoli&&sup&&c.d2<=P.raggioDisco?2:0); c.tD=0; } }
    // impegna ; spingi(D1, verso P2) ‖ spingi(D2, verso P1) fino a contatto(D1, D2) ; attacca(D1, D2) ; stacca(D1 ∪ D2, resto)
    // Dopo qualche tick in DISCHI (la marca si è diffusa) la cellula entra in SPINGI e smette di crescere e morire.
    const DOPO={SPINGI:1,ATTACCA:1,STACCA:1,FORA:1,MANTIENI:1}, ord={SPINGI:1,ATTACCA:2,STACCA:3,FORA:4,MANTIENI:5};
    for(const c of cells){ if(dying(c)) continue;
      // entrando in SPINGI la superficie diventa epitelio: un foglio coeso che il disco trascina con sé,
      // così il canale dietro il disco resta rivestito e l'interno non lo richiude
      if(c.fase==="DISCHI"&&(++c.tD>=P.attesaSpinta||c.nb.some(o=>DOPO[o.fase]))){ c.fase="SPINGI"; c.epi=c.state==="BORDO"||!!c.disco; c.giunzioni=null; log("impegna, inizia la spinta"); }
      // giunzioni: le vicine di epitelio al primo tick in SPINGI. La coesione forte vale solo per loro; un contatto nuovo
      // tra epiteli (le pareti opposte del canale) aderisce poco, come il lato apicale delle cellule vere
      if(c.fase==="SPINGI"&&c.epi&&!c.giunzioni) c.giunzioni=new Set(c.nb.filter(o=>!dying(o)&&o.epi));
      // sempre marca(D1, epitelio che tocca il disco da dietro): il foglio segue il disco e riveste il canale.
      // "da dietro" = ρ maggiore di quello della cellula di disco che tocca, così la marca risale la scia invece di
      // allargarsi sulla superficie esterna, dove i ρ sono tutti uguali. È marca con "sempre", non una primitiva nuova.
      if(P.reclutaDisco&&c.fase==="SPINGI"&&c.epi&&!c.disco){ let k=0,d0=0;
        for(const o of c.nb){ if(dying(o)||!o.disco||o.fase!=="SPINGI") continue; if(c.rho>o.rho+P.reclutaDietro*(E.r||6)){ k++; d0=o.disco; } }
        if(k>=P.reclutaTocchi){ c.disco=d0; c.u=null; } }
      if(DOPO[c.fase]) c.impegnata=true; }
    // direzione di spinta. Locale: il gradiente del campo dell'altro polo (−∇d₂ per D1), stimato dalle vicine,
    // poi allineato con le compagne di disco (polarità planare). Per simmetria la media del disco è l'asse.
    const grad=(c,k)=>{ let g=[0,0,0],n=0; for(const o of c.nb){ if(dying(o)||!isFinite(o[k])) continue; const d=[o.x[0]-c.x[0],o.x[1]-c.x[1],o.x[2]-c.x[2]], l2=d[0]*d[0]+d[1]*d[1]+d[2]*d[2]; if(l2<1e-6) continue;
      const s=(o[k]-c[k])/l2; g[0]+=s*d[0];g[1]+=s*d[1];g[2]+=s*d[2];n++; } const l=Math.hypot(...g); return n&&l>1e-9?[g[0]/l,g[1]/l,g[2]/l]:null; };
    const spinge=c=>!dying(c)&&c.fase==="SPINGI"&&c.disco, piana=c=>!dying(c)&&(c.fase==='SPINGI'||c.fase==='ATTACCA')&&c.disco&&c.u; // staccato, il cuscino è libero
    for(const c of cells){ if(!spinge(c)) continue; let loc=null;
      if(duePoli){ const gl=grad(c,c.disco===1?"d2":"d1"); loc=gl?[-gl[0],-gl[1],-gl[2]]:null; }
      else { // gastrula: il disco scende verso il centro virtuale, che ogni cellula conosce già (c.c − c.x)
        const d=[c.c[0]-c.x[0],c.c[1]-c.x[1],c.c[2]-c.x[2]], l=Math.hypot(...d); if(l>1e-6) loc=[d[0]/l,d[1]/l,d[2]/l]; }
      if(!c.u) c.u=loc||[0,0,0];
      let s=[c.u[0],c.u[1],c.u[2]]; for(const o of c.nb) if(o.disco===c.disco&&o.u&&!dying(o)){ s[0]+=o.u[0];s[1]+=o.u[1];s[2]+=o.u[2]; }
      if(loc){ s[0]+=P.pesoGradiente*loc[0];s[1]+=P.pesoGradiente*loc[1];s[2]+=P.pesoGradiente*loc[2]; }
      const l=Math.hypot(...s)||1; c.un=[s[0]/l,s[1]/l,s[2]/l]; }
    for(const c of cells) if(spinge(c)) c.u=c.un;
    // forze: spinta lungo u (solo in SPINGI e solo se non tocco l'altro disco), più planarità in tutte le fasi
    // (ognuna va alla quota media delle compagne di disco lungo u: il disco resta piatto)
    for(const c of cells){ if(!piana(c)) continue; const u=c.u; let h=0,k=0;
      const altro=c.nb.find(o=>!dying(o)&&o.disco&&o.disco!==c.disco);
      c.tocca=altro?(c.tocca||0)+1:0; if(altro&&!R.primoContatto){ R.primoContatto={a:c.x.slice(),b:altro.x.slice()}; log("primo contatto"); }
      for(const o of c.nb){ if(dying(o)) continue;
        if(o.disco===c.disco){ h+=(o.x[0]-c.x[0])*u[0]+(o.x[1]-c.x[1])*u[1]+(o.x[2]-c.x[2])*u[2]; k++; } }
      // spinge solo chi ha almeno due compagne di disco: un frammento isolato non parte da solo nel tessuto
      const m=(c.fase==="SPINGI"&&!altro&&k>=2?P.spinta:0)+(k?P.piatto*h/k:0); c.x=[c.x[0]+m*u[0],c.x[1]+m*u[1],c.x[2]+m*u[2]]; c.c=[c.c[0]+m*u[0],c.c[1]+m*u[1],c.c[2]+m*u[2]]; }
    // scava: una cellula del corpo che sta davanti a una cellula di disco (nella sua direzione di spinta) per P.mangia tick
    // muore. Il volume del canale deve sparire: spostato, scorrerebbe attorno al disco e richiuderebbe la scia (cavitazione).
    if(P.mangia) for(const c of cells){ if(dying(c)||c.disco||c.fase!=="SPINGI") continue; let davanti=false;
      for(const o of c.nb){ if(dying(o)||!o.disco||!o.u||o.fase!=="SPINGI") continue; const d=[c.x[0]-o.x[0],c.x[1]-o.x[1],c.x[2]-o.x[2]], l=Math.hypot(...d)||1;
        if((d[0]*o.u[0]+d[1]*o.u[1]+d[2]*o.u[2])/l>0.5){ davanti=true; break; } }
      c.davanti=davanti?(c.davanti||0)+1:0; if(c.davanti>=P.mangia*sc){ c.state="APOPTOSI"; R.scavate=(R.scavate||0)+1; log("scava"); } }
    // stampa(canale) ‖ svuota(canale): il disco marca il vuoto dietro di sé; una cellula di corpo che si trova nel
    // canale muore, così la bocca resta aperta invece di richiudersi. Sono le stesse primitive stampa e svuota,
    // applicate durante la deformazione invece che alla fine.
    if(P.canale) for(const c of cells){ if(dying(c)||!DOPO[c.fase]) continue;
      if(c.disco){ stampaVuoto(c.x); continue; }
      if(leggi(c.x)===2){ if((c.tCan=(c.tCan||0)+1)>=P.fuoriCanale){ c.state="APOPTOSI"; R.svuotate=(R.svuotate||0)+1; log("svuota il canale"); } } else c.tCan=0; }
    // evento di fine della spinta: un contatto che dura (un tocco di passaggio non basta). Poi la fase avanza come onda.
    // attacca: le cellule dei due dischi si incollano più forte (colla) e la spinta si ferma ovunque.
    // stacca, dopo durataAttacca tick: i dischi fusi smettono di aderire al resto e il cuscino si separa dall'anello.
    for(const c of cells){ if(dying(c)) continue;
      if(duePoli&&c.fase==="SPINGI"&&c.disco&&c.tocca>=P.tieni){ c.fase="ATTACCA"; log("contatto stabile: attacca"); }
      // gastrula: la spinta finisce quando il disco è sceso fino a ρ ≤ rhoFine·R. Poi si passa direttamente a mantieni.
      if(!duePoli&&c.fase==="SPINGI"&&c.disco&&c.rho<=P.rhoFine*S.R){ c.fase="MANTIENI"; log("fine della spinta"); }
      if(!duePoli&&c.fase==="SPINGI"&&c.nb.some(o=>o.fase==="MANTIENI")) c.fase="MANTIENI";
      if(c.fase==="SPINGI"&&c.nb.some(o=>ord[o.fase]>=2)) c.fase="ATTACCA";
      if(c.fase==="ATTACCA"){ c.tA=(c.tA||0)+1; if(c.tA>=P.durataAttacca){ c.fase="STACCA"; log("stacca"); } }
      if(c.fase==="ATTACCA"&&c.nb.some(o=>o.fase==="STACCA")) c.fase="STACCA";
      // fora: dopo lo stacco la membrana doppia (il cuscino) muore e il foro si apre, come la membrana buccofaringea
      if(c.fase==="STACCA"){ c.tS=(c.tS||0)+1; if(c.tS>=P.durataStacca){ c.fase="FORA"; log("fora"); } }
      if(c.fase==="STACCA"&&c.nb.some(o=>o.fase==="FORA")) c.fase="FORA";
      if(c.fase==="FORA"&&c.disco&&P.fora){ c.state="APOPTOSI"; } }
    // cucitura: durante la fusione una cellula del corpo che tocca entrambi i dischi è rimasta intrappolata nella linea
    // di fusione e muore (come l'epitelio della linea mediana nella fusione del palato)
    // "tra i dischi" vuol dire da parti opposte: le direzioni verso una vicina di D1 e una di D2 formano più di 120°
    for(const c of cells){ if(dying(c)||c.disco||ord[c.fase]<2) continue; let tra=false;
      for(const o of c.nb){ if(dying(o)||o.disco!==1) continue; const a=[o.x[0]-c.x[0],o.x[1]-c.x[1],o.x[2]-c.x[2]], la=Math.hypot(...a);
        for(const q of c.nb){ if(dying(q)||q.disco!==2) continue; const b=[q.x[0]-c.x[0],q.x[1]-c.x[1],q.x[2]-c.x[2]];
          if((a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/(la*Math.hypot(...b))<-0.5){ tra=true; break; } } if(tra) break; }
      if(tra){ c.state="APOPTOSI"; R.cuciture=(R.cuciture||0)+1; log("cucitura"); } }
    // mantieni(stampo): finito il foro, ogni cellula del corpo stampa la propria posizione nella matrice (una griglia di
    // occupazione che non cambia più). Da lì la forma è tenuta dalla stessa coppia di regole di sempre, lette dallo stampo:
    //   cresci: una cellula con pochi vicini mette una figlia in una direzione libera dove lo stampo è pieno;
    //   svuota: una cellula che si trova dove lo stampo è vuoto, per abbastanza tick, muore.
    // È mantieni(R) = sempre(cresci(R) ‖ svuota(¬R)) con R ricordata dalla matrice invece che calcolata da ρ.
    for(const c of cells){ if(dying(c)) continue;
      if(c.fase==="FORA"){ c.tF=(c.tF||0)+1; if(c.tF>=P.assesto){ c.fase="MANTIENI"; log("mantieni: stampo"); } }
      if(c.fase==="FORA"&&c.nb.some(o=>o.fase==="MANTIENI")) c.fase="MANTIENI";
      if(c.fase==="MANTIENI"&&!c.stampata&&!c.disco){ c.stampata=true; stampa(c.x); } }
    const nati=[];
    for(const c of cells){ if(dying(c)||c.fase!=="MANTIENI"||c.disco) continue;
      // svuota(¬stampo)
      if(leggi(c.x)!==1){ if((c.fuoriStampo=(c.fuoriStampo||0)+1)>P.fuoriStampo){ c.state="APOPTOSI"; R.morteStampo=(R.morteStampo||0)+1; continue; } } else c.fuoriStampo=0;
      // cresci(stampo)
      // cercare dove crescere costa (48 direzioni): ogni 3 tick, sfasato tra le cellule
      if(c._ph==null) c._ph=(Math.random()*3)|0; if((S.tick+c._ph)%3!==0&&c.timerT==null) continue;
      const nbl=c.nb.filter(o=>!dying(o)); if(nbl.length>=10) continue; let best=null,bm=0;
      for(const d of DIRS){ const p=[c.x[0]+d[0]*PASSO,c.x[1]+d[1]*PASSO,c.x[2]+d[2]*PASSO]; if(leggi(p)!==1) continue; let m=Infinity;
        for(const o of nbl){ const q=dist(p,o.x); if(q<m) m=q; } if(m<2*r0*0.9) continue; if(m>bm){ bm=m; best=p; } }
      if(!best){ c.timerT=null; continue; } if(c.timerT==null) c.timerT=55+Math.random()*60;
      if(--c.timerT<=0){ c.timerT=null; nati.push({x:best,c:c.c.slice(),Q:c.Q.slice(),rho:c.rho,r:r0,state:"QUIESCENTE",timer:999,migra:0,blk:0,out:0,age:0,g:0,fk:0,nb:[],md:null,
          fase:c.fase,tA:c.tA,tS:c.tS,tF:c.tF,impegnata:true,stampata:true,fuoriStampo:0,q:0,e1:{bv:-Infinity,ba:0,win:0},e2:{bv:-Infinity,ba:0,win:0},d1:Infinity,d1st:0,d2:Infinity,d2st:0,jit:Math.random(),polo:0,epi:false}); R.nascite=(R.nascite||0)+1; } }
    if(nati.length) S.cells.push(...nati);
    // adesione differenziale (letta dal motore per ogni coppia a contatto: dipende solo dalle due cellule).
    // attacca: D1–D2 incollati. stacca: i dischi non aderiscono più al corpo e si incollano forte fra loro, così il cuscino
    // si raccoglie e si ritira dalle pareti del foro
    S.adesione=(a,b)=>{ const st=ord[a.fase]>=3||ord[b.fase]>=3;
      if(a.disco&&b.disco){ if(st) return P.colla; if(a.disco!==b.disco&&ord[a.fase]>=2&&ord[b.fase]>=2) return P.colla; return P.coesioneEpi; }
      if((a.disco||b.disco)&&st) return 0;
      if(a.epi&&b.epi) return (a.giunzioni&&a.giunzioni.has(b))||(b.giunzioni&&b.giunzioni.has(a))?P.coesioneEpi:P.apicale; return 1; };
  }
  // osservatore
  function report(){ const L=S.cells.filter(c=>!dying(c)), C=E.centro(); const p1=L.filter(c=>c.polo===1), p2=L.filter(c=>c.polo===2);
    let ang=null; if(p1.length&&p2.length){ const a=[p1[0].x[0]-C[0],p1[0].x[1]-C[1],p1[0].x[2]-C[2]], b=[p2[0].x[0]-C[0],p2[0].x[1]-C[1],p2[0].x[2]-C[2]];
      ang=Math.acos((a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/(Math.hypot(...a)*Math.hypot(...b)))*180/Math.PI; }
    const fasi={}; L.forEach(c=>fasi[c.fase]=(fasi[c.fase]||0)+1);
    const D=[1,2].map(k=>{const A=L.filter(c=>c.disco===k);if(!A.length)return null;
      // spessore del disco lungo l'asse e raggio attorno all'asse
      const pa=k===1?p1:p2; const ax=p1.length&&p2.length?[p2[0].x[0]-p1[0].x[0],p2[0].x[1]-p1[0].x[1],p2[0].x[2]-p1[0].x[2]]:null; if(!ax) return {n:A.length};
      const al=Math.hypot(...ax), u=ax.map(v=>v/al); const hs=A.map(c=>(c.x[0]-C[0])*u[0]+(c.x[1]-C[1])*u[1]+(c.x[2]-C[2])*u[2]);
      const rs=A.map((c,i)=>{const w=[c.x[0]-C[0]-hs[i]*u[0],c.x[1]-C[1]-hs[i]*u[1],c.x[2]-C[2]-hs[i]*u[2]];return Math.hypot(...w);});
      return {n:A.length,hmin:Math.min(...hs),hmax:Math.max(...hs),rmax:Math.max(...rs)};});
    return {p1:p1.length,p2:p2.length,ang,fasi,D,eventi:R.eventi.slice(),r1:p1.length?dist(p1[0].x,C):null,r2:p2.length?dist(p2[0].x,C):null}; }
  return {step,report,P,R,leggi};
}
window.Primitive3D=Primitive;
})();
