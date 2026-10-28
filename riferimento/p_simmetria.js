/* Prova della legge di simmetria: elegge(biglietto=f) rompe la simmetria se e solo se f è
   costante sulle orbite del gruppo attuale.

   Previsione. In sequenza, P2 ha biglietto d₁, che varia sull'orbita e ha massimo determinato
   nell'antipodo: l'elezione è equivariante e l'angolo fra i due poli deve essere ~180°.
   In parallelo, entrambi i biglietti sono ρ, costante sull'orbita: decide il rumore, e i poli
   non hanno motivo di essere antipodali.

   Tre modi:
     SEQ  P1 con biglietto ρ, poi campo d₁, poi P2 con biglietto d₁
     PARU P1 e P2 in parallelo, biglietto identico ρ per entrambe (lettura letterale del ‖)
     PARN P1 e P2 in parallelo, biglietto ρ più un errore di misura indipendente per elezione
          (lettura più generosa: due processi che misurano la stessa cosa, ognuno col suo errore)
   Misura: l'angolo fra i vettori dal centro virtuale ai due poli. */
const RR=+(process.env.RR||4), PIANO=process.env.D3?false:true, NP=+(process.env.NP||8);
const RUM=+(process.env.RUM||1.5);          // errore di misura di ρ, dell'ordine di quello misurato
const MODI=(process.env.MODI||"SEQ,PARU,PARN").split(",");
const K=window.MORFO_COSTANTI;
const gauss=()=>{let u=0,v=0;while(!u)u=Math.random();while(!v)v=Math.random();
  return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);};

// un'elezione generica: biglietto congelato, massimo che viaggia un contatto per passo con
// un'età che scade, chiusura, e vittoria dopo 45 passi imbattuto. Identica nei tre modi.
function elezione(E,L,k,biglietto){
  const S=E.S;
  for(const c of L){ if(c[k]==null){ const vic=c.nb.filter(o=>!E.dying(o));
    const dentro=vic.length===0||vic.some(o=>o[k]&&o[k].chiuso);
    c[k]={bv:dentro?-Infinity:biglietto(c), ba:0, win:0, chiuso:dentro}; } }
  for(const c of L){ const e=c[k]; let bv=e.chiuso?-Infinity:e.bv, ba=e.ba;
    for(const o of c.nb){ const f=o[k]; if(E.dying(o)||!f) continue; if(f.ba+1>K.scadenza) continue;
      if(f.bv>bv){ bv=f.bv; ba=f.ba+1; } }
    e.nbv=bv; e.nba=ba; }
  for(const c of L){ const e=c[k]; const mio=e.nbv===e.bv&&e.nba===0;
    e.win=mio?e.win+1:0; e.bv=e.nbv; e.ba=e.nba; if(!mio) e.chiuso=true; }
  return L.find(c=>c[k].win>=45)||null;
}

const risultati={};
for(const modo of MODI){
  const ang=[], note=[], tcs=[];
  for(let rep=1;rep<=NP;rep++){
    const E=Motore3D({r:RR,piano:PIANO}); E.set({m:"disco",R:60,matOn:true}); E.reset();
    const S=E.S, live=()=>S.cells.filter(c=>!E.dying(c));
    const dist=(a,b)=>Math.hypot(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
    let P1=null,P2=null,fase="CRESCI",tF=0,stab=0,dich=0,tconcl=null;
    for(let t=0;t<30000;t++){
      E.step(); const L=live();
      if(fase==="CRESCI"){ if(S.calm>200&&L.length>40){ fase="ELEGGI"; tF=S.tick; } continue; }
      if(modo==="SEQ"||modo==="SEQL"){
        if(!P1){ P1=elezione(E,L,"e1",c=>c.rho); if(P1) P1.polo=1; continue; }
        // campo d₁ dal polo, fino a stabilità
        for(const c of L) if(c.d1==null) c.d1=Infinity;
        for(const c of L){ let d=c.polo===1?0:Infinity;
          for(const o of c.nb){ if(E.dying(o)||o.d1==null) continue; const v=o.d1+dist(o.x,c.x); if(v<d) d=v; }
          c.nd1=Math.min(d,9999); }
        let ok=true; for(const c of L){ if(Math.abs(c.nd1-c.d1)>0.5) ok=false; c.d1=c.nd1; }
        stab=ok?stab+1:0;
        if(stab<60) continue;
        if(!P2){
          if(modo==="SEQL"){
            // P2 senza apparato: nessun biglietto congelato, nessuna età, nessuna chiusura, nessuna
            // attesa di 45 passi. Solo un predicato locale: sono un massimo locale di d₁.
            const cand=L.filter(c=>c.d1<9999&&c.polo!==1
              &&!c.nb.some(o=>!E.dying(o)&&o.d1!=null&&o.d1>c.d1));
            if(cand.length){ P2=cand[0]; P2.polo=2; dich=cand.length; }
          } else { P2=elezione(E,L,"e2",c=>(c.d1<9999?c.d1:-Infinity)); if(P2) P2.polo=2; }
        }
      } else {
        // le due elezioni corrono insieme, sulla stessa configurazione
        const b1=modo==="PARU"?(c=>c.rho):(c=>c.rho+RUM*gauss());
        const b2=modo==="PARU"?(c=>c.rho):(c=>c.rho+RUM*gauss());
        if(!P1){ P1=elezione(E,L,"e1",b1); if(P1) P1.polo=1; }
        if(!P2){ P2=elezione(E,L,"e2",b2); if(P2) P2.polo=2; }
      }
      if(P1&&P2){ tconcl=S.tick-tF; break; }
      if(S.tick-tF>20000) break;
    }
    if(!P1||!P2){ note.push("#"+rep+" non conclude"); continue; }
    if(P1===P2){ ang.push(0); note.push("#"+rep+" stessa cellula"); continue; }
    const C=E.centro();
    const u=[P1.x[0]-C[0],P1.x[1]-C[1],P1.x[2]-C[2]], v=[P2.x[0]-C[0],P2.x[1]-C[1],P2.x[2]-C[2]];
    const lu=Math.hypot(...u), lv=Math.hypot(...v);
    const cos=Math.max(-1,Math.min(1,(u[0]*v[0]+u[1]*v[1]+u[2]*v[2])/(lu*lv)));
    if(tconcl!=null)tcs.push(tconcl);
    ang.push(Math.acos(cos)*180/Math.PI); if(dich>1) note.push('#'+rep+' '+dich+' candidati');
  }
  risultati[modo]=ang; if(tcs.length) console.log('      passi per concludere: '+tcs.join(' '));
  const m=ang.length?ang.reduce((a,b)=>a+b,0)/ang.length:NaN;
  const s=ang.length>1?Math.sqrt(ang.reduce((a,b)=>a+(b-m)**2,0)/(ang.length-1)):NaN;
  console.log(modo.padEnd(5)+' angolo fra i poli: '+ang.map(a=>a.toFixed(0)+'°').join(' ')
    +'  | media '+m.toFixed(0)+'° dev '+(isFinite(s)?s.toFixed(0):'-')+'°'
    +(note.length?'  ['+note.join('; ')+']':''));
}
