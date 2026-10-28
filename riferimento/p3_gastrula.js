/* Gastrulazione in 3D (stessa catena, stesso motore) sul motore unico. Niente spingi, niente stampo, niente tempi globali.
   Catena: cresci(anello) ; elegge(P1) ; campo(d1) ; marca(arco apicale) ; costringi(arco che si allarga)
   Il polo è eletto localmente: biglietto congelato = ρ, vince il più alto, con chiusura.
   L'arco si allarga mentre la piega scende — è l'anello che "scende e si allarga" forzando la ripiegatura.
   Nessuna adesione imposta fra le due labbra: se restano separate è per la polarità apicale. */
const RR=+(process.env.RR||6), RIN=+(process.env.RIN||35), RE=+(process.env.RE||60);
const KC=+(process.env.KC||0.4), AP0=+(process.env.AP0||25), CRESCE=+(process.env.CRESCE||1.2);
const ALONE=+(process.env.ALONE||2), TMAX=+(process.env.TMAX||9000), NP=+(process.env.NP||2);
const K=window.MORFO_COSTANTI;

for(let rep=1;rep<=NP;rep++){
  const E=Motore3D({r:RR}); E.set({m:"anello",R:RE,Rin:RIN,matOn:true}); E.reset();
  const S=E.S, live=()=>S.cells.filter(c=>!E.dying(c)), dist=(a,b)=>Math.hypot(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
  let fase="CRESCI", tFase=0, eventi=[], ax=null, apFine=0, stD=0, poloD=null, succ=0;
  const log=s=>eventi.push(s+" t="+S.tick);

  // --- misura: fin dove arriva un canale vuoto dall'esterno lungo l'asse del polo
  const canale=()=>{ const C=E.centro(), L=live(); let apert=0;
    for(let s=RE+8;s>-RE;s-=1){ const p=[C[0]+ax[0]*s,C[1]+ax[1]*s,C[2]+ax[2]*s];
      let m=Infinity; for(const c of L){const d=dist(c.x,p); if(d<m)m=d;}
      if(m<RR*1.1) break; apert=s; }
    return apert; };            // ρ a cui il canale si chiude: più è basso, più la piega è profonda
  const pezzi=()=>{ const L=live(), vis=new Set(); let k=0;
    for(const c of L){ if(vis.has(c))continue; k++; const st=[c]; vis.add(c);
      while(st.length){const a=st.pop(); for(const o of a.nb) if(!E.dying(o)&&!vis.has(o)){vis.add(o);st.push(o);}}}
    return k; };
  const misura=nome=>{ const L=live(), M=L.filter(c=>c.con>0), s=E.stats();
    const rm=M.length?Math.min(...M.map(c=>c.rho)):NaN;
    console.log(`  ${nome} t=${S.tick}: n ${L.length}, arco ${M.length} cellule, fronte a ρ=${rm.toFixed(0)}`
      +`, canale aperto fino a ρ=${canale().toFixed(0)}, bordo ${s.rmin.toFixed(0)}–${s.rmax.toFixed(0)}`
      +`, coppie ${S.cos.n}, dominio ${L.filter(c=>c.era).length}/${L.length}, pezzi ${pezzi()}`); };

  for(let t=0;t<TMAX;t++){
    E.step();
    const L=live();
    // --- elezione di P1: biglietto congelato = ρ. Chi ha ricevuto il vincitore smette di candidarsi.
    if(fase==="CRESCI"){ if(S.calm>200&&L.length>40){ fase="ELEGGI1"; tFase=S.tick; log("quiete, si elegge"); } }
    else if(fase==="ELEGGI1"){
      // chiusura ereditata: una cellula nata dove l'elezione è già decisa non si candida. Senza questo
      // il ricambio immette biglietti nuovi che battono il vincitore e nessuno arriva a 45 tick.
      for(const c of L){ if(c.bv==null){ const vic=c.nb.filter(o=>!E.dying(o));
        // senza vicine una cellula non sa nulla dell'elezione: non si candida. E se una vicina è già chiusa,
        // qui l'elezione è decisa. Senza questo i nati nuovi battono il vincitore e nessuno arriva a 45.
        const dentro=vic.length===0||vic.some(o=>o.chiuso);
        c.bv=dentro?-Infinity:c.rho; c.ba=0; c.win=0; c.chiuso=dentro; } }
      for(const c of L){ let bv=c.chiuso?-Infinity:c.bv, ba=c.ba;
        for(const o of c.nb){ if(E.dying(o)||o.bv==null) continue; if(o.ba+1>K.scadenza) continue;
          if(o.bv>bv){bv=o.bv;ba=o.ba+1;} }
        c.nbv=bv; c.nba=ba; }
      for(const c of L){ const mio=c.nbv===c.bv&&c.nba===0; c.win=mio?c.win+1:0; c.bv=c.nbv; c.ba=c.nba;
        if(!mio) c.chiuso=true; }
      const p1=L.find(c=>c.win>=45);
      if(p1){ p1.polo=1; fase="CAMPO"; tFase=S.tick; log("eletto P1"); }
      else if(S.tick-tFase>6000){ log("elezione non conclude, candidati aperti: "+L.filter(c=>!c.chiuso).length+", max win "+Math.max(...L.map(c=>c.win||0))); break; }
    }
    // --- campo d1 dal polo, e asse per la sola misura (le cellule non lo usano)
    else {
      // successione del polo: una marca non deve morire con la cellula che la porta. Se il polo muore la
      // passa alla vicina viva più prossima. Senza questo il campo d1 perde l'ancora e l'arco svanisce.
      for(const c of S.cells){ if(c.polo!==1||!E.dying(c)) continue;
        let best=null,bd=Infinity;
        for(const o of c.nb){ if(E.dying(o)||o.polo) continue; const q=dist(o.x,c.x); if(q<bd){bd=q;best=o;} }
        c.polo=0; if(best){ best.polo=1; best.d1=0; succ++; } }
      for(const c of L){ if(c.d1==null)c.d1=Infinity; }
      for(const c of L){ let d=c.polo===1?0:Infinity;
        for(const o of c.nb){ if(E.dying(o)||o.d1==null) continue; const v=o.d1+dist(o.x,c.x); if(v<d)d=v; }
        c.nd1=Math.min(d,9999); }
      let stab=true; for(const c of L){ if(Math.abs(c.nd1-c.d1)>0.5) stab=false; c.d1=c.nd1; }
      const P1=L.find(c=>c.polo===1);
      if(P1){ const C=E.centro(), a=[P1.x[0]-C[0],P1.x[1]-C[1],P1.x[2]-C[2]], l=Math.hypot(...a); if(l>2) ax=[a[0]/l,a[1]/l,a[2]/l]; }
      if(fase==="CAMPO"&&stab&&S.tick-tFase>60){ fase="COSTRINGI"; tFase=S.tick; log("campo stabile, si costringe"); misura("prima"); }
    }
    // --- marca + costringi: l'arco apicale entro un raggio che si allarga nel tempo
    if(fase==="COSTRINGI"){
      const ap=AP0+CRESCE*(S.tick-tFase)/100;      // raggio dell'arco, in unità di d1
      for(const c of S.cells){ if(E.dying(c)){c.con=0;continue;}
        // sulla superficie esterna: ha una faccia libera che guarda via dal centro virtuale. Due cose
        // che la cellula ha già — lo spazio vuoto attorno a sé, e il vettore al centro.
        const u=c.lib; const apicale=!!u&&(u[0]*(c.x[0]-c.c[0])+u[1]*(c.x[1]-c.c[1])+u[2]*(c.x[2]-c.c[2]))>0;
        c.con=(c.d1<ap&&apicale)?KC:0; if(c.con>0)c.era=true; }
      // rallenta: arresto mitotico nell'arco e in un alone di ALONE contatti
      for(const c of S.cells) c.i0=c.con>0;
      for(let g=0;g<ALONE;g++){ for(const c of S.cells) if(!c.i0&&c.nb.some(o=>o.i0)) c.i1=true;
        for(const c of S.cells){ if(c.i1){c.i0=true;c.i1=false;} } }
      for(const c of S.cells) c.ferma=!!c.i0;
      const dt=S.tick-tFase;
      if(dt===300||dt===1500||dt===3000||dt===5000) misura("dopo "+dt);
      if(dt>=(+(process.env.TC||3000))){ misura("fine costrizione");
        // rilascio: la costrizione smette e il tessuto torna a fare manutenzione. Se la coppa regge,
        // regge da sola: nessuna adesione imposta fra le labbra, nessuno stampo che la ricordi.
        apFine=ap; fase=process.env.RILASCIA?"RILASCIO":"MANTIENI"; tFase=S.tick; log(fase.toLowerCase());
        if(process.env.RILASCIA){
          if(process.env.GIUNZIONI){ const n=E.consolida(); log("consolidate "+n+" giunzioni"); }
          for(const c of S.cells){ c.con=0;
            // DOMINIO: le invaginate sono una sottopopolazione differenziata e post-mitotica. La loro regola
            // non è più la banda vecchia: è "resta come sei". Identità intrinseca, locale, permanente.
            c.ferma=!!(process.env.DOMINIO&&c.era);
            // il dominio passa subito sotto la propria regola: anche un solo tick sotto la banda madre
            // basta a far partire l'apoptosi di chi la piega ha portato sotto Rin, e l'apoptosi non torna.
            if(process.env.CAMPO&&c.era) c.impegnata=true;
            if(process.env.FREEZE) c.impegnata=true; } } }
    }
    else if(fase==="RILASCIO"){
      if(process.env.DOMINIO) for(const c of S.cells) c.ferma=!!c.era;
      const dt=S.tick-tFase;
      if(dt===200||dt===800||dt===2000||dt===4000) misura("rilasciata da "+dt);
      // ferita nella coppa: metà del dominio asportata. Con la sola regola "resta come sei" il dominio
      // non può ricrescere — è il prezzo di non avere un campo posizionale proprio.
      if(process.env.FERISCI&&dt===1000){ const C=E.centro(); let k=0;
        for(const c of S.cells){ if(E.dying(c)||!c.era) continue;
          const d=[c.x[0]-C[0],c.x[1]-C[1],c.x[2]-C[2]], h=d[0]*ax[0]+d[1]*ax[1]+d[2]*ax[2];
          const t=[d[0]-h*ax[0],d[1]-h*ax[1],d[2]-h*ax[2]];
          if(t[0]+t[1]+t[2]>0){ c.state="APOPTOSI"; c.era=false; k++; } }
        S.calm=0; log("ferita nella coppa, -"+k); misura("appena ferita"); }
      if(dt>=4000) break;
    }
    else if(fase==="MANTIENI"){
      // mantieni(costringi): l'arco non si allarga più, ma chi costringe continua. Nessuno stampo:
      // la forma è tenuta dalla stessa regola locale che l'ha fatta, lasciata accesa.
      // la regola resta accesa, non la marca: una cellula che non è più sulla superficie libera smette di
      // costringere. Tenere la marca farebbe costringere fra loro anche le due pareti del solco, e lacera.
      for(const c of S.cells){ if(E.dying(c)){c.con=0;continue;}
        const u=c.lib; const apicale=!!u&&(u[0]*(c.x[0]-c.c[0])+u[1]*(c.x[1]-c.c[1])+u[2]*(c.x[2]-c.c[2]))>0;
        c.con=(c.d1<apFine&&apicale)?KC:0; }
      for(const c of S.cells) c.i0=c.con>0;
      for(let g=0;g<ALONE;g++){ for(const c of S.cells) if(!c.i0&&c.nb.some(o=>o.i0)) c.i1=true;
        for(const c of S.cells){ if(c.i1){c.i0=true;c.i1=false;} } }
      for(const c of S.cells) c.ferma=!!c.i0;
      const dt=S.tick-tFase;
      if(dt===500||dt===2000||dt===4000) misura("mantenuta da "+dt);
      if(dt>=4000) break;
    }
  }
  console.log(`r=${RR} #${rep} | ${eventi.join(", ")} | successioni del polo: ${succ}`);
  if(process.env.DUMP&&ax){ const C=E.centro();
    require('fs').writeFileSync(process.env.DUMP+'_'+rep+'.json',JSON.stringify(live().map(c=>{
      const d=[c.x[0]-C[0],c.x[1]-C[1],c.x[2]-C[2]], h=d[0]*ax[0]+d[1]*ax[1]+d[2]*ax[2];
      const t=Math.hypot(d[0]-h*ax[0],d[1]-h*ax[1],d[2]-h*ax[2]);
      return [h, t, c.era?1:0];}))); }
}
