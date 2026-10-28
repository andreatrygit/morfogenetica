/* Simulatore: motore + primitive, fuori dalla grafica. Gira in un Web Worker (così il calcolo non ruba tempo al
   disegno) oppure, se il browser non lo permette, nella pagina stessa con lo stesso protocollo di messaggi.
   La pagina manda comandi (init, set, reset, cellule, corsa, ferita); il simulatore manda istantanee:
   posizioni e raggi in un Float32Array, stato e marche in due Uint8Array, ρ in un Float32Array, più le statistiche. */
(function(){
const STATI=["MIGRA","PROLIFERA","QUIESCENTE","BORDO","APOPTOSI","IMPALCATURA","CAVITA"], SI={}; STATI.forEach((s,i)=>SI[s]=i);

function Simulatore(invia){
  let E=null, Pr=null, prog=false, corre=true, velocita=2, par={r:6,m:"disco",R:60,Rin:35,rot:0,matOn:true}, quota=0, ultimo=0, conta=0, fermo=false;
  // prog: false, oppure il nome del programma ("toro", "gastrula")
  function nuovo(scalda){ E=window.Motore3D({r:par.r}); E.set({m:par.m,R:par.R,Rin:par.Rin,rot:par.rot,matOn:par.matOn}); E.reset(); Pr=window.Primitive3D(E);
    if(prog) Pr.P.programma=prog;
    if(scalda) for(let i=0;i<200;i++) E.step(); }
  function istantanea(){
    const cells=E.S.cells, n=cells.length, pos=new Float32Array(4*n), st=new Uint8Array(n), fl=new Uint8Array(n), rho=new Float32Array(n);
    // durante la metamorfosi il centro stimato non ha più senso: si inquadra il baricentro
    let C=E.centro(), impegnata=false; for(const c of cells) if(c.impegnata){ impegnata=true; break; }
    if(prog&&impegnata){ let s=[0,0,0],k=0; for(const c of cells){ if(E.dying(c)) continue; s[0]+=c.x[0];s[1]+=c.x[1];s[2]+=c.x[2];k++; } if(k) C=s.map(v=>v/k); }
    let a=[0,0,0],na=0;
    for(let i=0;i<n;i++){ const c=cells[i]; pos[4*i]=c.x[0]-C[0]; pos[4*i+1]=c.x[1]-C[1]; pos[4*i+2]=c.x[2]-C[2]; pos[4*i+3]=c.r;
      st[i]=SI[c.state]??2; fl[i]=(c.polo?c.polo:0)|((c.disco||0)<<2)|(c.epi?16:0)|(E.dying(c)?32:0); rho[i]=c.rho;
      // asse del toro (osservatore): direzione di spinta media dei dischi, u per il disco 1 e −u per il disco 2
      if(prog==="toro"&&c.disco&&c.u&&!E.dying(c)){ const sg=c.disco===1?1:-1; a[0]+=sg*c.u[0];a[1]+=sg*c.u[1];a[2]+=sg*c.u[2];na++; } }
    const m={tipo:"istantanea",n,pos,st,fl,rho,tick:E.S.tick,calm:E.S.calm,R:E.S.R};
    const la=Math.hypot(...a); if(na>=6&&la>na*0.5) m.asse=a.map(v=>v/la);
    if(conta++%6===0){ const s=E.stats(); m.stats={n:s.n,rmin:s.rmin,rmax:s.rmax,imin:s.imin,imax:s.imax,err:s.err,cnt:s.cnt,readers:s.readers,matOn:E.S.matOn};
      if(prog){ const f={}; for(const c of cells) if(c.fase) f[c.fase]=(f[c.fase]||0)+1; m.fasi=f; m.eventi=Pr.R.eventi.slice(); m.forato=Pr.R.tick0&&Pr.R.tick0.fora!=null; } }
    invia(m,[pos.buffer,st.buffer,fl.buffer,rho.buffer]); }
  function ciclo(){
    const t0=Date.now();
    if(corre&&E){ quota=Math.min(quota+velocita,velocita*4); while(quota>=1&&Date.now()-t0<30){ E.step(); if(prog) Pr.step(); quota--; } }
    if(E&&Date.now()-ultimo>=33){ ultimo=Date.now(); istantanea(); }
    if(!fermo) setTimeout(ciclo, Math.max(0,16-(Date.now()-t0)));
  }
  function ferita(d){ if(!E) return; const R=E.S.R;
    if(d.tipo==="calotta"){ const u=d.u, l=Math.hypot(...u); E.ablate(x=>(x[0]*u[0]+x[1]*u[1]+x[2]*u[2])/l>R*0.45,d.cancella); }
    else if(d.tipo==="centro") E.ablate(x=>Math.hypot(...x)<R*0.45,d.cancella);
    else if(d.tipo==="taglio") E.dividi(Math.round(Math.min(R*1.08,E.W-R-8)),d.cancella); }
  function ricevi(d){
    switch(d.cmd){
      case "init": Object.assign(par,d.par||{}); prog=d.prog||false; velocita=d.velocita||velocita; nuovo(true); ultimo=0; ciclo(); break;
      case "forma": par.m=d.m; prog=d.prog||false; nuovo(false); break;
      case "cellule": par.r=d.r; nuovo(false); break;
      case "reset": nuovo(false); break;
      case "set": Object.assign(par,d.par); if(E) E.set(d.par); break;
      case "corsa": if("corre" in d) corre=d.corre; if("velocita" in d) velocita=d.velocita; break;
      case "ferita": ferita(d); break;
    } }
  return ricevi;
}

if(typeof WorkerGlobalScope!=="undefined"&&self instanceof WorkerGlobalScope){
  self.window=self; importScripts("macchine.js","motore3d.js","primitive3d.js");
  const r=Simulatore((m,t)=>self.postMessage(m,t)); self.onmessage=e=>r(e.data);
} else window.Simulatore=Simulatore;
})();
