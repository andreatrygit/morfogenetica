const K=window.MORFO_COSTANTI;
const varianti=[['tutte le regole',()=>{}],
 ['senza polarità (niente apicale)',()=>{K.polarita=0;}],
 ['senza divisione planare',()=>{K.planare=9;}],
 ['senza deposizione basale',()=>{K.basale=0;}],
 ['senza adesione apicale debole',()=>{K.adesioneApicale=1;}]];
const base=JSON.parse(JSON.stringify(K));
for(const [nome,f] of varianti){ Object.assign(K,JSON.parse(JSON.stringify(base))); f();
 const out=[]; for(let k=0;k<3;k++){
  const E=Motore3D({r:4,piano:true});E.set({m:'anello',R:60,Rin:35,matOn:true});E.reset();
  let tq=null; for(let t=0;t<5000;t++){E.step(); if(tq==null&&E.S.calm>300&&E.S.cells.length>40)tq=E.S.tick;}
  const s=E.stats(), P=E.S.cells.filter(c=>!E.dying(c)&&c.apic).length;
  out.push('n '+s.n+' est '+s.rmin.toFixed(0)+'-'+s.rmax.toFixed(0)+' int '+s.imin.toFixed(0)+'-'+s.imax.toFixed(0)+' err '+s.err.toFixed(2)+' quiete '+(tq||'mai')+' pol '+P); }
 console.log(nome.padEnd(34)+' | '+out.join('  |  ')); }
