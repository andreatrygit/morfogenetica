const N=+(process.env.N||3);
const casi=[['disco',6,false,0,'sfera 3D',5000],['disco',4,true,0,'disco 2D',4000],
            ['anello',6,false,35,'guscio 3D',6000],['anello',4,true,35,'anello 2D',5000]];
for(const [mm,rr,pi,rin,lab,TT] of casi){
 for(let k=1;k<=N;k++){
  const E=Motore3D({r:rr,piano:pi});E.set({m:mm,R:60,Rin:rin,matOn:true});E.reset();
  let tq=null; for(let t=0;t<TT;t++){E.step(); if(tq==null&&E.S.calm>300&&E.S.cells.length>40)tq=E.S.tick;}
  const s=E.stats(), n0=s.n;
  const d=p=>p.toFixed(0);
  let r=lab+' #'+k+': n '+n0+' esterno '+d(s.rmin)+'-'+d(s.rmax)+(rin?' interno '+d(s.imin)+'-'+d(s.imax):'')
        +' err '+s.err.toFixed(2)+' quiete t='+(tq||'mai');
  // ferita: si asporta un quarto
  const tolte=E.ablate(v=>v[0]>0&&v[1]>0,true);
  let tr=null; for(let t=0;t<TT;t++){E.step(); if(tr==null&&E.S.calm>300)tr=E.S.tick;}
  const s2=E.stats();
  r+=' | ferita -'+tolte+' → n '+s2.n+' esterno '+d(s2.rmin)+'-'+d(s2.rmax)+(rin?' interno '+d(s2.imin)+'-'+d(s2.imax):'')
     +' err '+s2.err.toFixed(2)+' riparata t='+(tr||'mai');
  console.log(r);
 }}
