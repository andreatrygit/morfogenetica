for(const [rr,pi,lab] of [[4,true,'disco 2D'],[6,false,'sfera 3D']]){
const E=Motore3D({r:rr,piano:pi});E.set({m:'disco',R:60,matOn:true});E.reset();
for(let t=0;t<(pi?2500:3000);t++)E.step();
const C=E.centro(), L=E.S.cells.filter(c=>!E.dying(c));
const bin={}; for(const c of L){const r=Math.hypot(c.x[0]-C[0],c.x[1]-C[1],c.x[2]-C[2]); const k=Math.floor(r/10)*10;
  (bin[k]=bin[k]||[]).push(c.fr||0);}
console.log(lab+' — frazione libera per fascia di raggio (mediana, max, n):');
for(const k of Object.keys(bin).map(Number).sort((a,b)=>a-b)){const v=bin[k].sort((a,b)=>a-b);
 console.log('  r '+k+'-'+(k+10)+': mediana '+v[v.length>>1].toFixed(2)+' max '+v[v.length-1].toFixed(2)+'  n '+v.length);}
}
