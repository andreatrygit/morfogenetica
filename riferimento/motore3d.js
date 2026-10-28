/* Motore morfogenetico 3D. Usa le stesse macchine a stati di macchine.js (window.MACCHINE).
   Quattro canali, come in 2D:
   1 dentro la cellula: vettore verso il centro + orientamento del corpo (matrice Q), che ruota a caso
   2 tra cellule a contatto: ρ per differenze lungo la direzione del centro; senza matrice, media dei vettori
     letti con la propria bussola (errore = rotazione relativa dei due corpi)
   3 fuori, nella matrice: griglia volumetrica. Ogni voxel tiene la normale al piano delle fibre (come tensore
     n·nᵀ, perché n e −n sono la stessa fibra) e il ρ di chi l'ha deposta
   4 misura assoluta: la curvatura della superficie delle fibre. Spostandosi di s lungo il piano, la normale
     gira di s/ρ. Media sulla zona, correzione lenta, apoptosi per chi risulta fuori. */
(function(){
const M=window.MACCHINE, K=window.MORFO_COSTANTI, MX=K.matrice, CO=K.consenso;
const W=150;
const rnd=()=>Math.random()*2-1;
let UID=0;
// --- piccola algebra
const dist2=(a,b)=>{const x=a[0]-b[0],y=a[1]-b[1],z=a[2]-b[2];return Math.sqrt(x*x+y*y+z*z);};
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2], len=a=>Math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
const sub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]], add=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]], mul=(a,k)=>[a[0]*k,a[1]*k,a[2]*k];
const norm=a=>{const l=len(a)||1;return [a[0]/l,a[1]/l,a[2]/l];};
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const mv=(m,v)=>[m[0]*v[0]+m[1]*v[1]+m[2]*v[2],m[3]*v[0]+m[4]*v[1]+m[5]*v[2],m[6]*v[0]+m[7]*v[1]+m[8]*v[2]];
const mtv=(m,v)=>[m[0]*v[0]+m[3]*v[1]+m[6]*v[2],m[1]*v[0]+m[4]*v[1]+m[7]*v[2],m[2]*v[0]+m[5]*v[1]+m[8]*v[2]];
const mm=(a,b)=>{const r=new Array(9);for(let i=0;i<3;i++)for(let j=0;j<3;j++)r[3*i+j]=a[3*i]*b[j]+a[3*i+1]*b[3+j]+a[3*i+2]*b[6+j];return r;};
function rotM(u,t){const c=Math.cos(t),s=Math.sin(t),k=1-c,[x,y,z]=u;
  return [c+x*x*k, x*y*k-z*s, x*z*k+y*s, y*x*k+z*s, c+y*y*k, y*z*k-x*s, z*x*k-y*s, z*y*k+x*s, c+z*z*k];}
const randUnit=()=>{let v;do{v=[rnd(),rnd(),rnd()];}while(len(v)>1||len(v)<0.1);return norm(v);};
const tangent=n=>{const a=Math.abs(n[0])<0.9?[1,0,0]:[0,1,0];return norm(cross(n,a));};
// direzioni di prova per le figlie (spirale di Fibonacci)
const DIRS=[];{const N=48,g=Math.PI*(3-Math.sqrt(5));for(let i=0;i<N;i++){const y=1-2*(i+.5)/N,r=Math.sqrt(1-y*y);DIRS.push([Math.cos(g*i)*r,y,Math.sin(g*i)*r]);}}
// nel piano: stesse regole, direzioni sul cerchio. Il mondo è la sezione z = 0 dello stesso motore.
const DIRS2=[];{const N=32;for(let i=0;i<N;i++){const a=2*Math.PI*i/N;DIRS2.push([Math.cos(a),Math.sin(a),0]);}}
// direzioni di campionamento della matrice per la polarità: 8 nel piano, 14 nello spazio (assi + diagonali)
const POL3=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],
  [.577,.577,.577],[-.577,.577,.577],[.577,-.577,.577],[-.577,-.577,.577],
  [.577,.577,-.577],[-.577,.577,-.577],[.577,-.577,-.577],[-.577,-.577,-.577]];
const POL2=[[1,0,0],[.707,.707,0],[0,1,0],[-.707,.707,0],[-1,0,0],[-.707,-.707,0],[0,-1,0],[.707,-.707,0]];

function Motore(o){
  // risoluzione: il raggio della cellula è un parametro. Le grandezze locali (passo di lettura della matrice, caselle
  // della matrice, soglie vicino al centro, rumore) scalano con lui; R e le forme restano in unità assolute.
  o=o||{}; const r0=o.r||K.r, sc=r0/K.r, STEP=2*r0*0.85, VS=5*sc, VN=Math.ceil(2*W/VS), RL=(MX.raggioLettura3d||K.r)*sc;
  const piano=!!o.piano, DD=piano?DIRS2:DIRS, POL=piano?POL2:POL3;
  const LP=MX.passo*sc, LVC=MX.vicinoCentro*sc, LLD=MX.leggiDa*sc, LCR=MX.curvRhoMin*sc, LFM=MX.fuoriMargine*sc;
  let grid=null, nxt=new Int32Array(1024), PX=new Float64Array(1024), PY=new Float64Array(1024), PZ=new Float64Array(1024), PR=new Float64Array(1024), PD=new Uint8Array(1024);
  // costrizione: asse e intensità per cellula, copiati in array numerici a ogni tick
  let NX=new Float64Array(1024), NY=new Float64Array(1024), NZ=new Float64Array(1024), KC=new Float64Array(1024);
  // faccia apicale: direzione della superficie libera che non ha matrice oltre di sé
  let AX=new Float64Array(1024), AY=new Float64Array(1024), AZ=new Float64Array(1024), AP=new Uint8Array(1024);
  const S={m:M.disco,R:60,Rin:35,rot:0,noise:0.02,matOn:true,cells:[],tick:0,calm:0,killed:0};
  // matrice: densità, tensore delle normali (6 componenti), ρ
  const ND=VN*VN*VN, D=new Float32Array(ND), T=new Float32Array(ND*6), MR=new Float32Array(ND);
  const vi=x=>Math.floor((x+W)/VS);
  const vox=(i,j,k)=>(k*VN+j)*VN+i;
  const dying=c=>S.m.muoiono.includes(c.state);
  function forVox(p,rad,fn){const px=p[0],py=p[1],pz=p[2],r2=rad*rad,i0=Math.max(0,vi(px-rad)),i1=Math.min(VN-1,vi(px+rad)),j0=Math.max(0,vi(py-rad)),j1=Math.min(VN-1,vi(py+rad)),k0=Math.max(0,vi(pz-rad)),k1=Math.min(VN-1,vi(pz+rad));
    for(let k=k0;k<=k1;k++){const qz=-W+(k+.5)*VS-pz;for(let j=j0;j<=j1;j++){const qy=-W+(j+.5)*VS-py;for(let i=i0;i<=i1;i++){const qx=-W+(i+.5)*VS-px;if(qx*qx+qy*qy+qz*qz<=r2)fn((k*VN+j)*VN+i,qx,qy,qz);}}}}
  function mxClear(p,rad){ if(rad<0){D.fill(0);T.fill(0);MR.fill(0);return;} forVox(p,rad,v=>{D[v]=0;MR[v]=0;for(let a=0;a<6;a++)T[6*v+a]=0;}); }
  function mxDeposit(p,n,rho){ const o0=n[0]*n[0],o1=n[1]*n[1],o2=n[2]*n[2],o3=n[0]*n[1],o4=n[0]*n[2],o5=n[1]*n[2];
    forVox(p,r0,(v,qx,qy,qz)=>{const dd=D[v],a=MX.deposito*(1-dd)+MX.rimodella,b=1-a,t=6*v;
      T[t]=T[t]*b+a*o0;T[t+1]=T[t+1]*b+a*o1;T[t+2]=T[t+2]*b+a*o2;T[t+3]=T[t+3]*b+a*o3;T[t+4]=T[t+4]*b+a*o4;T[t+5]=T[t+5]*b+a*o5;
      const rq=rho+qx*n[0]+qy*n[1]+qz*n[2]; MR[v]=dd<1e-3?rq:MR[v]*b+a*rq; D[v]=dd+MX.deposito*(1-dd);}); }
  // lettura: normale dominante (iterazione di potenza a partire da una direzione), coerenza, ρ di memoria
  function mxRead(p,guess){ let sd=0,n=0,sr=0,A0=0,A1=0,A2=0,A3=0,A4=0,A5=0; const px=p[0],py=p[1],pz=p[2],rad=RL,r2=rad*rad;
    const i0=Math.max(0,vi(px-rad)),i1=Math.min(VN-1,vi(px+rad)),j0=Math.max(0,vi(py-rad)),j1=Math.min(VN-1,vi(py+rad)),k0=Math.max(0,vi(pz-rad)),k1=Math.min(VN-1,vi(pz+rad));
    for(let k=k0;k<=k1;k++){const qz=-W+(k+.5)*VS-pz;for(let j=j0;j<=j1;j++){const qy=-W+(j+.5)*VS-py, base=(k*VN+j)*VN;for(let i=i0;i<=i1;i++){const qx=-W+(i+.5)*VS-px;
      if(qx*qx+qy*qy+qz*qz>r2) continue; n++; const v=base+i, d=D[v]; if(d<=0) continue; sd+=d; sr+=d*MR[v]; const t=6*v;
      A0+=d*T[t];A1+=d*T[t+1];A2+=d*T[t+2];A3+=d*T[t+3];A4+=d*T[t+4];A5+=d*T[t+5];}}}
    const A=[A0,A1,A2,A3,A4,A5];
    if(sd<1e-6) return {dens:0}; for(let t=0;t<6;t++)A[t]/=sd;
    // direzione dominante: iterazione di potenza sul tensore, conti scalari
    const a00=A[0],a11=A[1],a22=A[2],a01=A[3],a02=A[4],a12=A[5];
    let ex=1,ey=0,ez=0; if(guess){const gl=Math.sqrt(guess[0]*guess[0]+guess[1]*guess[1]+guess[2]*guess[2]); if(gl>1e-6){ex=guess[0]/gl;ey=guess[1]/gl;ez=guess[2]/gl;}}
    for(let it=0;it<8;it++){ const x=a00*ex+a01*ey+a02*ez, y=a01*ex+a11*ey+a12*ez, z=a02*ex+a12*ey+a22*ez, l=Math.sqrt(x*x+y*y+z*z)||1; ex=x/l;ey=y/l;ez=z/l; }
    const lam=ex*(a00*ex+a01*ey+a02*ez)+ey*(a01*ex+a11*ey+a12*ez)+ez*(a02*ex+a12*ey+a22*ez), tr=a00+a11+a22;
    return {dens:sd/n, coh:tr>1e-9?(3*lam/tr-1)/2:0, n:[ex,ey,ez], rho:sr/sd}; } // coh = 1 se tutte le fibre hanno la stessa normale

  function newCell(p,c,Q,rho){return {id:++UID,x:p,c,Q,rho,r:r0,ev:-1,pausa:0,prs:0,anc:0,cmp:0,stc:0,state:"MIGRA",timer:55+Math.random()*60,migra:30,blk:0,out:0,age:0,g:0,fk:0,nb:[],md:null};}
  function reset(){ S.tick=0;S.calm=0;S.killed=0; mxClear(null,-1);
    UID=0; S.cells=[{id:++UID,x:[0,0,0],c:[0,0,0],Q:[1,0,0,0,1,0,0,0,1],rho:0,r:r0,ev:-1,pausa:0,prs:0,anc:0,cmp:0,stc:0,state:"PROLIFERA",timer:20,migra:0,blk:0,out:0,age:0,g:0,fk:0,nb:[],md:null}]; }

  function step(){
    const cells=S.cells, n=cells.length, R=S.R;
    // la banda di crescita: [Rin, R] per le forme che dichiarano un lume imposto, [0, R] per quelle il cui
    // lume emerge dalla profondità. Rin resta il parametro, ma diventa lo spessore della parete.
    // la banda di crescita: [Rin, R]. I due bordi sono la stessa regola applicata due volte.
    const rin=S.m.parametri.includes("Rin")?S.Rin:0;

    S._emax=-1; S._es=0; S._en=0;
    // vicini con griglia spaziale (lista concatenata su interi: niente stringhe, niente array temporanei)
    const H=2.5*r0, G=Math.ceil(2*W/H)+1;
    if(!grid||grid.G!==G) grid={G,head:new Int32Array(G*G*G)};
    const head=grid.head; head.fill(-1); if(nxt.length<n){ const m=Math.max(2*n,1024); nxt=new Int32Array(m); PX=new Float64Array(m); PY=new Float64Array(m); PZ=new Float64Array(m); PR=new Float64Array(m); PD=new Uint8Array(m); NX=new Float64Array(m); NY=new Float64Array(m); NZ=new Float64Array(m); KC=new Float64Array(m); AX=new Float64Array(m); AY=new Float64Array(m); AZ=new Float64Array(m); AP=new Uint8Array(m); }
    const gi=v=>{const q=Math.floor((v+W)/H); return q<0?0:q>=G?G-1:q;};
    S.cos={n:0,s:0};
    for(let i=0;i<n;i++){const c=cells[i]; c.nb.length=0; c.prs=0; c.anc=0; if(!c.f) c.f=[0,0,0]; c.f[0]=0;c.f[1]=0;c.f[2]=0;
      PX[i]=c.x[0];PY[i]=c.x[1];PZ[i]=c.x[2];PR[i]=c.r;PD[i]=dying(c)?1:0;
      // costringi: l'asse è "verso l'esterno", cioè via dal centro stimato, se la primitiva non ne indica un altro
      AP[i]=c.apic?1:0; if(AP[i]){AX[i]=c.apic[0];AY[i]=c.apic[1];AZ[i]=c.apic[2];}
      KC[i]=c.con||0; if(KC[i]){ let ax=c.conN; if(!ax){ ax=[c.x[0]-c.c[0],c.x[1]-c.c[1],c.x[2]-c.c[2]]; }
        const al=Math.sqrt(ax[0]*ax[0]+ax[1]*ax[1]+ax[2]*ax[2]); if(al>1e-9){NX[i]=ax[0]/al;NY[i]=ax[1]/al;NZ[i]=ax[2]/al;} else KC[i]=0; }
      c._gx=gi(c.x[0]);c._gy=gi(c.x[1]);c._gz=gi(c.x[2]); const k=(c._gz*G+c._gy)*G+c._gx; nxt[i]=head[k]; head[k]=i;}
    for(let i=0;i<n;i++){const a=cells[i],ax=PX[i],ay=PY[i],az=PZ[i],ar=PR[i],da=PD[i];
      const z0=Math.max(0,a._gz-1),z1=Math.min(G-1,a._gz+1),y0=Math.max(0,a._gy-1),y1=Math.min(G-1,a._gy+1),x0=Math.max(0,a._gx-1),x1=Math.min(G-1,a._gx+1);
      for(let z=z0;z<=z1;z++)for(let y=y0;y<=y1;y++)for(let x=x0;x<=x1;x++)for(let j=head[(z*G+y)*G+x];j>=0;j=nxt[j]){ if(j<=i) continue;
        const dx=PX[j]-ax, dy=PY[j]-ay, dz=PZ[j]-az, rs0=ar+PR[j], lim=rs0*1.25, d2=dx*dx+dy*dy+dz*dz; if(d2>=lim*lim) continue;
        const b=cells[j], dd=Math.sqrt(d2)||0.01; a.nb.push(b); b.nb.push(a);
        // costringi: due cellule che costringono la stessa faccia accorciano il contatto fra loro. La fila che
        // costringe diventa più corta di quella sotto, e il foglio si incurva verso la faccia accorciata.
        // La modifica è della coppia, non della singola cellula: le forze restano appaiate per costruzione.
        // Un contatto consolidato ha adottato la lunghezza accorciata come lunghezza di riposo:
        // rimodellamento giunzionale. È per contatto, non per cellula, e non stocca energia.
        const cop=KC[i]&&KC[j];
        let fatt=cop?0.5*(KC[i]+KC[j]):0;
        if(!cop&&a.gnz&&b.gnz){ const f=a.gnz.get(b.id); if(f!=null&&b.gnz.get(a.id)!=null) fatt=f; }
        const rs=fatt?rs0*(1-fatt):rs0;
        if(fatt){ S.cos.n++; S.cos.s+=(dd-rs0)/rs0; }
        // adesione apicale: se il contatto cade sulla faccia apicale di entrambe, non si incollano. È quello che
        // tiene aperto un lume: le pareti che lo guardano non si attaccano, e chi ci galleggia dentro non ha appigli.
        let apap=1;
        if(AP[i]&&AP[j]){ const ux=dx/dd,uy=dy/dd,uz=dz/dd;
          if(ux*AX[i]+uy*AY[i]+uz*AZ[i]>0.5 && -(ux*AX[j]+uy*AY[j]+uz*AZ[j])>0.5) apap=K.adesioneApicale; }
        let f=0; if(dd<rs) f=-(rs-dd)*0.3; else if(!da&&!PD[j]) f=(dd-rs)*0.12*apap*(S.adesione?S.adesione(a,b):1); // adesione differenziale: le primitive possono cambiarla
        if(!da&&!PD[j]){ const e=(dd-rs0)/rs0; if(e>S._emax) S._emax=e; S._es+=e; S._en++;
          if(e<0){ a.prs-=e; b.prs-=e; } // sovrapposizione: quanto le due si compenetrano, in frazione del raggio
          const ad=apap*(S.adesione?S.adesione(a,b):1); a.anc+=ad; b.anc+=ad; } // ancoraggio: quanta adesione totale ha intorno
        if(f){const k=f/dd; a.f[0]+=dx*k;a.f[1]+=dy*k;a.f[2]+=dz*k; b.f[0]-=dx*k;b.f[1]-=dy*k;b.f[2]-=dz*k;} } }
    // fisica + rotazione del corpo (il vettore ruota con la cellula)
    for(const c of cells){const nz=(c.state==="MIGRA"?0.18:0.04)*sc*0.5, dim=piano?2:3;
      if(piano){c.f[2]=0;c.x[2]=0;c.c[2]=0;}
      for(let a=0;a<dim;a++){const old=c.x[a]; let v=old+c.f[a]+rnd()*nz; v=v<-W+r0?-W+r0:v>W-r0?W-r0:v; c.x[a]=v; c.c[a]+=v-old;}
      if(S.rot>0&&!dying(c)){const Rm=rotM(piano?[0,0,1]:randUnit(),S.rot*rnd()*(c.state==="MIGRA"?3:1));c.Q=mm(Rm,c.Q);c.c=add(c.x,mv(Rm,sub(c.c,c.x)));}}
    // canale 3: lettura della matrice
    let readers=0,live=0;
    // la matrice cambia lentamente: ogni cellula la rilegge ogni 3 tick (sfasata), e intanto usa l'ultima lettura
    for(const c of cells){if(dying(c)){c.md=null;continue;}live++;
      if(!S.matOn||c.rho<LLD){c.md=null;c._mdT=null;continue;}
      if(c._ph==null) c._ph=(Math.random()*3)|0;
      if(c._mdT!=null&&(S.tick+c._ph)%3!==0){ if(c.md) readers++; continue; }
      c._mdT=S.tick; c.md=null; const v=sub(c.c,c.x);
      const q=mxRead(c.x,v); if(q.dens<MX.soglia||q.coh<MX.coerenza)continue;
      let nrm=q.n; c._nrm=nrm; let sg=0,rks=0,rka=0,rkn=0;
      // curvatura: due tangenti, la normale letta un passo più in là
      if(MX.curvLocale){
        // curvatura locale: confronto la mia normale con quella che la vicina legge sotto di sé. L'angolo diviso la
        // distanza è la curvatura, il suo inverso è ρ. Niente letture fuori dal proprio corpo.
        for(const o of c.nb){ if(dying(o)||!o._nrm) continue; let b=o._nrm; const sgn=dot(b,nrm)<0?-1:1; b=mul(b,sgn);
          const ang=Math.acos(Math.max(-1,Math.min(1,dot(b,nrm)))); if(ang<MX.curvMin) continue;
          const dv=sub(o.x,c.x), L=len(dv)||1, t=mul(dv,1/L); const dlt=sub(b,nrm);
          // somma lunghezze e angoli e divide alla fine: la media di L/ang ha una coda pesante, questa no
          sg+=dot(dlt,t)>0?-1:1; rks+=L; rka+=ang; rkn++; }
      } else {
        // vecchia lettura: campiona la matrice a un passo di distanza, fuori dal proprio corpo
        const t1=piano?norm([-nrm[1],nrm[0],0]):tangent(nrm), t2=piano?null:cross(nrm,t1);
        for(const t of (t2?[t1,t2]:[t1])){const qb=mxRead(add(c.x,mul(t,LP)),nrm);if(qb.dens<MX.soglia)continue;let b=qb.n;if(dot(b,nrm)<0)b=mul(b,-1);
          const dlt=sub(b,nrm),ang=Math.acos(Math.min(1,dot(b,nrm)));if(ang<MX.curvMin)continue; sg+=dot(dlt,t)>0?-1:1; rks+=LP; rka+=ang; rkn++;}
      }
      // dlt·t > 0: la normale si apre allontanandosi, quindi nrm punta verso l'esterno e il centro è dalla parte opposta
      let s=sg>0?1:sg<0?-1:(dot(nrm,v)>=0?1:-1);
      c.md={n:mul(nrm,s),w:Math.min(1,q.dens)*q.coh,rk:rka>1e-6?rks/rka:0,mr:q.rho,sg:sg!==0}; readers++;}
    S.readers=live?readers/live:0;
    // canale 4: media della curvatura sulla zona
    for(const c of cells){if(dying(c))continue;let s=c.g,k=1;for(const o of c.nb)if(!dying(o)){s+=o.g;k++;}let gn=s/k;
      if(c.md&&c.md.rk>0&&c.rho>LCR){const m=Math.max(-1,Math.min(1,Math.log(c.md.rk/c.rho)-MX.curvBias));gn=(1-MX.curvMedia)*gn+MX.curvMedia*m;}c.gn=gn;}
    for(const c of cells)if(!dying(c))c.g=c.gn;
    // canale 2 (+1, +3, +4): consenso su ρ e sul centro
    const verso=(o,ref)=>{const md=o.md;if(!md)return null;let d=md.n;if(!md.sg&&dot(d,ref)<0)d=mul(d,-1);return d;};
    for(const c of cells){c.nc=c.c;if(dying(c))continue;const v=sub(c.c,c.x),vm=verso(c,v);
      if(vm){let sr=c.rho,k=1; if(c.md.mr>0){const wm=MX.pesoMemoria*c.md.w;sr+=wm*c.md.mr;k+=wm;}
        // ρ_mio = ρ_vicino − (spostamento da lui a me) proiettato sulla direzione del centro;
        // direzione: la media dei centri indicati da me e dalle vicine (ognuna lungo la propria fibra).
        // Una normale letta da sola sbaglia di qualche grado; la media su una decina di vicine molto meno. (conti scalari)
        const vx=vm[0],vy=vm[1],vz=vm[2], cx=c.x[0],cy=c.x[1],cz=c.x[2]; let Px=cx+vx*c.rho,Py=cy+vy*c.rho,Pz=cz+vz*c.rho,kP=1;
        for(const o of c.nb){if(dying(o)||!o.md)continue; const m=o.md.n; let ox=m[0],oy=m[1],oz=m[2]; if(!o.md.sg&&ox*vx+oy*vy+oz*vz<0){ox=-ox;oy=-oy;oz=-oz;}
          let nx=vx+ox,ny=vy+oy,nz=vz+oz; const nl=Math.sqrt(nx*nx+ny*ny+nz*nz)||1; nx/=nl;ny/=nl;nz/=nl; const s1=1+S.noise*rnd();
          sr+=o.rho*(1+S.noise*rnd()*0.1)-((cx-o.x[0])*nx+(cy-o.x[1])*ny+(cz-o.x[2])*nz)*s1; k++;
          Px+=o.x[0]+ox*o.rho; Py+=o.x[1]+oy*o.rho; Pz+=o.x[2]+oz*o.rho; kP++;}
        let rn=Math.max(0,c.rho+CO.alfa*(sr/k-c.rho)); if(c.g)rn*=1+MX.curvGuadagno*c.g;
        let dx=Px/kP-cx,dy=Py/kP-cy,dz=Pz/kP-cz; const dl=Math.sqrt(dx*dx+dy*dy+dz*dz); let dir=dl>1e-6?[dx/dl,dy/dl,dz/dl]:vm;
        c.nc=add(c.x,mul(dir,rn)); continue;}
      let p=[0,0,0],kp=0;
      for(const o of c.nb){if(dying(o)||!o.md)continue;const d=verso(o,v);p=add(p,add(o.x,mul(d,o.rho)));kp++;}
      if(kp){const tgt=mul(p,1/kp);c.nc=add(c.c,mul(sub(tgt,c.c),CO.alfa));continue;}
      let s=c.c,k=1;
      for(const o of c.nb){if(dying(o))continue;
        const Rrel=mm(c.Q,[o.Q[0],o.Q[3],o.Q[6],o.Q[1],o.Q[4],o.Q[7],o.Q[2],o.Q[5],o.Q[8]]); // Q_c · Q_oᵀ
        const Rn=rotM(randUnit(),S.noise*rnd()); const vo=mv(Rn,mv(Rrel,sub(o.c,o.x)));
        s=add(s,add(add(c.x,mul(sub(o.x,c.x),1+S.noise*rnd())),vo));k++;}
      c.nc=add(c.c,mul(sub(mul(s,1/k),c.c),CO.alfa));}
    for(const c of cells){c.c=c.nc;c.rho=len(sub(c.c,c.x));}
    // canale 3: deposito (normale = direzione dal centro)
    // La lamina si depone sul lato basale, non attorno. È questo che tiene il lume permanentemente senza matrice:
    // nessuna cellula viva ne versa dentro, e chi ci finisce non trova appigli.
    if(S.matOn)for(const c of cells){if(dying(c)||c.rho<LVC)continue;
      if(c.fr>K.stacco&&c.age>K.giovane)continue;
      const q=c.apic?[c.x[0]-c.apic[0]*r0*K.basale,c.x[1]-c.apic[1]*r0*K.basale,c.x[2]-c.apic[2]*r0*K.basale]:c.x;
      mxDeposit(q,norm(sub(c.x,c.c)),c.rho);}
    // decisioni
    const born=[];let act=0;
    for(const c of cells){if(dying(c)){c.r-=0.18*sc;continue;}c.age++;
      // impegna: la cellula non decide più (niente divisioni, niente apoptosi). Serve durante le metamorfosi,
      // altrimenti la manutenzione della forma disfa le deformazioni (legge di conflitto).
      if(c.impegnata){if(c.state==="MIGRA"||c.state==="PROLIFERA")c.state="QUIESCENTE";c.migra=0;continue;}
      // Fuori fascia, i due bordi, con la stessa regola applicata due volte. Sopra R: subito. Sotto Rin: solo se
      // non è più giovane. Una cellula nata sotto Rin è il fronte che sale verso la fascia e deve poterci arrivare —
      // è la stessa condizione che le permette di mettere una figlia più in fuori. Passata la giovinezza il fronte
      // è arrivato, e chi è rimasto indietro è fuori fascia come chi ha superato R.
      c.out=(c.rho>R+r0/2||(rin>0&&c.rho<rin-r0/2&&c.age>=K.giovinezza))?c.out+1:0;
      const fuoriCurv=S.matOn&&c.md&&c.g&&c.rho*Math.exp(c.g)>R+LFM;c.fk=fuoriCurv?c.fk+1:0;
      // direzioni libere: la figlia non deve sovrapporsi a nessun vicino
      // la ricerca delle direzioni libere costa (48 direzioni per vicina): si rifà ogni 3 tick, sfasata tra le cellule
      let free=0,best=null,bestScore=-Infinity,gx=0,gy=0,gz=0; const nbl=c.nb.filter(o=>!dying(o)), lim2=(2*r0*0.9)**2;
      if(c._ph==null) c._ph=(Math.random()*3)|0; const rifai=c._fc==null||(S.tick+c._ph)%3===0;
      if(!rifai){ free=c._fc.free; best=c._fc.best; }
      else if(nbl.length<K.pienoDa) for(const d of DD){const px=c.x[0]+d[0]*STEP,py=c.x[1]+d[1]*STEP,pz=c.x[2]+d[2]*STEP;let mind=Infinity;
        for(const o of nbl){const ex=px-o.x[0],ey=py-o.x[1],ez=pz-o.x[2],q=ex*ex+ey*ey+ez*ez;if(q<mind){mind=q;if(mind<lim2)break;}}
        if(mind<lim2)continue;mind=Math.sqrt(mind);free++; gx+=d[0];gy+=d[1];gz+=d[2]; const p=[px,py,pz];const rd=len(sub(c.c,p));
        if(!(rd<=R-r0*(S.margine||0)&&(rd>=rin||(rd>c.rho+1&&c.age<K.giovinezza))))continue;
        // Nessun vincolo di matrice preesistente sulla crescita: un epitelio che si espande depone lamina nuova
        // man mano, e pretendere la matrice prima della figlia impedisce al tessuto di allargarsi (misurato: il
        // bordo rientra di una dozzina di unità). L'ancoraggio vincola la sopravvivenza, non l'espansione.

        // divisione planare: nel piano del foglio. Una figlia spinta oltre la faccia apicale finirebbe nel lume,
        // dove non ha né matrice né giunzioni: è lì che un epitelio tiene aperta una cavità.
        if(c.apic&&dot(d,c.apic)>K.planare) continue;
        const score=c.rho<rin?rd:mind; if(score>bestScore){bestScore=score;best={d,rd};}}
      if(rifai) c._fc={free,best,g:[gx,gy,gz]};
      // direzione della faccia libera (baricentro dello spazio vuoto). Per una cellula con un solo lato
      // libero — la fila esterna di una parete — dice da che parte è la superficie. Serve a marcare.
      if(rifai){ const G=c._fc.g, l=Math.sqrt(G[0]*G[0]+G[1]*G[1]+G[2]*G[2]);
        c.lib=(free>0&&l>1e-6)?[G[0]/l,G[1]/l,G[2]/l]:null; }
      // Polarità apico-basale dal gradiente della matrice. Il baricentro dello spazio libero non basta: in un
      // monostrato la cellula ha facce libere sui due lati e la media si annulla. Qui legge la lamina attorno a sé
      // e punta l'apicale dove la lamina manca — è l'informazione che una cellula vera usa (contatto con la lamina
      // basale). Confronto relativo, senza soglie assolute, e si rinforza da sé perché la deposizione è basale.
      if(rifai){ c.apic=null;
        if(S.matOn&&nbl.length>=(piano?3:5)){ const D=r0*1.4; let g0=0,g1=0,g2=0,sm=0;
          for(const u of POL){ const q=mxRead([c.x[0]+u[0]*D,c.x[1]+u[1]*D,c.x[2]+u[2]*D],null).dens;
            g0+=u[0]*q; g1+=u[1]*q; g2+=u[2]*q; sm+=q; }
          const gl=Math.sqrt(g0*g0+g1*g1+g2*g2);
          if(sm>1e-9&&gl>K.polarita*sm*0.25) c.apic=[-g0/gl,-g1/gl,-g2/gl]; } }
      // figliata non vitale: una cellula la cui figlia è morta appena nata smette di dividersi per un po'.
      // Senza questa regola il genitore rifà la stessa figlia nello stesso posto all'infinito (ciclo futile).
      if(c.pausa>0){ c.pausa--; best=null; }
      c.blk=best?0:c.blk+1;
      // rallenta (c.ferma): arresto mitotico locale nella zona che si deforma. La cellula continua a percepire dove si
      // trova e in che stato è; smette solo di dividersi, e non muore per essere uscita dalla banda, perché lì
      // l'informazione posizionale sta venendo riscritta. Diverso da impegna, che spegne tutto.
      if(rifai) c.fr=free/DD.length;

      // compressa: schiacciata dalle vicine oltre soglia, per abbastanza tick → estrusione (affollamento)
      c.cmp=(c.prs>K.compressione&&!c.ferma)?c.cmp+1:0;
      // staccata: quasi nessuna adesione intorno → anoikis. Non vale per le appena nate, che hanno pochi contatti.
      c.stc=((c.anc<K.ancoraggioMin||c.fr>K.stacco)&&c.age>K.giovane&&!c.ferma)?c.stc+1:0;
      const x={rho:c.rho,R,Rin:rin,r:r0,esposta:free>0,puoCrescere:!!best&&!c.ferma,bloccata:c.blk>K.pazienza&&!c.ferma,
               fuori:(c.out>K.fuoriTick||c.fk>MX.fuoriTick)&&!c.ferma,
               compressa:c.cmp>K.compressioneTick, staccata:c.stc>K.anoikisTick};
      const prima=c.state; c.state=S.m.decide(c,x);
      // un evento reale: una cellula comincia a morire. Lo segna anche alle vicine, che "vedono" succedere qualcosa
      // una cellula che muore avvisa le vicine (le morenti rilasciano segnali). Chi si è diviso da poco e riceve
      // l'avviso smette di dividersi per un po': non serve sapere di chi era figlia.
      if(dying(c)&&!S.m.muoiono.includes(prima)){ c.ev=S.tick; const neonata=c.age<K.giovane;
        for(const o of c.nb){ o.ev=S.tick;
          if(neonata&&!dying(o)&&S.tick-(o.divisaT||-9999)<K.giovane) o.pausa=K.pausaFigliata; } }
      if(dying(c)&&!S.m.muoiono.includes(prima)&&S.matOn){ mxClear(c.x,r0); S.killed++; } // il vuoto lasciato è vuoto davvero
      if(c.migra>0)c.migra--;
      const divides=S.m.dividono.includes(c.state)&&best; if(divides||c.state==="MIGRA")act++;
      if(divides&&--c.timer<=0){c.timer=55+Math.random()*60;const p=add(c.x,mul(best.d,STEP));const nc=newCell(p,c.c.slice(),c.Q.slice(),best.rd);nc.g=c.g;c.ev=S.tick;nc.ev=S.tick;c.divisaT=S.tick;born.push(nc);}}
    S.cells=cells.filter(c=>c.r>sc).concat(born); S.calm=act===0?S.calm+1:0; S.tick++; S.act=act;
    S.tensione={media:S._en?S._es/S._en:0, max:S._emax, n:S._en}; return act;}

  // osservatore: centro condiviso = mediana dei centri stimati
  function centro(){const L=S.cells.filter(c=>!dying(c));if(!L.length)return [0,0,0];const med=k=>{const a=L.map(c=>c.c[k]).sort((p,q)=>p-q);return a[a.length>>1];};return [med(0),med(1),med(2)];}
  function stats(){const L=S.cells.filter(c=>!dying(c)),C=centro();if(!L.length)return{n:0};
    const rs=L.map(c=>len(sub(c.x,C))); const bins=new Array(32).fill(0),binsIn=new Array(32).fill(1e9);
    const NB=piano?16:26;
    L.forEach((c,i)=>{const d=norm(sub(c.x,C));let bi=0,bd=-2;for(let k=0;k<NB;k++){const e=DD[k];const v=dot(d,e);if(v>bd){bd=v;bi=k;}}bins[bi]=Math.max(bins[bi],rs[i]);binsIn[bi]=Math.min(binsIn[bi],rs[i]);});
    const err=L.reduce((s,c,i)=>s+Math.abs(rs[i]-c.rho),0)/L.length; const cnt={};L.forEach(c=>cnt[c.state]=(cnt[c.state]||0)+1);
    const filled=bins.filter(b=>b>0);
    return {n:L.length,rmin:Math.min(...filled),rmax:Math.max(...filled),imin:Math.min(...binsIn),imax:Math.max(...binsIn.filter(b=>b<1e9)),err,cnt,center:C,readers:S.readers,killed:S.killed};}
  // ferite
  function ablate(f,clear){const C=centro();let h=0;for(const c of S.cells){if(dying(c))continue;if(f(sub(c.x,C))){c.state="APOPTOSI";h++;if(clear){mxClear(c.x,r0+2);if(S.onClear)S.onClear(c.x,r0+2);}}}S.calm=0;return h;}
  function dividi(shift,clear){const C=centro();for(const c of S.cells){if(dying(c))continue;const dx=c.x[0]-C[0];if(Math.abs(dx)<MX.striscia){c.state="APOPTOSI";continue;}
      const s=dx>0?shift:-shift;c.x=[c.x[0]+s,c.x[1],c.x[2]];c.c=[c.c[0]+s,c.c[1],c.c[2]];}
    const sh=Math.round(shift/VS); const D0=D.slice(),T0=T.slice(),R0=MR.slice(); D.fill(0);T.fill(0);MR.fill(0);
    if(!clear)for(let k=0;k<VN;k++)for(let j=0;j<VN;j++)for(let i=0;i<VN;i++){const px=-W+(i+.5)*VS-C[0];if(Math.abs(px)<MX.striscia)continue;const ni=i+(px>0?sh:-sh);if(ni<0||ni>=VN)continue;
      const a=vox(i,j,k),b=vox(ni,j,k);D[b]=D0[a];MR[b]=R0[a];for(let t=0;t<6;t++)T[6*b+t]=T0[6*a+t];}
    S.calm=0;}
  reset();
  // consolida(): ogni coppia che sta costringendo fissa il proprio fattore come nuova lunghezza di riposo
  function consolida(){ let n=0; for(const c of S.cells){ if(dying(c)||!c.con) continue;
    for(const o of c.nb){ if(dying(o)||!o.con) continue;
      if(!c.gnz) c.gnz=new Map(); if(!o.gnz) o.gnz=new Map();
      const f=0.5*(c.con+o.con); c.gnz.set(o.id,f); o.gnz.set(c.id,f); n++; } }
    return n; }
  // genera(c,d): una figlia nella direzione d. Serve a un programma che governa una sottopopolazione
  // con una regola propria, invece di lasciarla alla macchina della forma madre.
  function genera(c,d){ const p=[c.x[0]+d[0]*STEP,c.x[1]+d[1]*STEP,c.x[2]+d[2]*STEP];
    const nc=newCell(p,c.c.slice(),c.Q.slice(),len(sub(c.c,p))); nc.g=c.g; S.cells.push(nc); S.calm=0; return nc; }
  return {S,step,reset,stats,centro,ablate,dividi,dying,W,r:r0,piano,consolida,genera,DD,STEP,
    set(o){if(o.m)S.m=M[o.m];for(const k of ["R","Rin","rot","noise","matOn","margine"])if(k in o)S[k]=o[k];},
    density:p=>{const q=mxRead(p,null);return q.dens||0;}};
}
window.Motore3D=Motore;
})();
