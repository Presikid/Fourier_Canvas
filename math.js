/* Fourier Studio mathematical core. No dependencies. */
(function(root){
'use strict';
const TAU=2*Math.PI;
function fft(re,im,inverse=false){
 const n=re.length;if(n!==im.length||n<1||(n&(n-1)))throw Error('FFT length must be a power of two');
 for(let i=1,j=0;i<n;i++){let b=n>>1;for(;j&b;b>>=1)j^=b;j^=b;if(i<j){[re[i],re[j]]=[re[j],re[i]];[im[i],im[j]]=[im[j],im[i]];}}
 for(let len=2;len<=n;len*=2){const a=(inverse?1:-1)*TAU/len,wr=Math.cos(a),wi=Math.sin(a);for(let i=0;i<n;i+=len){let ur=1,ui=0;for(let j=0;j<len/2;j++){const p=i+j,q=p+len/2,tr=re[q]*ur-im[q]*ui,ti=re[q]*ui+im[q]*ur;re[q]=re[p]-tr;im[q]=im[p]-ti;re[p]+=tr;im[p]+=ti;const t=ur*wr-ui*wi;ui=ur*wi+ui*wr;ur=t;}}}
 if(inverse)for(let i=0;i<n;i++){re[i]/=n;im[i]/=n;}return {re,im};
}
function normalize(strokes){const pts=strokes.flat();if(!pts.length)return [];let minX=Infinity,minY=Infinity,maxX=-Infinity,maxY=-Infinity;for(const p of pts){minX=Math.min(minX,p.x);minY=Math.min(minY,p.y);maxX=Math.max(maxX,p.x);maxY=Math.max(maxY,p.y);}const cx=(minX+maxX)/2,cy=(minY+maxY)/2,s=500/Math.max(maxX-minX,maxY-minY,1);return strokes.filter(s=>s.length>1).map(a=>a.map(p=>({x:(p.x-cx)*s,y:(p.y-cy)*s})));}
function resample(strokes,n,close=true){
 const segs=[];let total=0;
 function add(a,b,pen){const d=Math.hypot(b.x-a.x,b.y-a.y);if(d>1e-9){segs.push({a,b,pen,start:total,len:d});total+=d;}}
 const valid=strokes.filter(s=>s.length>1);if(!valid.length)throw Error('请先画出一条曲线');
 valid.forEach((s,i)=>{for(let j=1;j<s.length;j++)add(s[j-1],s[j],true);if(i<valid.length-1)add(s.at(-1),valid[i+1][0],false);});add(valid.at(-1).at(-1),valid[0][0],close&&valid.length===1);
 if(total<1e-9)throw Error('曲线长度太短');let idx=0;const points=[],pen=[];for(let i=0;i<n;i++){const d=i*total/n;while(idx<segs.length-1&&d>=segs[idx].start+segs[idx].len)idx++;const q=segs[idx],t=(d-q.start)/q.len;points.push({x:q.a.x+t*(q.b.x-q.a.x),y:q.a.y+t*(q.b.y-q.a.y)});pen.push(q.pen);}return {points,pen,length:total};
}
function coefficients(points){const n=points.length,re=points.map(p=>p.x),im=points.map(p=>p.y);fft(re,im);return re.map((r,k)=>({freq:k<=n/2?k:k-n,re:r/n,im:im[k]/n,amp:Math.hypot(r,im[k])/n,phase:Math.atan2(im[k],r)}));}
function order(cs,mode='amplitude'){return [...cs].sort(mode==='frequency'?(a,b)=>Math.abs(a.freq)-Math.abs(b.freq)||a.freq-b.freq:(a,b)=>b.amp-a.amp);}
function evaluate(cs,t){let x=0,y=0;for(const c of cs){const a=TAU*c.freq*t,co=Math.cos(a),si=Math.sin(a);x+=c.re*co-c.im*si;y+=c.re*si+c.im*co;}return {x,y};}
function reconstruct(cs,n){const re=new Array(n).fill(0),im=new Array(n).fill(0);for(const c of cs){const k=(c.freq%n+n)%n;re[k]+=c.re*n;im[k]+=c.im*n;}fft(re,im,true);return re.map((x,i)=>({x,y:im[i]}));}
function metrics(all,selected){const energy=all.reduce((s,c)=>s+c.amp*c.amp,0),kept=selected.reduce((s,c)=>s+c.amp*c.amp,0);return {energy:energy?kept/energy:1,rms:Math.sqrt(Math.max(0,energy-kept))};}
const api={TAU,fft,normalize,resample,coefficients,order,evaluate,reconstruct,metrics};if(typeof module!=='undefined')module.exports=api;root.Fourier=api;
})(typeof globalThis==='undefined'?window:globalThis);
