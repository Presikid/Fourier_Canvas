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
function displayIndices(n,pen,limit=6000){if(n<=0)return[];const keep=new Set([0,n-1]),step=Math.max(1,Math.ceil(n/limit));for(let i=0;i<n;i+=step)keep.add(i);if(pen)for(let i=0;i<n-1;i++)if(pen[i]!==pen[i+1]){keep.add(i);keep.add(i+1);keep.add(Math.min(n-1,i+2));}return [...keep].sort((a,b)=>a-b);}
function packedView(data,stride,fields,indices=null){const n=indices?indices.length:data.length/stride,get=i=>{if(i<0)i+=n;if(i<0||i>=n)return undefined;const k=(indices?indices[i]:i)*stride,p={};fields.forEach((f,j)=>p[f]=data[k+j]);return p;};const target={length:n,at:get,slice:(start=0,end=n)=>{start=Math.max(0,start<0?n+start:start);end=Math.min(n,end<0?n+end:end);return Array.from({length:Math.max(0,end-start)},(_,i)=>get(start+i));},map:fn=>Array.from({length:n},(_,i)=>fn(get(i),i)),every:fn=>{for(let i=0;i<n;i++)if(!fn(get(i),i))return false;return true;},[Symbol.iterator]:function*(){for(let i=0;i<n;i++)yield get(i);}};return new Proxy(target,{get:(o,k)=>typeof k==='string'&&/^\d+$/.test(k)?get(Number(k)):Reflect.get(o,k)});}
function packPoints(points){const a=new Float64Array(points.length*2);for(let i=0;i<points.length;i++){a[i*2]=points[i].x;a[i*2+1]=points[i].y;}return a;}
function createEngine(){let samples=[],pen=[],all=[],sorted={},last=null,display=[],signalMax=1;
 function analyze(strokes,n,closed,settings){const r=resample(strokes,n,closed);samples=r.points;pen=r.pen;all=coefficients(samples);sorted={};display=displayIndices(n,pen,strokes.reduce((sum,s)=>sum+s.length,0)>20000?24000:6000);signalMax=1;for(const p of samples)signalMax=Math.max(signalMax,Math.abs(p.x),Math.abs(p.y));return select(settings,true);}
 function select(settings,includeSource=false){if(!all.length)throw Error('尚无曲线可供计算');const mode=settings.sort==='frequency'?'frequency':'amplitude';if(!sorted[mode])sorted[mode]=order(all,mode);const ordered=sorted[mode];let count=Math.max(1,Math.min(all.length,Math.round(settings.harmonics)||1));if(settings.autoEnergy){const energy=all.reduce((sum,c)=>sum+c.amp*c.amp,0);let kept=0;count=0;while(count<all.length&&(count===0||kept<energy*.99))kept+=ordered[count++].amp**2;}
 const selected=ordered.slice(0,count),reconstruction=reconstruct(selected,samples.length),selectedIndices=new Uint32Array(count),kept=new Set();for(let i=0;i<count;i++){const k=(selected[i].freq+all.length)%all.length;selectedIndices[i]=k;kept.add(k);}
 const bins=Math.min(100,all.length),bars=[],n=all.length,first=-Math.floor((n-1)/2);let max=1e-6;for(const c of all)max=Math.max(max,c.amp);for(let b=0;b<bins;b++){let peak=null;for(let j=Math.floor(b*n/bins);j<Math.floor((b+1)*n/bins);j++){const freq=first+j,c=all[(freq+n)%n];if(!peak||c.amp>peak.amp)peak=c;}bars.push({amp:peak.amp,retained:kept.has((peak.freq+n)%n)});}
 const rows=ordered.slice(0,100).map(c=>({...c,retained:kept.has((c.freq+n)%n)})),m={energy:1,rms:0};let total=0,omitted=0;for(let k=0;k<n;k++){const power=all[k].amp**2;total+=power;if(!kept.has(k))omitted+=power;}m.energy=total?Math.max(0,1-omitted/total):1;m.rms=Math.sqrt(omitted);
 last={selected,reconstruction,kept,settings:{...settings}};
 const result={selectedIndices,reconstruction:packPoints(reconstruction),metrics:m,count,spectrum:{bars,max,first,last:first+n-1},rows,signalMax};
 if(includeSource){const packed=new Float64Array(n*5);all.forEach((c,i)=>{packed[i*5]=c.freq;packed[i*5+1]=c.re;packed[i*5+2]=c.im;packed[i*5+3]=c.amp;packed[i*5+4]=c.phase;});result.coefficients=packed;result.samples=packPoints(samples);result.pen=Uint8Array.from(pen);result.display=Uint32Array.from(display);}
 return result;}
 function exportFile(kind,style){if(!last)throw Error('请先完成计算');if(kind==='csv'){let text='frequency,real,imaginary,amplitude,phase_radians,retained\n',n=all.length,first=-Math.floor((n-1)/2);const lines=[];for(let j=0;j<n;j++){const k=(first+j+n)%n,c=all[k];lines.push([c.freq,c.re,c.im,c.amp,c.phase,last.kept.has(k)?1:0].join(','));}return{text:text+lines.join('\n'),name:'fourier-coefficients.csv',mime:'text/csv'};}
 if(kind!=='svg')throw Error('未知导出格式');const pts=last.reconstruction,n=pts.length,commands=[];let minX=Infinity,minY=Infinity,maxX=-Infinity,maxY=-Infinity;for(let i=0;i<n;i++){const p=pts[i];minX=Math.min(minX,p.x);maxX=Math.max(maxX,p.x);minY=Math.min(minY,p.y);maxY=Math.max(maxY,p.y);}for(let i=0;i<=n;i++){const j=i%n,p=pts[j],visible=!style.penLift||pen[(j-1+n)%n];commands.push((i===0||!visible?'M':'L')+p.x.toFixed(3)+' '+p.y.toFixed(3));}return{name:'fourier-curve.svg',mime:'image/svg+xml',text:'<svg xmlns="http://www.w3.org/2000/svg" viewBox="'+[minX-15,minY-15,maxX-minX+30,maxY-minY+30].join(' ')+'"><path d="'+commands.join(' ')+'" fill="none" stroke="'+style.color+'" stroke-width="'+style.lineWidth+'" stroke-linecap="round" stroke-linejoin="round"/></svg>'};}
 return{analyze,select,exportFile};}
const api={TAU,fft,normalize,resample,coefficients,order,evaluate,reconstruct,metrics,displayIndices,packedView,packPoints,createEngine};if(typeof module!=='undefined')module.exports=api;root.Fourier=api;
})(typeof globalThis==='undefined'?window:globalThis);
