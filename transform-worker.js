/* Background-only Fourier processing. Transfer numeric buffers, not point objects. */
'use strict';
importScripts('math.js?v=worker-1');
const engine=Fourier.createEngine();
self.onmessage=({data})=>{
 const {id,type,settings}=data;
 try{
  if(type==='export'){const file=engine.exportFile(data.kind,data.style),bytes=new TextEncoder().encode(file.text);self.postMessage({id,type,name:file.name,mime:file.mime,bytes},[bytes.buffer]);return;}
  const result=type==='analyze'?engine.analyze(data.strokes,data.n,data.closed,settings):engine.select(settings);
  const transfers=[];for(const value of Object.values(result))if(ArrayBuffer.isView(value))transfers.push(value.buffer);
  self.postMessage({id,type,result},transfers);
 }catch(error){self.postMessage({id,type,error:error.message});}
};
