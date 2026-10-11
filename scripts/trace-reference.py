from pathlib import Path
from PIL import Image,ImageDraw
import numpy as np
from scipy.ndimage import gaussian_filter,sobel,binary_dilation,label
import json,base64,io,sys
root=Path(__file__).resolve().parent.parent

source=sys.argv[1] if len(sys.argv)>1 else io.BytesIO(base64.b64decode((root/'assets/detail-reference.svg').read_text().split('base64,')[1].split(chr(34))[0]))
im=Image.open(source).convert('RGB');rgb=np.asarray(im).astype(float)/255
gray=rgb@np.array([.299,.587,.114]);blur=gaussian_filter(gray,.75)
gx=sobel(blur,axis=1)/8;gy=sobel(blur,axis=0)/8;mag=np.hypot(gx,gy);angle=(np.arctan2(gy,gx)*180/np.pi)%180
nms=np.zeros_like(mag)
for mask,dy,dx in [((angle<22.5)|(angle>=157.5),0,1),((angle>=22.5)&(angle<67.5),1,1),((angle>=67.5)&(angle<112.5),1,0),((angle>=112.5)&(angle<157.5),1,-1)]:
 good=mask&(mag>=np.roll(mag,(dy,dx),(0,1)))&(mag>=np.roll(mag,(-dy,-dx),(0,1)));nms[good]=mag[good]
strong=nms>.065;weak=nms>.025;labs,n=label(weak,np.ones((3,3)));keep=np.unique(labs[strong]);keep=keep[keep>0];edges=np.isin(labs,keep);edges[[0,-1],:]=False;edges[:,[0,-1]]=False
# Thin nonmaximum plateaus, then split the skeleton into independent pen strokes.
for it in range(40):
 changed=False
 for phase in (0,1):
  a=np.pad(edges,1);p=[a[:-2,1:-1],a[:-2,2:],a[1:-1,2:],a[2:,2:],a[2:,1:-1],a[2:,:-2],a[1:-1,:-2],a[:-2,:-2]]
  count=sum(q.astype(int) for q in p);trans=sum((~p[i])&p[(i+1)%8] for i in range(8))
  cond=(~(p[0]&p[2]&p[4])&~(p[2]&p[4]&p[6])) if phase==0 else (~(p[0]&p[2]&p[6])&~(p[0]&p[4]&p[6]))
  remove=edges&(count>=2)&(count<=6)&(trans==1)&cond
  if remove.any():edges[remove]=False;changed=True
 if not changed:break
coords={tuple(x) for x in np.argwhere(edges)};adj={}
for y,x in coords:
 neighbors=[]
 for dy in (-1,0,1):
  for dx in (-1,0,1):
   if (dy or dx) and (y+dy,x+dx) in coords:
    if dy and dx and ((y+dy,x) in coords or (y,x+dx) in coords):continue
    neighbors.append((y+dy,x+dx))
 adj[(y,x)]=neighbors
used=set();paths=[]
def key(a,b):return tuple(sorted((a,b)))
def walk(a,b):
 out=[a,b];used.add(key(a,b));prev,cur=a,b
 while len(adj[cur])==2:
  nex=next(q for q in adj[cur] if q!=prev)
  if key(cur,nex) in used:break
  out.append(nex);used.add(key(cur,nex));prev,cur=cur,nex
 return out
for a in sorted(coords,key=lambda p:(len(adj[p])==2,p)):
 for b in adj[a]:
  if key(a,b) not in used:
   stroke=walk(a,b)
   if len(stroke)>=4:paths.append(stroke)
# Keep subpixel-sized details, join adjacent fragments by nearest end to reduce pen-up travel.
strokes=[np.array([(x,y) for y,x in p],float) for p in paths]
ordered=[];cur=np.array([0,0])
while strokes:
 ends=np.array([(s[0],s[-1]) for s in strokes]);dist=((ends-cur)**2).sum(2);i,j=np.unravel_index(dist.argmin(),dist.shape);s=strokes.pop(i)
 if j:s=s[::-1]
 ordered.append(s);cur=s[-1]
data=[[{'x':float(x),'y':float(y)} for x,y in s] for s in ordered]
(root/'detail-preset.js').write_text('/* Pixel-traced contours of the user supplied 512 × 512 reference. */\nFourierPresets.detail='+json.dumps(data,separators=(',',':'))+';\n')
commands=['M'+' L'.join(f'{x:g},{y:g}' for x,y in s) for s in ordered]
svg='<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512"><rect width="512" height="512" fill="white"/><path d="'+' '.join(commands)+'" fill="none" stroke="#222" stroke-width=".85" stroke-linecap="round" stroke-linejoin="round"/></svg>'
(root/'assets/detail-contours.svg').write_text(svg)
preview=Image.new('RGB',(512,512),'white');d=ImageDraw.Draw(preview)
for s in ordered:d.line([tuple(p) for p in s],fill='#222222',width=1)
# The SVG is the permanent vector preview; raster preview stays optional.
# Embed original reference in an SVG so it remains self-contained and reviewable.
import io
buf=io.BytesIO();im.save(buf,format='PNG')
(root/'assets/detail-reference.svg').write_text('<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 512 512"><image width="512" height="512" xlink:href="data:image/png;base64,'+base64.b64encode(buf.getvalue()).decode()+'"/></svg>')
print('strokes',len(data),'points',sum(map(len,data)),'edge pixels',len(coords))
