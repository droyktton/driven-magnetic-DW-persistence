import numpy as np, json, glob
from PIL import Image
from scipy import ndimage as ndi
from scipy.special import erf
from scipy.optimize import curve_fit
from skimage.filters import threshold_otsu
SRC='/home/koltona/MOKEdata_local/TbCo_Paris/Datos_DW_raw/AV-6517_1_0_#2/TbCo_1_0_#2_%05d.tif'
R='/home/koltona/Codigos/driven-magnetic-DW-persistence/TbCo_Paris_Nov2023/AV6517_1_0_#2_0.83nm/'
J=json.load(open(R+'pulses.json')); sc={int(k):v for k,v in J['scale'].items()}
ld=lambda k: np.asarray(Image.open(SRC%k)).astype(np.float32)/sc[k]
G=lambda x: ndi.gaussian_filter(x,2)
st=np.stack([G(ld(k)) for k in range(2,56)])
B=np.percentile(st,80,0); D=np.percentile(st,20,0); del st
np.save('B.npy',B.astype(np.float32)); np.save('D.npy',D.astype(np.float32))
mean_img=np.mean([ld(k) for k in range(1,56,5)],0); g8=ndi.gaussian_filter(mean_img,8)
dia=ndi.binary_fill_holes(g8>threshold_otsu(g8)); lab,n=ndi.label(dia); dia=lab==1+np.argmax(ndi.sum(dia,lab,range(1,n+1)))
core=ndi.binary_erosion(dia,iterations=30); np.save('core.npy',core)
W=np.load(R+'walls.npz')
res={}
for k in (12,20,30,36,48,49,50,51):
    nraw=(ld(k)-D)/(B-D)
    # noise: pixels far (> 25 px) from any wall of this file, inside the field
    pts=np.concatenate([W[f] for f in W.files if f.startswith(f'f{k}_')])
    far=np.ones_like(core); far[np.clip(pts[:,0].astype(int),0,2047),np.clip(pts[:,1].astype(int),0,2047)]=False
    far=ndi.distance_transform_edt(far)>25
    lev=G(nraw)
    nb=(far&core&(lev>0.5)); nd=(far&core&(lev<0.5))
    hp=nraw-ndi.gaussian_filter(nraw,6)                              # high-pass: white-noise part
    sig_n=float(np.std(hp[nb|nd]))*1/np.sqrt(1-1/(4*np.pi*36)*0)     # ~ std of the pixel noise
    # edge profiles along the normal of the wall (every 10th point of pieces)
    prof=[]; dd=np.arange(-15,15.25,0.25)
    for f in [f for f in W.files if f.startswith(f'f{k}_')]:
        c=W[f]
        for i in range(10,len(c)-10,10):
            t=c[i+4]-c[i-4]; t/=np.hypot(*t); nrm=np.array([-t[1],t[0]])
            yy=c[i,0]+dd*nrm[0]; xx=c[i,1]+dd*nrm[1]
            p=ndi.map_coordinates(nraw,[yy,xx],order=1)
            if p[:20].mean()>p[-20:].mean(): p=p[::-1]
            prof.append(p)
    P=np.mean(prof,0)
    f_=lambda d,a,b,x0,s: a+b*0.5*(1+erf((d-x0)/(np.sqrt(2)*s)))
    po,_=curve_fit(f_,dd,P,p0=[0,1,0,2])
    res[k]=dict(noise_sd=sig_n,edge_sigma_px=abs(po[3]),step=po[1],n_prof=len(prof))
    print(k,res[k],flush=True)
json.dump(res,open('measured.json','w'))
