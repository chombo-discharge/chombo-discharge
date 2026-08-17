"""Band width for the ACTUAL chombo-discharge kernels (overlap-integral CIC and TSC)."""
import numpy as np
from mirror_test import CASES, DX, N_VALID, PAD, kappa_field, sample_random
SHAPE=(N_VALID+2*PAD,)*3; LO=-PAD; SIDE=N_VALID*DX

def support_offsets(L):
    """Integer cell offsets iv-base whose overlap weight can be nonzero, per direction."""
    # a = iv - pos/dx ; nonzero requires a < L/2 and a+1 > -L/2
    return int(np.ceil(L/2.0))+1

def fires(pos, kappa, L):
    """Does the CIC/TSC cloud of width L centred at pos touch a cell with kappa>0?"""
    base=np.floor(pos/DX).astype(np.int64)
    rad=support_offsets(L)
    out=np.zeros(pos.shape[0],bool)
    for ox in range(-rad,rad+1):
        for oy in range(-rad,rad+1):
            for oz in range(-rad,rad+1):
                iv=base+np.array([ox,oy,oz])
                w=np.ones(pos.shape[0])
                for d in range(3):
                    a=iv[:,d]-pos[:,d]/DX
                    b=a+1.0
                    al=np.maximum(a,-0.5*L); be=np.minimum(b,0.5*L)
                    ok=al<be
                    w*=np.where(ok,(be-al)-(be*np.abs(be)-al*np.abs(al))/L,0.0)
                sh=iv-LO
                ins=np.all((sh>=0)&(sh<np.array(SHAPE)),axis=-1)
                s=np.clip(sh,0,np.array(SHAPE)-1)
                kap=np.zeros(pos.shape[0]); kap[ins]=kappa[s[ins,0],s[ins,1],s[ins,2]]
                out |= (w>1e-14)&ins&(kap>0)
    return out

print("### Measured band, with the deposition kernels the code actually uses")
print("### CIC: cloud width L=1 (reach dx).  TSC: L=2 (reach 1.5 dx).")
print(f"\n    {'case':<18} {'sum|n|':>7} {'3 s_max':>8} {'CIC max d':>10} "
      f"{'4 s_max':>8} {'TSC max d':>10} {'plan band':>10}  TSC verdict")
rng=np.random.default_rng(11)
for name,raw,p0 in CASES:
    n=raw/np.linalg.norm(raw)
    kap=kappa_field(n,p0,SHAPE,LO)
    pos=sample_random(n,p0,900_000,SIDE,rng)
    d=(pos-p0)@n
    img=pos-2.0*d[:,None]*n[None,:]
    sn=float(np.abs(n).sum()); s3=1.5*DX*sn; s4=2.0*DX*sn
    res={}
    for lab,L in (("CIC",1.0),("TSC",2.0)):
        f=(d>0)&fires(img,kap,L)
        res[lab]=float(d[f].max()) if f.any() else float("nan")
    verd = "plan band SHORT by %.3f dx"%(res['TSC']-s3) if res['TSC']>s3+1e-9 else "ok"
    print(f"    {name:<18} {sn:7.3f} {s3:8.3f} {res['CIC']:10.3f} "
          f"{s4:8.3f} {res['TSC']:10.3f} {s3:10.3f}  {verd}")
