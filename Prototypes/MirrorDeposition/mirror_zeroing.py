"""Claim 5's rationale: on a planar EB, do the IMAGES ALONE read phi = n in cut cells?"""
import numpy as np
from mirror_test import CASES, DX, KAPPA_BINS, N_VALID, PAD, box_fluid_volume, kappa_field, sample_random, scoreable_mask
import mirror_planar_ext as mpe
SHAPE=(N_VALID+2*PAD,)*3; LO=-PAD; SIDE=N_VALID*DX
rng=np.random.default_rng(20260817)
print("### Planar EB. 'both' = real + images (the correct build).")
print("### 'images only' = the build that zeroed the real particles' deposit.\n")
print(f"    {'case':<18} {'kappa bin':>13} {'both':>9} {'images only':>13}")
worst=0.0
for name,raw,p0 in CASES:
    n=raw/np.linalg.norm(raw)
    kap=kappa_field(n,p0,SHAPE,LO)
    pos=sample_random(n,p0,600_000,SIDE,rng)
    mass=box_fluid_volume(n,p0,SIDE)/pos.shape[0]
    d=(pos-p0)@n
    img=pos-2.0*d[:,None]*n[None,:]
    sel=(d>0)&(d<=1.5*DX*np.abs(n).sum())
    base=(pos,np.full(pos.shape[0],mass))
    imgs=(img[sel],np.full(int(sel.sum()),mass))
    fb=mpe.deposit_sets(kap,[base,imgs]); fi=mpe.deposit_sets(kap,[imgs])
    mask=scoreable_mask(kap,LO,SHAPE,SIDE)
    for lo,hi in KAPPA_BINS:
        s=mask&(kap>lo)&(kap<=hi); nn=int(s.sum())
        if nn<3: continue
        b=float(np.mean(fb[s]-1.0)); i=float(np.mean(fi[s]-1.0))
        worst=max(worst,abs(i))
        print(f"    {name:<18} {lo:5.3f}-{hi:5.3f} {b:+9.4f} {i:+13.4f}")
print(f"\n    Largest |deviation| of the images-only build over all planar cut-cell bins: {100*worst:.1f}%")
