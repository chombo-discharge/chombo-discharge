"""Tight cavities: the Jacobian blows up long before it becomes singular."""
import numpy as np
import mirror_sphere as ms, mirror_sphere_ext as mse
from mirror_sphere import BINS, DX, LO, N_VALID, SHAPE, kappa_sphere

BANDMAX = 1.5*np.sqrt(3.0)
print("### Analytic: the band reaches d = %.3f dx, so a cavity of radius R gives\n"
      "### max image weight amplification J = ((R+d)/(R-d))^2 at d = min(band, R)\n" % BANDMAX)
print(f"    {'cavity R/dx':>12} {'max d in band':>14} {'max J':>14}  note")
for R in (10.,8.,6.,5.,4.,3.5,3.0,2.8,2.6,2.4):
    d = min(BANDMAX, R*0.999)
    J = ((R+d)/(R-d))**2
    note = "singular INSIDE band" if R <= BANDMAX else ("J > 20" if J>20 else "")
    print(f"    {R:12.1f} {d:14.3f} {J:14.1f}  {note}")

def run(R, npart=4_000_000, seed=17, margin=3):
    centre=np.array([8.0+R,8.,8.])
    kappa,_=kappa_sphere(centre,R,False,nsub=16); mse.KAPPA=kappa
    rng=np.random.default_rng(seed); side=N_VALID*DX
    t=rng.uniform(0,side,size=(npart,3)); rt=np.linalg.norm(t-centre,axis=-1)
    pos=t[rt<R]
    if pos.shape[0]<5000:
        print(f"  concave R={R}: only {pos.shape[0]} particles, skipped"); return
    mass=side**3/npart
    v=pos-centre; r=np.linalg.norm(v,axis=-1); u=v/r[:,None]
    d=R-r; img=centre+(2*R-r)[:,None]*u; n=-u
    sel=(d>0)&(d<=1.5*DX*np.abs(n).sum(-1))
    ds=d[sel]; J=((1+ds/R)/(1-ds/R))**2
    base=(pos,np.full(pos.shape[0],mass))
    fld=mse.deposit([base,(img[sel],mass*J)])
    ax=[np.arange(SHAPE[q])+LO for q in range(3)]
    ii,jj,kk=np.meshgrid(*ax,indexing="ij"); inter=np.ones(SHAPE,bool)
    for a in (ii,jj,kk): inter&=(a>=margin)&(a<=N_VALID-1-margin)
    m=inter&(kappa>0)
    rows=[]
    for lo,hi in BINS:
        s=m&(kappa>lo)&(kappa<=hi); nn=int(s.sum())
        if nn>=3: rows.append((lo,hi,nn,float(np.mean(fld[s]-1.0)),float(np.std(fld[s]-1.0))))
    w=max(abs(x[3]) for x in rows) if rows else float("nan")
    tot=pos.shape[0]*mass
    print(f"  concave R={R:4.1f}  particles {pos.shape[0]:8d}  reflected {100*sel.mean():5.1f}%"
          f"  Jmax {J.max():9.1f}  worst |dev| {100*w:7.2f}%"
          f"  mass {float(np.sum(fld*kappa)*DX**3)/tot:.5f}")
    for lo,hi,nn,mu,sd in rows:
        print(f"      kappa {lo:5.3f}-{hi:5.3f}  cells {nn:4d}  dev {mu:+9.4f}  cell-to-cell sd {sd:8.4f}")

print("\n### Deposition, exact reflection + EXACT Jacobian, tight cavities\n")
for R in (8.0, 6.0, 4.0, 3.5, 3.0):
    run(R)
