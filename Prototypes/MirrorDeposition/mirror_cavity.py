"""Tight cavities: the Jacobian blows up long before it becomes singular.

Revision 6. The revision-3..5 version of this file selected reflecting particles with
the SUPERSEDED distance band `d <= 1.5*dx*sum|n_i|` (BANDMAX = 1.5*sqrt(3) = 2.598) and
reported the analytic amplification at that d. PLAN.md 1.2 no longer states a band as a
distance: the criterion is "reflect iff the image's cloud overlaps at least one kappa > 0
cell", with the Chebyshev-2 mask as an early-out only. This file now applies that
criterion, and reports the measured max d it admits rather than an assumed one.

The deposition runs are CIC (mirror_sphere_ext.deposit is a CIC harness). The analytic
header table is given for BOTH kernels' planar reach, since PLAN.md 3.2's guard has to be
sized for the larger one.
"""
import numpy as np
import mirror_sphere as ms, mirror_sphere_ext as mse
from mirror_sphere import BINS, DX, LO, N_VALID, SHAPE, kappa_sphere

# Planar reach of a reflecting particle, MEASURED by reach_cells.py over ten planar cases.
REACH = {"CIC": 2.570, "TSC": 3.390}
print("### Analytic: on a PLANE the deepest reflecting particle sits at d = %.3f dx (CIC)\n"
      "### and %.3f dx (TSC). A cavity of radius R amplifies an image by\n"
      "### J = ((R+d)/(R-d))^2, singular at d = R.\n" % (REACH["CIC"], REACH["TSC"]))
print(f"    {'cavity R/dx':>12} {'J at d=2.570':>14} {'J at d=3.390':>14}   note")
for R in (10.,8.,6.,5.,4.,3.5,3.0,2.8,2.6,2.4):
    out = []
    for k in ("CIC", "TSC"):
        d = REACH[k]
        out.append(float("inf") if d >= R else ((R+d)/(R-d))**2)
    note = [f"{k} singular" for k, J in zip(("CIC", "TSC"), out) if J == float("inf")]
    print(f"    {R:12.1f} {out[0]:14.1f} {out[1]:14.1f}   {', '.join(note)}")
print("\n### The planar reach is an UPPER bound here: in a cavity the image of a deep\n"
      "### particle lands far inside the solid, where its cloud touches no kappa > 0 cell,\n"
      "### so the criterion prunes it. The measured max d below is what actually happens.\n")

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
    # PLAN.md 1.2's criterion: reflect iff d > 0 AND the image's cloud overlaps a kappa > 0 cell.
    idx, wt = mse.cic_stencil(img)
    sh = idx - LO
    ins = np.all((sh >= 0) & (sh < np.array(SHAPE)), axis=-1)
    s = np.clip(sh, 0, np.array(SHAPE) - 1)
    kap = np.zeros(wt.shape)
    kap[ins] = kappa[s[..., 0][ins], s[..., 1][ins], s[..., 2][ins]]
    reaches = np.any((wt > 0) & ins & (kap > 0), axis=1)
    sel = (d > 0) & reaches
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
          f"  max d {ds.max():6.3f}  Jmax {J.max():9.1f}  worst |dev| {100*w:7.2f}%"
          f"  mass {float(np.sum(fld*kappa)*DX**3)/tot:.5f}")
    for lo,hi,nn,mu,sd in rows:
        print(f"      kappa {lo:5.3f}-{hi:5.3f}  cells {nn:4d}  dev {mu:+9.4f}  cell-to-cell sd {sd:8.4f}")

print("### Deposition (CIC), exact reflection + EXACT Jacobian, tight cavities\n")
for R in (8.0, 6.0, 4.0, 3.5, 3.0):
    run(R)


def reach_scan(R, npart=4_000_000, seed=17):
    """Kernel-resolved: the deepest d the PLAN's criterion admits in a cavity, and the
    resulting Jacobian. Support-only (a cell is reached iff the cloud overlaps it), so it
    needs no per-kernel weight formula -- CIC and TSC share their support machinery."""
    centre = np.array([8.0 + R, 8., 8.])
    kappa, _ = kappa_sphere(centre, R, False, nsub=16)
    rng = np.random.default_rng(seed)
    side = N_VALID * DX
    t = rng.uniform(0, side, size=(npart, 3))
    rt = np.linalg.norm(t - centre, axis=-1)
    pos = t[rt < R]
    if pos.shape[0] < 5000:
        return None
    v = pos - centre
    r = np.linalg.norm(v, axis=-1)
    u = v / r[:, None]
    d = R - r
    img = centre + (2 * R - r)[:, None] * u
    out = {}
    for name, half in (("CIC", 0.5), ("TSC", 1.0)):
        lo = np.floor((img - half * DX) / DX).astype(int)
        hi = np.floor((img + half * DX) / DX).astype(int)
        hit = np.zeros(img.shape[0], bool)
        span = int(round(2 * half)) + 1
        for a in range(span):
            for b in range(span):
                for c in range(span):
                    cell = np.minimum(lo + np.array([a, b, c]), hi)
                    sh = cell - LO
                    ins = np.all((sh >= 0) & (sh < np.array(SHAPE)), axis=-1)
                    s = np.clip(sh, 0, np.array(SHAPE) - 1)
                    k = np.zeros(img.shape[0])
                    k[ins] = kappa[s[:, 0][ins], s[:, 1][ins], s[:, 2][ins]]
                    hit |= ins & (k > 0)
        sel = (d > 0) & hit
        ds = d[sel]
        out[name] = (float(sel.mean()), float(ds.max()), float(((1 + ds / R) / (1 - ds / R)).max() ** 2))
    return out


print("\n### Kernel-resolved reach under the PLAN's criterion (support-only scan)\n")
print(f"    {'cavity R/dx':>12} {'CIC refl%':>10} {'CIC max d':>10} {'CIC Jmax':>10}"
      f" {'TSC refl%':>10} {'TSC max d':>10} {'TSC Jmax':>10}")
for R in (8.0, 6.0, 4.0, 3.5, 3.0):
    o = reach_scan(R)
    if o is None:
        continue
    print(f"    {R:12.1f} {100*o['CIC'][0]:10.1f} {o['CIC'][1]:10.3f} {o['CIC'][2]:10.1f}"
          f" {100*o['TSC'][0]:10.1f} {o['TSC'][1]:10.3f} {o['TSC'][2]:10.1f}")
