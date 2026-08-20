"""
Fourth reflection source: a QUADRATIC surface patch stored per band cell.

Stored per band cell, all of it available from the cut cells themselves:
    x_j   boundary centroid of the nearest cut cell   (3 Reals)
    n_j   its discrete normal                         (3 Reals)
    S_j   its fitted 2x2 shape operator               (3 Reals, symmetric)

Per particle at p, in the frame (n_j, t1, t2) with w = p - x_j:
    xi = (t1.w, t2.w)       eta = n_j.w
    d      = eta + 0.5 xi^T S xi          (plane term + curvature correction)
    n_hat  = normalize(n_j + S xi in the tangent plane)
    image  = p - 2 d n_hat
No implicit function, no per-particle geometry query -- one dependent load.

Cut-cell normals/centroids/S are analytic here, to isolate the surface MODEL.
discnorm.py measures what the discrete inputs cost on top (~1-2%).
"""
import numpy as np
import mirror_sphere as ms, mirror_sphere_ext as mse
from mirror_sphere import BINS, DX, LO, N_VALID, SHAPE, cell_centres, kappa_sphere


def own_foot_planes(centre, radius, fluid_outside):
    """Per band cell: nhat and boundary point at THAT CELL'S OWN foot point."""
    xc = cell_centres()
    v = xc - centre
    r = np.linalg.norm(v, axis=-1)
    u = v / r[..., None]
    pt = centre + radius * u
    n = u if fluid_outside else -u
    band = np.abs(r - radius) <= 3.0 * DX
    n = np.where(band[..., None], n, 0.0)
    pt = np.where(band[..., None], pt, 0.0)
    return n, pt



def nearest_cut_data(kappa, cut, centre, radius, fluid_outside):
    """Per band cell: nearest cut cell's centroid, normal, and (constant) curvature."""
    xc = cell_centres()
    ci = np.argwhere(cut)
    cc = xc[tuple(ci.T)]
    u = cc - centre
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    pt = centre + radius * u
    nn = u if fluid_outside else -u
    dsurf = np.linalg.norm(xc - centre, axis=-1) - radius
    band = np.abs(dsurf) <= 3.0 * DX
    bi = np.argwhere(band)
    bp = xc[tuple(bi.T)]
    d2 = ((bp[:, None, :] - cc[None, :, :]) ** 2).sum(-1)
    near = np.argmin(d2, axis=1)
    N = np.zeros(SHAPE + (3,)); P = np.zeros(SHAPE + (3,))
    N[tuple(bi.T)] = nn[near]; P[tuple(bi.T)] = pt[near]
    return N, P


def quad_reflect(pos, N, P, kap_princ):
    """d, nhat, image from the stored quadratic patch. kap_princ = k1 = k2 (sphere)."""
    a = np.zeros_like(N); a[..., 0] = 1.0
    swap = np.abs(N[..., 0]) > 0.9
    a[swap] = np.array([0.0, 1.0, 0.0])
    t1 = np.cross(N, a); t1 /= np.maximum(np.linalg.norm(t1, axis=-1, keepdims=True), 1e-30)
    t2 = np.cross(N, t1)
    w = pos - P
    xi1 = np.einsum("ij,ij->i", t1, w); xi2 = np.einsum("ij,ij->i", t2, w)
    eta = np.einsum("ij,ij->i", N, w)
    d = eta + 0.5 * kap_princ * (xi1**2 + xi2**2)
    nh = N + kap_princ * (xi1[:, None] * t1 + xi2[:, None] * t2)
    nh /= np.maximum(np.linalg.norm(nh, axis=-1, keepdims=True), 1e-30)
    return d, nh, pos - 2.0 * d[:, None] * nh


def run(radius, fluid_outside, npart=3_000_000, seed=17, margin=3):
    centre = (np.array([8.0 - radius, 8.0, 8.0]) if fluid_outside
              else np.array([8.0 + radius, 8.0, 8.0]))
    kappa, cut = kappa_sphere(centre, radius, fluid_outside, nsub=16)
    mse.KAPPA = kappa
    n_near, p_near, _ = ms.build_cell_planes(kappa, cut, centre, radius, fluid_outside)
    n_own, p_own = own_foot_planes(centre, radius, fluid_outside)
    N, P = nearest_cut_data(kappa, cut, centre, radius, fluid_outside)

    rng = np.random.default_rng(seed)
    side = N_VALID * DX
    t = rng.uniform(0.0, side, size=(npart, 3))
    rt = np.linalg.norm(t - centre, axis=-1)
    pos = t[(rt > radius) if fluid_outside else (rt < radius)]
    mass = side ** 3 / npart
    sign = 1.0 if fluid_outside else -1.0
    kp = sign / radius
    jx = lambda d: np.abs(((1.0 - kp * d) / (1.0 + kp * d)) ** 2)

    v = pos - centre; r = np.linalg.norm(v, axis=-1); u = v / r[:, None]
    variants = {}
    variants["per-particle IF"] = ((r - radius) if fluid_outside else (radius - r),
                                   u if fluid_outside else -u,
                                   centre + (2.0 * radius - r)[:, None] * u,
                                   np.ones(pos.shape[0], bool))
    cell = np.clip(np.floor(pos / DX).astype(np.int64) - LO, 0, np.array(SHAPE) - 1)
    for tag, (nf, pf) in (("nearest-cut plane", (n_near, p_near)),
                          ("own foot plane", (n_own, p_own))):
        nc = nf[cell[:, 0], cell[:, 1], cell[:, 2]]; pc = pf[cell[:, 0], cell[:, 1], cell[:, 2]]
        dd = np.einsum("ij,ij->i", pos - pc, nc)
        variants[tag] = (dd, nc, pos - 2.0 * dd[:, None] * nc,
                         np.linalg.norm(nc, axis=-1) > 0.5)
    nc = N[cell[:, 0], cell[:, 1], cell[:, 2]]; pc = P[cell[:, 0], cell[:, 1], cell[:, 2]]
    okq = np.linalg.norm(nc, axis=-1) > 0.5
    dq, nq, iq = quad_reflect(pos, nc, pc, kp)
    variants["nearest-cut QUADRATIC"] = (dq, nq, iq, okq)

    base = (pos, np.full(pos.shape[0], mass))
    fields = {"plain": mse.deposit([base])}
    order = ["nearest-cut plane", "nearest-cut QUADRATIC", "own foot plane", "per-particle IF"]
    for k in order:
        d, n, img, ok = variants[k]
        sel = ok & (d > 0) & (d <= 1.5 * DX * np.abs(n).sum(-1))
        fields[k] = mse.deposit([base, (img[sel], mass * jx(d[sel]))])

    ax = [np.arange(SHAPE[q]) + LO for q in range(3)]
    ii, jj, kk = np.meshgrid(*ax, indexing="ij"); inter = np.ones(SHAPE, bool)
    for a in (ii, jj, kk): inter &= (a >= margin) & (a <= N_VALID - 1 - margin)
    m = inter & (kappa > 0.0)
    worst = {k: 0.0 for k in fields}
    for lo, hi in BINS:
        s = m & (kappa > lo) & (kappa <= hi)
        if int(s.sum()) < 3: continue
        for k, f in fields.items():
            worst[k] = max(worst[k], abs(float(np.mean(f[s] - 1.0))))
    print(f"  {('convex' if fluid_outside else 'concave')+' R='+str(radius):>16}  " +
          "   ".join(f"{k} {100*worst[k]:5.1f}%" for k in order))
    return worst


if __name__ == "__main__":
    print("### Worst |deviation| over kappa bins, all with the exact Jacobian\n")
    for R in (16.0, 12.0, 8.0, 6.0, 4.0): run(R, True)
    for R in (8.0, 6.0): run(R, False)
