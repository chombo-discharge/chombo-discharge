#!/usr/bin/env python3
"""
Curvature from the DISCRETE normal, not from the level set.

The normal is reconstructed the way EBISBox does it -- from the face area
fractions via the divergence theorem:

    n_hat  =  (alpha_hi,x - alpha_lo,x, ...) / | . |        (pointing into the fluid)

It exists only on cut cells, i.e. on a codimension-1 scatter of points (the
boundary centroids), so grad(n_hat) cannot be taken with a Cartesian stencil.
The shape operator is fitted instead, in the tangent plane of the cell's own
normal:

    u_j = (t1.dx_j, t2.dx_j)      dx_j = centroid_j - centroid_i
    v_j = (t1.dn_j, t2.dn_j)      dn_j = n_j - n_i
    S~  = argmin sum_j | v_j - S~ u_j |^2          (2x2, least squares)
    2H  = tr S~        K = det S~

which needs no principal directions and no eigen-decomposition, exactly like the
level-set route -- only the source of S changes.

Reference: sphere of radius R, dx = 1, so 2H = 2/R and K = 1/R^2 everywhere.
"""
import sys
import numpy as np

DX = 1.0
NSUB = 16          # sub-samples per direction for kappa and face fractions
NSURF = 1_200_000  # surface samples for the boundary centroid


def grid(shape, lo):
    a = [np.arange(shape[d]) + lo for d in range(3)]
    return np.meshgrid(*a, indexing="ij")


def geometry(centre, R, fluid_outside, shape, lo):
    """kappa, discrete normal, cut mask, boundary centroids. Sub-samples cut cells only."""
    ii, jj, kk = grid(shape, lo)
    base = np.stack([ii, jj, kk], -1).astype(float) * DX
    hi = base + DX

    nearest = np.clip(centre, base, hi)
    dmin = np.linalg.norm(nearest - centre, axis=-1)
    far = np.where(np.abs(base - centre) > np.abs(hi - centre), base, hi)
    dmax = np.linalg.norm(far - centre, axis=-1)

    kappa = np.full(shape, np.nan)
    kappa[dmax <= R] = 0.0 if fluid_outside else 1.0
    kappa[dmin >= R] = 1.0 if fluid_outside else 0.0
    cut = np.isnan(kappa)
    kappa[~cut] = kappa[~cut]

    cidx = np.argwhere(cut)
    cbase = base[cut]

    t = (np.arange(NSUB) + 0.5) / NSUB * DX
    ox, oy, oz = np.meshgrid(t, t, t, indexing="ij")
    off = np.stack([ox.ravel(), oy.ravel(), oz.ravel()], -1)
    kv = np.empty(cbase.shape[0])
    for s0 in range(0, cbase.shape[0], 400):
        p = cbase[s0:s0+400, None, :] + off
        r = np.linalg.norm(p - centre, axis=-1)
        kv[s0:s0+400] = ((r > R) if fluid_outside else (r < R)).mean(-1)
    kappa[cut] = kv

    t2 = (np.arange(NSUB) + 0.5) / NSUB * DX
    p2, q2 = np.meshgrid(t2, t2, indexing="ij")
    n = np.zeros(shape + (3,))
    comp = np.empty((cbase.shape[0], 3))
    for d in range(3):
        oth = [x for x in range(3) if x != d]
        vals = []
        for side in (0, 1):
            o = np.zeros((p2.size, 3))
            o[:, oth[0]] = p2.ravel()
            o[:, oth[1]] = q2.ravel()
            o[:, d] = side * DX
            fr = np.empty(cbase.shape[0])
            for s0 in range(0, cbase.shape[0], 800):
                fp = cbase[s0:s0+800, None, :] + o
                rr = np.linalg.norm(fp - centre, axis=-1)
                fr[s0:s0+800] = ((rr > R) if fluid_outside else (rr < R)).mean(-1)
            vals.append(fr)
        comp[:, d] = vals[1] - vals[0]
    nn = np.linalg.norm(comp, axis=-1, keepdims=True)
    n[cut] = np.divide(comp, nn, out=np.zeros_like(comp), where=nn > 1e-12)

    rng = np.random.default_rng(5)
    u = rng.normal(size=(NSURF, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    sp = centre + R * u
    idx = np.floor(sp / DX).astype(np.int64) - lo
    ok = np.all((idx >= 0) & (idx < np.array(shape)), axis=1)
    flat = np.ravel_multi_index(idx[ok].T, shape)
    cnt = np.bincount(flat, minlength=int(np.prod(shape))).astype(float)
    acc = np.stack([np.bincount(flat, weights=sp[ok][:, d], minlength=int(np.prod(shape)))
                    for d in range(3)], -1)
    ctr = np.divide(acc, cnt[:, None], out=np.zeros_like(acc), where=cnt[:, None] > 0)
    ctr = ctr.reshape(shape + (3,))
    have = (cnt.reshape(shape) > 20) & cut
    return kappa, n, cut, ctr, have


def fit_shape_operator(n, ctr, have, shape, radius_cells, symmetrize):
    """Least-squares 2x2 shape operator per cut cell over its cut neighbourhood."""
    idx = np.argwhere(have)
    lut = -np.ones(shape, np.int64)
    lut[have] = np.arange(idx.shape[0])
    two_h = np.full(shape, np.nan)
    kk = np.full(shape, np.nan)
    nnb = np.zeros(shape, np.int64)

    offs = [np.array([a, b, c])
            for a in range(-radius_cells, radius_cells + 1)
            for b in range(-radius_cells, radius_cells + 1)
            for c in range(-radius_cells, radius_cells + 1)
            if (a, b, c) != (0, 0, 0)]

    for cell in idx:
        ni = n[tuple(cell)]
        xi = ctr[tuple(cell)]
        # tangent basis
        a = np.array([1.0, 0.0, 0.0])
        if abs(ni[0]) > 0.9:
            a = np.array([0.0, 1.0, 0.0])
        t1 = np.cross(ni, a)
        t1 /= np.linalg.norm(t1)
        t2 = np.cross(ni, t1)
        U, V = [], []
        for o in offs:
            c2 = cell + o
            if np.any(c2 < 0) or np.any(c2 >= np.array(shape)) or lut[tuple(c2)] < 0:
                continue
            dx = ctr[tuple(c2)] - xi
            dn = n[tuple(c2)] - ni
            U.append([t1 @ dx, t2 @ dx])
            V.append([t1 @ dn, t2 @ dn])
        if len(U) < 4:
            continue
        U = np.array(U)
        V = np.array(V)
        s, *_ = np.linalg.lstsq(U, V, rcond=None)   # U @ s ~ V, so s = S~^T
        S = s.T
        if symmetrize:
            S = 0.5 * (S + S.T)
        two_h[tuple(cell)] = np.trace(S)
        kk[tuple(cell)] = np.linalg.det(S)
        nnb[tuple(cell)] = len(U)
    return two_h, kk, nnb


def report(R, fluid_outside, radius_cells=1, symmetrize=True):
    side = "convex" if fluid_outside else "concave"
    n_valid = max(24, int(2.6 * R))
    pad = 2
    shape = (n_valid + 2 * pad,) * 3
    lo = -pad
    c = np.array([n_valid / 2.0] * 3) * DX
    kappa, nrm, cut, ctr, have = geometry(c, R, fluid_outside, shape, lo)
    two_h, kk, nnb = fit_shape_operator(nrm, ctr, have, shape, radius_cells, symmetrize)

    sgn = 1.0 if fluid_outside else -1.0
    t2h, tk = 2.0 * sgn / R, 1.0 / R ** 2
    m = np.isfinite(two_h)

    # normal accuracy, for reference
    xc = (np.stack(grid(shape, lo), -1).astype(float) + 0.5) * DX
    u = xc - c
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    ntrue = u * sgn
    ang = np.degrees(np.arccos(np.clip(np.sum(nrm[have] * ntrue[have], -1), -1, 1)))

    print(f"\n--- sphere R = {R:g} dx, {side}   stencil radius {radius_cells}"
          f"   {'symmetrized' if symmetrize else 'raw'}"
          f"   cut cells fitted {int(m.sum())}/{int(have.sum())}"
          f"   median neighbours {int(np.median(nnb[m]))}")
    print(f"    discrete normal vs analytic: median {np.median(ang):.3f} deg, "
          f"p95 {np.percentile(ang,95):.3f} deg, max {ang.max():.3f} deg")
    print(f"    {'kappa bin':>12} {'cells':>6} {'2H med':>9} {'2H p05':>9} {'2H p95':>9}"
          f" {'K med':>9} {'K p05':>9} {'K p95':>9}")
    print(f"    {'truth':>12} {'':>6} {t2h:9.5f} {'':>9} {'':>9} {tk:9.5f}")
    for los, his in ((0.0, 0.05), (0.05, 0.15), (0.15, 0.30), (0.30, 0.50),
                     (0.50, 0.75), (0.75, 1.0)):
        s = m & (kappa > los) & (kappa <= his)
        nc = int(s.sum())
        if nc < 4:
            continue
        h_, k_ = two_h[s], kk[s]
        print(f"    {los:5.2f}-{his:5.2f} {nc:6d} {np.median(h_):9.4f} "
              f"{np.percentile(h_,5):9.4f} {np.percentile(h_,95):9.4f} "
              f"{np.median(k_):9.4f} {np.percentile(k_,5):9.4f} {np.percentile(k_,95):9.4f}")
    # how bad is J built from these?
    ds = np.linspace(0.2, 2.6, 25)
    errs = []
    for d in ds:
        jt = abs((1 - t2h * d + tk * d * d) / (1 + t2h * d + tk * d * d))
        num = 1 - two_h[m] * d + kk[m] * d * d
        den = 1 + two_h[m] * d + kk[m] * d * d
        jm = np.abs(num / np.where(np.abs(den) < 1e-9, np.nan, den))
        errs.append(np.nanmedian(np.abs(jm - jt) / jt))
    jlin = [abs(np.maximum(1 - 2 * d * t2h, 0) - abs((1 - t2h*d + tk*d*d)/(1 + t2h*d + tk*d*d)))
            / abs((1 - t2h*d + tk*d*d)/(1 + t2h*d + tk*d*d)) for d in ds]
    print(f"    median |dJ/J| over d in [0.2, 2.6]:  discrete-normal fit "
          f"{np.mean(errs):.1%}   (linearized-with-exact-2H reference {np.mean(jlin):.1%})")
    return two_h, kk, kappa, m


def perturb(n, have, deg, rng):
    """Rotate every discrete normal by a random angle of the given scale."""
    out = n.copy()
    for c in np.argwhere(have):
        v = n[tuple(c)]
        r = rng.normal(size=3)
        r -= (r @ v) * v
        nr = np.linalg.norm(r)
        if nr < 1e-12:
            continue
        a = np.radians(deg) * rng.normal()
        w = v * np.cos(a) + (r / nr) * np.sin(a)
        out[tuple(c)] = w / np.linalg.norm(w)
    return out


def j_error(two_h, kk, m, t2h, tk):
    e = []
    for d in np.linspace(0.2, 2.6, 25):
        jt = abs((1 - t2h * d + tk * d * d) / (1 + t2h * d + tk * d * d))
        num = 1 - two_h[m] * d + kk[m] * d * d
        den = 1 + two_h[m] * d + kk[m] * d * d
        jm = np.abs(num / np.where(np.abs(den) < 1e-9, np.nan, den))
        e.append(np.nanmedian(np.abs(jm - jt) / jt))
    return float(np.mean(e))


def main_noise():
    """How much discrete-normal error can the fit absorb before J is worse than linearized?"""
    print("### Sensitivity of the fitted shape operator to discrete-normal error\n")
    for R, fo in ((6.0, True), (6.0, False)):
        n_valid = max(24, int(2.6 * R))
        pad = 2
        shape = (n_valid + 2 * pad,) * 3
        lo = -pad
        c = np.array([n_valid / 2.0] * 3) * DX
        kappa, nrm, cut, ctr, have = geometry(c, R, fo, shape, lo)
        sgn = 1.0 if fo else -1.0
        t2h, tk = 2.0 * sgn / R, 1.0 / R ** 2
        lin = []
        for d in np.linspace(0.2, 2.6, 25):
            jt = abs((1 - t2h * d + tk * d * d) / (1 + t2h * d + tk * d * d))
            lin.append(abs(max(1 - 2 * d * t2h, 0) - jt) / jt)
        print(f"--- sphere R={R:g}, {'convex' if fo else 'concave'}   "
              f"linearized-with-exact-2H reference: {np.mean(lin):.1%}")
        print(f"    {'normal noise':>13} {'stencil 3^3':>13} {'stencil 5^3':>13} {'stencil 7^3':>13}")
        rng = np.random.default_rng(9)
        for deg in (0.0, 0.5, 1.0, 2.0, 5.0, 10.0):
            pn = nrm if deg == 0 else perturb(nrm, have, deg, rng)
            row = []
            for rad in (1, 2, 3):
                th, kk_, _ = fit_shape_operator(pn, ctr, have, shape, rad, True)
                row.append(f"{j_error(th, kk_, np.isfinite(th), t2h, tk):13.1%}")
            print(f"    {deg:10.1f} deg " + " ".join(row))
        print()


def main_fit():
    print("### Curvature invariants from the DISCRETE (area-fraction) normal")
    print("### Reconstruction and neighbourhood fit only -- the implicit function is")
    print("### used to build the geometry and then never touched again.")
    for R in (8.0, 6.0, 4.0):
        report(R, True, radius_cells=1)
    report(6.0, False, radius_cells=1)
    print("\n\n### Wider stencil (5^3 neighbourhood)")
    for R in (8.0, 6.0, 4.0):
        report(R, True, radius_cells=2)
    report(6.0, False, radius_cells=2)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "fit"
    {"fit": main_fit, "noise": main_noise}[what]()
