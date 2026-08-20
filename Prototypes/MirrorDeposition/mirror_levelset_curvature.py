#!/usr/bin/env python3
"""
Claim 3 at DEPOSITION level, and with a self-consistent reference.

Two questions the continuum harness cannot answer:

  (1) Does the invariant route reproduce the SAME IF's own |det grad R|?  That is
      the recipe-fidelity question, and it is answerable for a faceted or noisy
      surface, where "the analytic sphere's J" is not the right target at all.

  (2) With (2H, K) built by differencing at a realisable step h, does the exact
      Jacobian still beat the linearized one IN THE DEPOSIT?

dx = 1, sphere R = 6 dx.
"""
import sys

import numpy as np
sys.path.insert(0, "/home/robertm/Projects/chombo-discharge/Prototypes/MirrorDeposition")
import mirror_sphere as ms
import mirror_sphere_ext as mse
from mirror_sphere import BINS, DX, LO, N_VALID, SHAPE, kappa_sphere

R = 6.0
STEPS = [1e-4, 1e-2, 0.05, 0.1, 0.25, 0.5, 1.0]


def make_ifs(centre, sign):
    def sph(x):
        return sign * (np.linalg.norm(x - centre, axis=-1) - R)

    def sph_f32(x):
        v = np.asarray(x, np.float32) - np.float32(centre)
        return (np.float32(sign) * (np.sqrt((v * v).sum(-1)) - np.float32(R))).astype(np.float64)

    def sph_quad(x):
        v = x - centre
        return sign * ((v * v).sum(-1) - R * R)

    def noisy(x):
        v = x - centre
        k = 2.0 * np.pi / 2.0
        return sign * (np.linalg.norm(v, axis=-1) - R
                       - 0.05 * np.sin(k * v[..., 0]) * np.sin(k * v[..., 1]) * np.sin(k * v[..., 2]))

    return {"analytic SDF (double)": sph, "analytic SDF (FLOAT32)": sph_f32,
            "non-SDF r^2-R^2": sph_quad, "SDF + 0.05dx noise, lam=2dx": noisy}


def vgrad(f, x, h):
    g = np.empty_like(x)
    for i in range(3):
        e = np.zeros(3)
        e[i] = h
        g[..., i] = (f(x + e) - f(x - e)) / (2.0 * h)
    return g


def vnhat(f, x, h):
    g = vgrad(f, x, h)
    return g / np.linalg.norm(g, axis=-1, keepdims=True)


def vinv(f, p, h):
    """(2H, K) at points p, central differences of nhat with step h."""
    cols = []
    for i in range(3):
        e = np.zeros(3)
        e[i] = h
        cols.append((vnhat(f, p + e, h) - vnhat(f, p - e, h)) / (2.0 * h))
    m = np.stack(cols, axis=-1)                       # (...,3,3)  M[:,i,j] = d n_i / d x_j
    tr = np.trace(m, axis1=-2, axis2=-1)
    tr2 = np.trace(m @ m, axis1=-2, axis2=-1)
    return tr, 0.5 * (tr * tr - tr2)


def vdet(f, x, h):
    """|det grad R| of the reflection map built from this IF. The ground truth."""
    def refl(y):
        g = vgrad(f, y, h)
        gn = np.linalg.norm(g, axis=-1, keepdims=True)
        return y - 2.0 * (f(y)[..., None] / gn) * (g / gn)

    cols = []
    for i in range(3):
        e = np.zeros(3)
        e[i] = h
        cols.append((refl(x + e) - refl(x - e)) / (2.0 * h))
    return np.abs(np.linalg.det(np.stack(cols, axis=-1)))


def fidelity():
    print("### (1) Recipe fidelity: J from (2H,K) at step h  vs  |det grad R| of the SAME IF")
    print("###     Both from the same implicit function, so this isolates the recipe.\n")
    centre = np.zeros(3)
    rng = np.random.default_rng(3)
    u = rng.normal(size=(400, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    d = rng.uniform(0.2, 2.6, 400)
    pts = centre + (R + d)[:, None] * u
    for name, f in make_ifs(centre, 1.0).items():
        print(f"--- {name}")
        print(f"    {'h/dx':>7} {'med 2H':>9} {'med K':>9} {'med |dJ/J|':>11} {'p90':>9} {'worst':>10}")
        for h in STEPS:
            fv = f(pts)
            g = vgrad(f, pts, h)
            gn = np.linalg.norm(g, axis=-1, keepdims=True)
            foot = pts - (fv[:, None] / gn) * (g / gn)
            dd = np.abs(fv) / gn[:, 0]
            two_h, k = vinv(f, foot, h)
            num, den = 1.0 - two_h * dd + k * dd**2, 1.0 + two_h * dd + k * dd**2
            jm = np.abs(num / np.where(np.abs(den) < 1e-9, np.nan, den))
            jt = vdet(f, pts, h)
            ok = np.isfinite(jm) & np.isfinite(jt) & (jt > 1e-6)
            e = np.abs(jm[ok] - jt[ok]) / jt[ok]
            print(f"    {h:7.0e} {np.median(two_h):9.5f} {np.median(k):9.5f} "
                  f"{np.median(e):11.2%} {np.percentile(e,90):9.2%} {e.max():10.2%}")
        print()


def deposit_test(fluid_outside, npart=3_000_000, seed=17, margin=3):
    sign = 1.0 if fluid_outside else -1.0
    centre = (np.array([8.0 - R, 8.0, 8.0]) if fluid_outside else np.array([8.0 + R, 8.0, 8.0]))
    kappa, _ = kappa_sphere(centre, R, fluid_outside, nsub=16)
    mse.KAPPA = kappa
    rng = np.random.default_rng(seed)
    side = N_VALID * DX
    trial = rng.uniform(0.0, side, size=(npart, 3))
    rt = np.linalg.norm(trial - centre, axis=-1)
    pos = trial[(rt > R) if fluid_outside else (rt < R)]
    mass = side**3 / npart

    vec = pos - centre
    r = np.linalg.norm(vec, axis=-1)
    unit = vec / r[:, None]
    d = (r - R) if fluid_outside else (R - r)
    img = centre + (2.0 * R - r)[:, None] * unit
    n = unit if fluid_outside else -unit
    sel = (d > 0) & (d <= 1.5 * DX * np.abs(n).sum(-1))
    base = (pos, np.full(pos.shape[0], mass))
    foot = centre + R * unit[sel]
    ds = d[sel]

    axes = [np.arange(SHAPE[q]) + LO for q in range(3)]
    ii, jj, kk = np.meshgrid(*axes, indexing="ij")
    interior = np.ones(SHAPE, bool)
    for a in (ii, jj, kk):
        interior &= (a >= margin) & (a <= N_VALID - 1 - margin)
    mask = interior & (kappa > 0.0)

    def worst(field):
        w = 0.0
        for lo, hi in BINS:
            s = mask & (kappa > lo) & (kappa <= hi)
            if int(s.sum()) >= 3:
                w = max(w, abs(float(np.mean(field[s] - 1.0))))
        return w

    side_s = "convex" if fluid_outside else "concave"
    print(f"\n--- sphere R = {R} dx, {side_s}: worst |dev| over kappa bins")
    tha, ka = 2.0 * sign / R, 1.0 / R**2
    rows = [("no Jacobian", np.ones(ds.shape)),
            ("linearized, ANALYTIC 2H", np.maximum(1.0 - 2.0 * ds * tha, 0.0)),
            ("exact, ANALYTIC 2H,K", np.abs((1 - tha*ds + ka*ds**2) / (1 + tha*ds + ka*ds**2)))]
    for label, jv in rows:
        print(f"    {label:<40} {100*worst(mse.deposit([base, (img[sel], mass*jv)])):6.2f}%")
    for name, f in make_ifs(centre, sign).items():
        for h in (1e-4, 0.1, 0.5, 1.0):
            two_h, k = vinv(f, foot, h)
            den = 1.0 + two_h * ds + k * ds**2
            jv = np.abs((1.0 - two_h * ds + k * ds**2) / den)
            jv = np.where(np.isfinite(jv) & (np.abs(den) > 1e-6), jv, 1.0)
            jv = np.clip(jv, 0.0, 50.0)
            lab = f"exact, MEASURED from {name}, h={h}"
            print(f"    {lab:<40} {100*worst(mse.deposit([base, (img[sel], mass*jv)])):6.2f}%")


def endtoend(fluid_outside, ifname, h, npart=3_000_000, seed=17, margin=3):
    """d, the image position AND (2H,K) all taken from the IF -- nothing analytic."""
    sign = 1.0 if fluid_outside else -1.0
    centre = np.array([8.0 - R, 8.0, 8.0]) if fluid_outside else np.array([8.0 + R, 8.0, 8.0])
    kappa, _ = kappa_sphere(centre, R, fluid_outside, nsub=16)
    mse.KAPPA = kappa
    rng = np.random.default_rng(seed)
    side = N_VALID * DX
    t = rng.uniform(0, side, size=(npart, 3))
    rt = np.linalg.norm(t - centre, axis=-1)
    pos = t[(rt > R) if fluid_outside else (rt < R)]
    mass = side ** 3 / npart
    f = make_ifs(centre, sign)[ifname]
    fv = f(pos)
    g = vgrad(f, pos, h)
    gn = np.linalg.norm(g, axis=-1, keepdims=True)
    n = g / gn
    d = fv / gn[:, 0]
    img = pos - 2.0 * d[:, None] * n
    sel = (d > 0) & (d <= 1.5 * DX * np.abs(n).sum(-1))
    foot = pos[sel] - d[sel, None] * n[sel]
    th, k = vinv(f, foot, h)
    ds = d[sel]
    den = 1.0 + th * ds + k * ds * ds
    nsing = int(np.sum(np.abs(den) < 0.05))
    jv = np.abs((1.0 - th * ds + k * ds * ds) / np.where(np.abs(den) < 1e-9, np.nan, den))
    jv = np.clip(np.where(np.isfinite(jv), jv, 1.0), 0.0, 50.0)
    base = (pos, np.full(pos.shape[0], mass))
    fld = mse.deposit([base, (img[sel], mass * jv)])
    ax = [np.arange(SHAPE[q]) + LO for q in range(3)]
    ii, jj, kk = np.meshgrid(*ax, indexing="ij")
    inter = np.ones(SHAPE, bool)
    for a in (ii, jj, kk):
        inter &= (a >= margin) & (a <= N_VALID - 1 - margin)
    m = inter & (kappa > 0)
    w = max(abs(float(np.mean(fld[s] - 1.0))) for lo, hi in BINS
            if int((s := (m & (kappa > lo) & (kappa <= hi))).sum()) >= 3)
    tot = pos.shape[0] * mass
    print(f"  {'convex' if fluid_outside else 'concave':>7}  {ifname:<28} h={h:<7} "
          f"worst {100 * w:6.2f}%   mass {float(np.sum(fld * kappa) * DX ** 3) / tot:.5f}   "
          f"|den|<0.05 for {100 * nsing / max(1, ds.size):.3f}% of images  (Jmax {jv.max():.1f})")


def main_endtoend():
    print("### End to end: d, image position AND (2H,K) all from the IF, nothing analytic\n")
    for fo in (True, False):
        for nm in ("analytic SDF (double)", "non-SDF r^2-R^2", "SDF + 0.05dx noise, lam=2dx"):
            for h in (1e-4, 0.1, 1.0):
                endtoend(fo, nm, h)
        print()


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "endtoend"):
        main_endtoend()
    if what not in ("all", "fidelity", "deposit"):
        raise SystemExit(0)
    fidelity()
    print("\n### (2) Deposition-level: does a MEASURED (2H,K) still beat the linearized form?\n")
    deposit_test(True)
    deposit_test(False)
