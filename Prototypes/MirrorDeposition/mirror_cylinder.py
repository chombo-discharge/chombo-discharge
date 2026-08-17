#!/usr/bin/env python3
"""
Cylindrical EB: does the curvature Jacobian generalize, or was the sphere a fluke?

A cylinder has principal curvatures (1/R, 0). The product form
    J = PROD_i (1 - k_i d)/(1 + k_i d)
therefore predicts J = rho'/rho -- exponent ONE -- whereas the sphere gave
(r'/r)^2. If the exponent were simply D-1 = 2, the cylinder would be
over-corrected. This test discriminates the two.

Everything is Cartesian: uniform cubic cells, tensor-product CIC, cube-cut-by-
cylinder volume fractions. Only the boundary shape is cylindrical.
"""

import math
import numpy as np

N_VALID, N_GHOST, DX = 16, 2, 1.0
PAD = N_GHOST + 4
SHAPE = (N_VALID + 2 * PAD,) * 3
LO = -PAD
AXIS = np.array([8.0, 8.0])          # cylinder axis position in (x, y)


def kappa_cylinder(R, fluid_outside, nsub=400):
    """
    Volume fraction. The cylinder is z-invariant, so kappa of a cell equals the
    2-D area fraction of its (x,y) square -- computed once per column and
    broadcast along z. Corner extrema settle the all-in / all-out columns.
    """
    a = [np.arange(SHAPE[d]) + LO for d in range(2)]
    ii, jj = np.meshgrid(*a, indexing="ij")
    lo = np.stack([ii, jj], axis=-1).astype(float) * DX
    hi = lo + DX

    nearest = np.clip(AXIS, lo, hi)
    dmin = np.linalg.norm(nearest - AXIS, axis=-1)
    far = np.where(np.abs(lo - AXIS) > np.abs(hi - AXIS), lo, hi)
    dmax = np.linalg.norm(far - AXIS, axis=-1)

    col = np.full(lo.shape[:2], np.nan)
    col[dmax <= R] = 0.0 if fluid_outside else 1.0
    col[dmin >= R] = 1.0 if fluid_outside else 0.0

    cut = np.isnan(col)
    if np.any(cut):
        t = (np.arange(nsub) + 0.5) / nsub * DX
        ox, oy = np.meshgrid(t, t, indexing="ij")
        off = np.stack([ox.ravel(), oy.ravel()], axis=-1)
        pts = lo[cut][:, None, :] + off[None, :, :]
        rho = np.linalg.norm(pts - AXIS, axis=-1)
        col[cut] = np.mean(rho > R if fluid_outside else rho < R, axis=-1)

    return np.repeat(col[:, :, None], SHAPE[2], axis=2)


def cic_stencil(pos):
    g = pos / DX - 0.5
    base = np.floor(g).astype(np.int64)
    frac = g - base
    idx = np.empty((pos.shape[0], 8, 3), dtype=np.int64)
    wt = np.ones((pos.shape[0], 8))
    for c, (ox, oy, oz) in enumerate(np.ndindex(2, 2, 2)):
        for dim, off in enumerate((ox, oy, oz)):
            idx[:, c, dim] = base[:, dim] + off
            wt[:, c] *= frac[:, dim] if off else (1.0 - frac[:, dim])
    return idx, wt


def accumulate(field, idx, wt, mass, kappa):
    sh = idx - LO
    ins = np.all((sh >= 0) & (sh < np.array(SHAPE)), axis=-1)
    s = np.clip(sh, 0, np.array(SHAPE) - 1)
    kap = np.zeros(wt.shape)
    kap[ins] = kappa[s[..., 0][ins], s[..., 1][ins], s[..., 2][ins]]
    w = (wt > 0) & ins & (kap > 0)
    fl = np.ravel_multi_index((s[..., 0], s[..., 1], s[..., 2]), SHAPE, mode="clip")
    np.add.at(field.reshape(-1), fl.ravel(),
              np.where(w, wt * mass[:, None], 0.0).ravel())


BINS = [(0.0, 0.05), (0.05, 0.15), (0.15, 0.30), (0.30, 0.50), (0.50, 0.75)]


def run(R, fluid_outside, npart=12_000_000, seed=23, margin=3):
    kappa = kappa_cylinder(R, fluid_outside)
    rng = np.random.default_rng(seed)
    side = N_VALID * DX
    b = rng.uniform(0.0, side, size=(npart, 3))
    rho0 = np.linalg.norm(b[:, :2] - AXIS, axis=-1)
    pos = b[(rho0 > R) if fluid_outside else (rho0 < R)]
    m = side ** 3 / npart

    v = pos[:, :2] - AXIS
    rho = np.linalg.norm(v, axis=-1)
    u = v / rho[:, None]
    d = (rho - R) if fluid_outside else (R - rho)
    rho_img = 2 * R - rho
    img = pos.copy()
    img[:, :2] = AXIS + rho_img[:, None] * u          # z untouched: no curvature there

    nhat = np.zeros_like(pos)
    nhat[:, :2] = u if fluid_outside else -u
    band = DX * (1.0 + np.abs(nhat).sum(-1) / 2)
    sel = (d > 0) & (d <= band)

    ratio = np.abs(rho_img) / rho
    variants = {
        "no J":            np.ones(pos.shape[0]),
        "J = ratio^1":     ratio,           # product form: k = (1/R, 0)
        "J = ratio^2":     ratio ** 2,      # wrong: assumes D-1 like the sphere
    }

    fields = {}
    base_idx, base_wt = cic_stencil(pos)
    img_idx, img_wt = cic_stencil(img[sel])
    plain = np.zeros(SHAPE)
    accumulate(plain, base_idx, base_wt, np.full(pos.shape[0], m), kappa)
    fields["plain"] = plain / DX ** 3
    for label, J in variants.items():
        f = plain.copy()
        accumulate(f, img_idx, img_wt, m * J[sel], kappa)
        fields[label] = f / DX ** 3

    a = [np.arange(SHAPE[q]) + LO for q in range(3)]
    ii, jj, kk = np.meshgrid(*a, indexing="ij")
    interior = np.ones(SHAPE, dtype=bool)
    for q in (ii, jj):
        interior &= (q >= margin) & (q <= N_VALID - 1 - margin)
    interior &= (kk >= margin) & (kk <= N_VALID - 1 - margin)
    mask = interior & (kappa > 0.0)

    print(f"\n=== cylinder R = {R} dx, fluid {'outside' if fluid_outside else 'inside'}"
          f"   particles {pos.shape[0]}   reflected {100*sel.mean():.1f}%")
    print(f"    {'kappa bin':>14} {'cells':>7} {'plain':>10} {'no J':>10}"
          f" {'J=ratio^1':>12} {'J=ratio^2':>12}")
    for lo_k, hi_k in BINS:
        sel_k = mask & (kappa > lo_k) & (kappa <= hi_k)
        n = int(sel_k.sum())
        if n < 5:
            continue
        row = [np.mean(fields[k][sel_k] - 1.0)
               for k in ("plain", "no J", "J = ratio^1", "J = ratio^2")]
        print(f"    {lo_k:5.3f}-{hi_k:5.3f} {n:7d} {row[0]:+10.4f} {row[1]:+10.4f}"
              f" {row[2]:+12.4f} {row[3]:+12.4f}")


if __name__ == "__main__":
    # Cartesian check of the Jacobian eigenvalues for the cylindrical map
    R = 5.0
    def T(x):
        v = x[:2] - AXIS
        r = np.linalg.norm(v)
        out = x.copy()
        out[:2] = AXIS + (2 * R - r) * v / r
        return out
    rng = np.random.default_rng(1)
    print("Cartesian Jacobian of the cylindrical reflection (finite differences)")
    for _ in range(3):
        th = rng.uniform(0, 2 * np.pi)
        r = rng.uniform(R + 0.05, R + 1.7)
        x = np.array([AXIS[0] + r * np.cos(th), AXIS[1] + r * np.sin(th),
                      rng.uniform(2, 14)])
        h = 1e-6
        J = np.stack([(T(x + h * e) - T(x - h * e)) / (2 * h)
                      for e in np.eye(3)], axis=1)
        ev = np.sort_complex(np.linalg.eigvals(J)).real
        print(f"   rho={r:6.3f}   eigenvalues {ev.round(5)}   |det| = "
              f"{abs(np.linalg.det(J)):.6f}   rho'/rho = {(2*R-r)/r:.6f}")

    for R in (12.0, 8.0, 5.0):
        run(R, True)
    run(6.0, False, margin=1)
