#!/usr/bin/env python3
"""
Curved-surface test for mirrored deposition: does the CELL-CONSTANT PLANE hold up?

Everything measured so far used a planar EB, where every cut cell trivially agrees
on one plane, so the plane approximation was never exercised. Here the EB is a
sphere and the three reflections diverge:

  plain    no mirror at all (baseline)
  cell     each band cell stores ONE plane, taken from the NEAREST CUT CELL's
           normal and boundary centroid -- this is the algorithm we plan to ship
  exact    reflect radially through the sphere: r -> 2R - r along the same ray

'exact' is the ideal the mirror is defined against. 'cell' minus 'exact' is the
error the cell-constant plane introduces, which is the open question.

SCOPE: the normals here are analytic (radial at the projected boundary point).
A real EBISBox normal is reconstructed from discrete area fractions and is noisier,
worst as kappa -> 0. This harness therefore isolates the GEOMETRIC approximation
of a per-cell plane; it does not measure discrete-normal reconstruction error.
That one can only be measured in the code.
"""

import math
import numpy as np

N_VALID, N_GHOST, DX = 16, 2, 1.0
PAD = N_GHOST + 4
SHAPE = (N_VALID + 2 * PAD,) * 3
LO = -PAD


def cell_centres():
    a = [np.arange(SHAPE[d]) + LO for d in range(3)]
    ii, jj, kk = np.meshgrid(*a, indexing="ij")
    return (np.stack([ii, jj, kk], axis=-1).astype(float) + 0.5) * DX


def kappa_sphere(centre, R, fluid_outside, nsub=16):
    """Exact-ish fluid fraction. Corner extrema decide all-in / all-out cells."""
    a = [np.arange(SHAPE[d]) + LO for d in range(3)]
    ii, jj, kk = np.meshgrid(*a, indexing="ij")
    lo = np.stack([ii, jj, kk], axis=-1).astype(float) * DX
    hi = lo + DX

    # min and max |x - centre| over the cell
    nearest = np.clip(centre, lo, hi)
    dmin = np.linalg.norm(nearest - centre, axis=-1)
    far = np.where(np.abs(lo - centre) > np.abs(hi - centre), lo, hi)
    dmax = np.linalg.norm(far - centre, axis=-1)

    kap = np.full(SHAPE, np.nan)
    inside_all = dmax <= R
    outside_all = dmin >= R
    kap[inside_all] = 0.0 if fluid_outside else 1.0
    kap[outside_all] = 1.0 if fluid_outside else 0.0

    cut = np.isnan(kap)
    if np.any(cut):
        t = (np.arange(nsub) + 0.5) / nsub * DX
        ox, oy, oz = np.meshgrid(t, t, t, indexing="ij")
        off = np.stack([ox.ravel(), oy.ravel(), oz.ravel()], axis=-1)
        base = lo[cut]                                    # (Ncut, 3)
        pts = base[:, None, :] + off[None, :, :]          # (Ncut, nsub^3, 3)
        r = np.linalg.norm(pts - centre, axis=-1)
        frac = np.mean(r > R if fluid_outside else r < R, axis=-1)
        kap[cut] = frac
    return kap, cut


def build_cell_planes(kappa, cut, centre, R, fluid_outside):
    """
    Per band cell: the plane of the NEAREST CUT CELL. Its normal is the analytic
    surface normal at that cut cell's boundary point, and the plane passes through
    that point. Returns (nhat, point) per cell, valid where band is True.
    """
    xc = cell_centres()
    cut_idx = np.argwhere(cut)
    cut_ctr = xc[cut[..., ] if False else tuple(cut_idx.T)]      # (Ncut, 3)

    u = cut_ctr - centre
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    cut_pt = centre + R * u                                       # boundary point
    cut_n = u if fluid_outside else -u                            # into the fluid

    # band = within 3 cells of the surface, generous
    dsurf = np.linalg.norm(xc - centre, axis=-1) - R
    band = np.abs(dsurf) <= 3.0 * DX

    band_idx = np.argwhere(band)
    bp = xc[tuple(band_idx.T)]                                    # (Nband, 3)
    # nearest cut cell by centre-to-centre distance
    d2 = ((bp[:, None, :] - cut_ctr[None, :, :]) ** 2).sum(-1)
    near = np.argmin(d2, axis=1)

    nhat = np.zeros(SHAPE + (3,))
    point = np.zeros(SHAPE + (3,))
    nhat[tuple(band_idx.T)] = cut_n[near]
    point[tuple(band_idx.T)] = cut_pt[near]
    return nhat, point, band


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


def deposit(pos, m, kappa, images=None):
    field = np.zeros(SHAPE)
    sets = [cic_stencil(pos)]
    if images is not None and images.shape[0]:
        sets.append(cic_stencil(images))
    for idx, wt in sets:
        sh = idx - LO
        ins = np.all((sh >= 0) & (sh < np.array(SHAPE)), axis=-1)
        s = np.clip(sh, 0, np.array(SHAPE) - 1)
        kap = np.zeros(wt.shape)
        kap[ins] = kappa[s[..., 0][ins], s[..., 1][ins], s[..., 2][ins]]
        w = (wt > 0) & ins & (kap > 0)
        fl = np.ravel_multi_index((s[..., 0], s[..., 1], s[..., 2]), SHAPE, mode="clip")
        np.add.at(field.reshape(-1), fl.ravel(), np.where(w, wt * m, 0.0).ravel())
    return field / DX ** 3


BINS = [(0.0, 0.05), (0.05, 0.15), (0.15, 0.30), (0.30, 0.50),
        (0.50, 0.75), (0.75, 0.999)]


def run(name, centre, R, fluid_outside, npart=3_000_000, seed=17, margin=3):
    kappa, cut = kappa_sphere(centre, R, fluid_outside)
    nhat, point, band = build_cell_planes(kappa, cut, centre, R, fluid_outside)
    xc = cell_centres()

    rng = np.random.default_rng(seed)
    side = N_VALID * DX
    Vbox = side ** 3

    # Each trial point represents Vbox/n_try of volume, so the kept (fluid) points
    # carry m = Vbox/n_try apiece and the density is exactly 1 by construction.
    # Deriving the fluid volume this way is robust to the sphere hanging outside
    # the box -- Vbox - Vsphere is only correct when it lies wholly inside, and
    # is wildly negative for a large radius.
    n_try = npart
    b = rng.uniform(0.0, side, size=(n_try, 3))
    r = np.linalg.norm(b - centre, axis=-1)
    pos = b[(r > R) if fluid_outside else (r < R)]
    m = Vbox / n_try
    if pos.shape[0] < 1000:
        print(f"\n=== {name}: only {pos.shape[0]} fluid particles, skipping")
        return

    # --- exact radial reflection
    v = pos - centre
    r = np.linalg.norm(v, axis=-1, keepdims=True)
    u = v / r
    d_exact = (r[:, 0] - R) if fluid_outside else (R - r[:, 0])
    img_exact = centre + (2 * R - r) * u

    # --- cell-constant plane, from the nearest cut cell
    pc = np.floor(pos / DX).astype(np.int64) - LO
    pc = np.clip(pc, 0, np.array(SHAPE) - 1)
    n_c = nhat[pc[:, 0], pc[:, 1], pc[:, 2]]
    p_c = point[pc[:, 0], pc[:, 1], pc[:, 2]]
    has = np.linalg.norm(n_c, axis=-1) > 0.5
    d_cell = np.einsum("ij,ij->i", pos - p_c, n_c)
    img_cell = pos - 2.0 * d_cell[:, None] * n_c

    banded = lambda d, n: (d > 0) & (d <= DX * (1.0 + np.abs(n).sum(-1) / 2))
    sel_exact = banded(d_exact, u)
    sel_cell = has & banded(d_cell, n_c)

    fields = {
        "plain": deposit(pos, m, kappa),
        "cell":  deposit(pos, m, kappa, img_cell[sel_cell]),
        "exact": deposit(pos, m, kappa, img_exact[sel_exact]),
    }

    a = [np.arange(SHAPE[d]) + LO for d in range(3)]
    ii, jj, kk = np.meshgrid(*a, indexing="ij")
    interior = np.ones(SHAPE, dtype=bool)
    for q in (ii, jj, kk):
        interior &= (q >= margin) & (q <= N_VALID - 1 - margin)
    mask = interior & (kappa > 0.0)

    print(f"\n=== {name}   R = {R} dx   fluid {'outside' if fluid_outside else 'inside'}"
          f"   particles {pos.shape[0]}   reflected: cell {100*sel_cell.mean():.1f}%"
          f"  exact {100*sel_exact.mean():.1f}%")
    print(f"    {'kappa bin':>14} {'cells':>6} {'plain':>10} {'cell-plane':>13}"
          f" {'exact':>11} {'cell - exact':>14}")
    for lo_k, hi_k in BINS:
        sel = mask & (kappa > lo_k) & (kappa <= hi_k)
        n = int(sel.sum())
        if n < 3:
            continue
        dv = {k: f[sel] - 1.0 for k, f in fields.items()}
        sem = np.std(dv["cell"], ddof=1) / math.sqrt(n)
        print(f"    {lo_k:5.3f}-{hi_k:5.3f} {n:6d} {dv['plain'].mean():+10.4f}"
              f" {dv['cell'].mean():+10.4f}+-{sem:4.4f} {dv['exact'].mean():+11.4f}"
              f" {dv['cell'].mean()-dv['exact'].mean():+14.4f}")
    tot = pos.shape[0] * m
    for k, f in fields.items():
        print(f"    mass {k:>6}: sum(kappa phi dV)/sum(m) = "
              f"{float(np.sum(f*kappa)*DX**3)/tot:.5f}")


if __name__ == "__main__":
    run("convex sphere (electrode)", np.array([8., 8., 8.]), 4.0, True)
    run("convex sphere, coarser",    np.array([8., 8., 8.]), 2.5, True)
    run("concave cavity",            np.array([8., 8., 8.]), 6.0, False, margin=1)
