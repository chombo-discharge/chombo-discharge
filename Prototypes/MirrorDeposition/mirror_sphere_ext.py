#!/usr/bin/env python3
"""
Extended sphere harness: the Jacobian, the plane source, and what the mass check
actually measures.

mirror_sphere.py compares three reflections (none / cell-constant plane / exact
radial) at three radii and applies no Jacobian at all. That is not enough to
choose a design, because the two error sources it mixes are of comparable size
and only one of them is visible while the other is uncorrected. This harness
crosses them:

    reflection source   plain, cell-constant plane, exact radial
    Jacobian            none, linearized (1 - 2 d div n), exact product form

and sweeps radius, curvature sign, kappa sub-sampling and scoring margin.

Three things it establishes that mirror_sphere.py cannot:

1. The Jacobian is needed, and it is needed independently of the plane. With
   EXACT reflection the uncorrected error is still +40% at R = 8 dx and +117%
   at R = 4 dx, so no improvement to the plane removes it.

2. The cell-constant plane is NOT a minor term. mirror_sphere.py's conclusion
   that it is was drawn from the kappa <= 0.05 column while the uncorrected
   curvature error ran to +1.17 and drowned it. Once the Jacobian is applied
   the ordering reverses below R ~ 12 dx.

3. sum(kappa phi dV)/sum(m) is an IDENTITY, not an independent check:

       mass ratio == mean_p [ kappatilde(x_p) + J_p * kappatilde(R(x_p)) ]

   where kappatilde is the CIC interpolant of the kappa field. This holds to
   machine precision (verified by `identity` below). Two consequences: on a
   plane it is 1 by construction and proves nothing, and on a curved surface it
   is the kappa-weighted volume average of the same density error the binned
   tables report. It is NOT a quadrature artifact of sub-sampled kappa -- the
   `nsub` sweep moves it by <= 3e-4 while the discrepancy is 1e-2 to 1e-1.

Usage:
    python3 mirror_sphere_ext.py radii      # the main cross, convex and concave
    python3 mirror_sphere_ext.py nsub       # is the mass ratio a kappa artifact?
    python3 mirror_sphere_ext.py identity   # verify the mass identity
    python3 mirror_sphere_ext.py margin     # scoring-margin sensitivity
"""

import sys

import numpy as np

from mirror_sphere import (
    BINS,
    DX,
    LO,
    N_VALID,
    SHAPE,
    build_cell_planes,
    cic_stencil,
    kappa_sphere,
)


# --------------------------------------------------------------------------- #
# Deposition
# --------------------------------------------------------------------------- #
def deposit(sets):
    """sets = [(positions, per-particle mass), ...]. Covered cells are skipped."""
    field = np.zeros(SHAPE)
    for pos, mass in sets:
        if pos.shape[0] == 0:
            continue
        idx, wt = cic_stencil(pos)
        sh = idx - LO
        ins = np.all((sh >= 0) & (sh < np.array(SHAPE)), axis=-1)
        s = np.clip(sh, 0, np.array(SHAPE) - 1)
        kap = np.zeros(wt.shape)
        kap[ins] = KAPPA[s[..., 0][ins], s[..., 1][ins], s[..., 2][ins]]
        live = (wt > 0) & ins & (kap > 0)
        flat = np.ravel_multi_index((s[..., 0], s[..., 1], s[..., 2]), SHAPE, mode="clip")
        np.add.at(field.reshape(-1), flat.ravel(), np.where(live, wt * mass[:, None], 0.0).ravel())
    return field / DX**3


def interpolated_kappa(pos):
    """sum_j kappa_j W(x_j - x): the CIC interpolant of kappa, evaluated at pos."""
    idx, wt = cic_stencil(pos)
    sh = idx - LO
    ins = np.all((sh >= 0) & (sh < np.array(SHAPE)), axis=-1)
    s = np.clip(sh, 0, np.array(SHAPE) - 1)
    kap = np.zeros(wt.shape)
    kap[ins] = KAPPA[s[..., 0][ins], s[..., 1][ins], s[..., 2][ins]]
    return np.sum(kap * wt * ins, axis=1)


KAPPA = None  # module-level so deposit() and interpolated_kappa() see one field


# --------------------------------------------------------------------------- #
# One configuration
# --------------------------------------------------------------------------- #
def build(radius, fluid_outside, npart, seed, nsub):
    """
    Places the surface near x = 8 whatever the radius is, so R = 40 or 200 are
    runnable inside a 16^3 box. mirror_sphere.py centres the sphere and so can
    only run radii that fit, which is why its table stops at R = 4.
    """
    global KAPPA
    centre = np.array([8.0 - radius, 8.0, 8.0]) if fluid_outside else np.array(
        [8.0 + radius, 8.0, 8.0]
    )
    KAPPA, cut = kappa_sphere(centre, radius, fluid_outside, nsub=nsub)
    normals, points, _ = build_cell_planes(KAPPA, cut, centre, radius, fluid_outside)

    rng = np.random.default_rng(seed)
    side = N_VALID * DX
    trial = rng.uniform(0.0, side, size=(npart, 3))
    r_trial = np.linalg.norm(trial - centre, axis=-1)
    pos = trial[(r_trial > radius) if fluid_outside else (r_trial < radius)]
    mass = side**3 / npart  # each trial point carries Vbox/npart, so density == 1

    # --- exact radial reflection
    vec = pos - centre
    r = np.linalg.norm(vec, axis=-1)
    unit = vec / r[:, None]
    d_exact = (r - radius) if fluid_outside else (radius - r)
    img_exact = centre + (2.0 * radius - r)[:, None] * unit
    n_exact = unit if fluid_outside else -unit

    # --- cell-constant plane, taken from the nearest cut cell
    cell = np.clip(np.floor(pos / DX).astype(np.int64) - LO, 0, np.array(SHAPE) - 1)
    n_cell = normals[cell[:, 0], cell[:, 1], cell[:, 2]]
    p_cell = points[cell[:, 0], cell[:, 1], cell[:, 2]]
    has_plane = np.linalg.norm(n_cell, axis=-1) > 0.5
    d_cell = np.einsum("ij,ij->i", pos - p_cell, n_cell)
    img_cell = pos - 2.0 * d_cell[:, None] * n_cell

    # The corrected band: (3/2)*dx*sum|n_i|, not dx*(1 + sum|n_i|/2). See
    # mirror_planar_ext.py band, which measures why.
    def in_band(d, n):
        return (d > 0) & (d <= 1.5 * DX * np.abs(n).sum(-1))

    sel_exact = in_band(d_exact, n_exact)
    sel_cell = has_plane & in_band(d_cell, n_cell)

    # Curvature is a property of the surface, so the same principal curvatures
    # apply to both reflections; only the distance differs.
    sign = 1.0 if fluid_outside else -1.0
    two_h = 2.0 * sign / radius

    def jac_exact(d):
        return abs(((1.0 - sign * d / radius) / (1.0 + sign * d / radius)) ** 2)

    def jac_linear(d):
        return np.maximum(1.0 - 2.0 * d * two_h, 0.0)

    return dict(
        pos=pos,
        mass=mass,
        img_exact=img_exact,
        sel_exact=sel_exact,
        d_exact=d_exact,
        img_cell=img_cell,
        sel_cell=sel_cell,
        d_cell=d_cell,
        jac_exact=jac_exact,
        jac_linear=jac_linear,
    )


VARIANTS = [
    ("plain", None, None),
    ("cell, no J", "cell", None),
    ("cell + Jlin", "cell", "linear"),
    ("cell + Jexact", "cell", "exact"),
    ("exact, no J", "exact", None),
    ("exact + Jlin", "exact", "linear"),
    ("exact + Jexact", "exact", "exact"),
]


def fields_for(cfg, which=None):
    pos, mass = cfg["pos"], cfg["mass"]
    base = (pos, np.full(pos.shape[0], mass))
    out = {}
    for label, source, jac in VARIANTS:
        if which is not None and label not in which:
            continue
        if source is None:
            out[label] = deposit([base])
            continue
        img, sel, dist = cfg[f"img_{source}"], cfg[f"sel_{source}"], cfg[f"d_{source}"]
        if jac is None:
            w = np.full(int(sel.sum()), mass)
        else:
            w = mass * cfg[f"jac_{jac}"](dist[sel])
        out[label] = deposit([base, (img[sel], w)])
    return out


def score_mask(margin):
    axes = [np.arange(SHAPE[q]) + LO for q in range(3)]
    ii, jj, kk = np.meshgrid(*axes, indexing="ij")
    interior = np.ones(SHAPE, dtype=bool)
    for a in (ii, jj, kk):
        interior &= (a >= margin) & (a <= N_VALID - 1 - margin)
    return interior & (KAPPA > 0.0)


# --------------------------------------------------------------------------- #
# Drivers
# --------------------------------------------------------------------------- #
def run(radius, fluid_outside, npart=3_000_000, seed=17, margin=3, nsub=16, which=None):
    cfg = build(radius, fluid_outside, npart, seed, nsub)
    fields = fields_for(cfg, which)
    mask = score_mask(margin)
    total = cfg["pos"].shape[0] * cfg["mass"]

    side = "convex" if fluid_outside else "concave"
    print(
        f"\n=== sphere R = {radius} dx, {side}   particles {cfg['pos'].shape[0]}"
        f"   reflected {100 * cfg['sel_exact'].mean():.1f}%   margin {margin}   nsub {nsub}"
    )
    header = " ".join(f"{k:>14}" for k in fields)
    print(f"    {'kappa bin':>13} {'cells':>6} {header}")
    for lo, hi in BINS:
        sel = mask & (KAPPA > lo) & (KAPPA <= hi)
        n = int(sel.sum())
        if n < 3:
            continue
        row = " ".join(f"{np.mean(f[sel] - 1.0):+14.4f}" for f in fields.values())
        print(f"    {lo:5.3f}-{hi:5.3f} {n:6d} {row}")
    masses = {k: float(np.sum(f * KAPPA) * DX**3) / total for k, f in fields.items()}
    print("    mass  " + "   ".join(f"{k} {v:.5f}" for k, v in masses.items()))
    return cfg, fields, masses


def main_radii():
    print(
        "### The cross. Read two things off it:\n"
        "###   'exact, no J' vs 'exact + Jexact'  -- the Jacobian is needed on its own\n"
        "###   'cell + Jexact' vs 'exact + Jexact' -- the plane is NOT a minor term\n"
    )
    for radius in (40.0, 16.0, 12.0, 8.0, 6.0, 4.0):
        run(radius, True)
    for radius in (8.0, 6.0):
        run(radius, False)


def main_nsub():
    print(
        "### Is the curved-geometry mass ratio a quadrature artifact of sub-sampled kappa?\n"
        "### If it were, a 512x increase in sub-sample density would move it.\n"
    )
    which = ["plain", "exact, no J", "exact + Jexact"]
    for radius, outside in ((4.0, True), (6.0, False)):
        side = "convex" if outside else "concave"
        print(f"  {side} R={radius}")
        for nsub in (8, 16, 32, 64):
            cfg = build(radius, outside, 2_000_000, 17, nsub)
            fields = fields_for(cfg, which)
            total = cfg["pos"].shape[0] * cfg["mass"]
            vals = {k: float(np.sum(f * KAPPA) * DX**3) / total for k, f in fields.items()}
            print(f"    nsub {nsub:3d}   " + "   ".join(f"{k} {v:.5f}" for k, v in vals.items()))
        # and the volume the kappa field actually represents
        axes = [np.arange(SHAPE[q]) + LO for q in range(3)]
        ii, jj, kk = np.meshgrid(*axes, indexing="ij")
        inbox = np.ones(SHAPE, dtype=bool)
        for a in (ii, jj, kk):
            inbox &= (a >= 0) & (a <= N_VALID - 1)
        v_kappa = float(np.sum(KAPPA[inbox]) * DX**3)
        v_mc = cfg["pos"].shape[0] * cfg["mass"]
        print(f"    sum(kappa)dV / MC fluid volume = {v_kappa / v_mc:.6f}\n")


def main_identity():
    print(
        "### mass ratio == mean_p [ kappatilde(x_p) + J_p kappatilde(R(x_p)) ] ?\n"
        "### If it is an identity, the mass check is not independent evidence.\n"
    )
    for radius, outside in ((200.0, True), (8.0, True), (6.0, False)):
        cfg = build(radius, outside, 2_000_000, 17, 16)
        fields = fields_for(cfg, ["exact, no J", "exact + Jexact"])
        total = cfg["pos"].shape[0] * cfg["mass"]
        side = "convex" if outside else "concave"
        print(f"  {side} R={radius}")
        for label, jac in (("exact, no J", None), ("exact + Jexact", "exact")):
            sel, img, dist = cfg["sel_exact"], cfg["img_exact"], cfg["d_exact"]
            own = interpolated_kappa(cfg["pos"])
            image = np.zeros(cfg["pos"].shape[0])
            factor = 1.0 if jac is None else cfg["jac_exact"](dist[sel])
            image[sel] = interpolated_kappa(img[sel]) * factor
            predicted = float(np.mean(own + image))
            measured = float(np.sum(fields[label] * KAPPA) * DX**3) / total
            print(
                f"    {label:>15}  predicted {predicted:.6f}   measured {measured:.6f}"
                f"   diff {predicted - measured:+.2e}"
            )
        print()


def main_margin():
    print(
        "### Scoring-margin sensitivity on a curved surface.\n"
        "### mirror_sphere.py runs its concave case at margin 1 and the others at 3.\n"
    )
    for radius, outside in ((8.0, True), (6.0, False)):
        for margin in (1, 2, 3, 4):
            run(radius, outside, npart=2_000_000, margin=margin, which=["exact + Jexact"])


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "radii"
    {"radii": main_radii, "nsub": main_nsub, "identity": main_identity, "margin": main_margin}[
        what
    ]()
