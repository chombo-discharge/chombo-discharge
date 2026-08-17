#!/usr/bin/env python3
"""
Single-patch 3-D prototype for mirrored cut-cell deposition (chombo-discharge issue #29).

Domain
------
16^3 valid cells, 2 ghost cells on every side, dx = 1, single patch, serial.
The EB is a plane given by a centroid p0 and a unit normal n. Fluid is the side
with (x - p0).n > 0. Volume fractions kappa are computed analytically (exact
plane-cube intersection), not sampled.

What is compared
----------------
  ngp     NGP-in-cut-cells (today's ItoKMC default). Predicts phi/n = kappa.
  plain   Plain CIC. Weight landing in covered cells is lost.
  mirror  CIC plus one reflected image per particle, covered cells skipped,
          no per-particle normalization.

Instrumentation
---------------
  * uniformity  : phi/n per cell against the exact answer 1, binned by kappa.
                  Only cells whose full CIC support lies inside the sampled
                  region are scored, so the domain edge does not contaminate
                  the EB measurement.
  * containment : how far, in cells, a deposit reaches from the particle's own
                  cell, and whether any write with kappa > 0 lands outside the
                  nominal 2-ghost region. That is the DEBUG assert from the plan.
  * mass        : sum(kappa_j phi_j dx^3) against sum(m_p), the conservation
                  statement in the true-density convention.

Note on scope: for a planar EB every cut cell shares the same plane, so the
per-cell "nearest cut-cell normal" is exact here. This test therefore exercises
the CIC/kappa/ghost mechanics and unbiasedness, and does NOT exercise the
cell-constant-plane approximation. A curved geometry is needed for that.
"""

import math
from dataclasses import dataclass
from itertools import combinations

import numpy as np

N_VALID = 16          # valid cells per direction
N_GHOST = 2           # nominal ghost cells per side
DX = 1.0
PAD = N_GHOST + 4     # internal pad, so out-of-range writes are observed not crashed


# --------------------------------------------------------------------------- #
# Geometry: exact volume fraction of a cube cut by a plane
# --------------------------------------------------------------------------- #
def unit_cube_solid_fraction(nhat, d_local):
    """
    Fraction of the unit cube [0,1]^3 satisfying nhat . y <= d_local.

    Inclusion-exclusion over the nonzero components of nhat, so axis-aligned
    and edge-aligned normals are handled exactly rather than by perturbation.
    d_local may be an array.
    """
    d_local = np.asarray(d_local, dtype=float)
    a = np.abs(nhat)
    # y_i -> 1 - y_i for negative components maps the problem to all-positive a
    d = d_local + np.sum(np.where(np.asarray(nhat) < 0.0, a, 0.0))

    idx = np.flatnonzero(a > 1e-12)
    k = idx.size
    if k == 0:
        return np.where(d >= 0.0, 1.0, 0.0)

    aa = a[idx]
    acc = np.zeros_like(d)
    for r in range(k + 1):
        for S in combinations(range(k), r):
            shifted = d - sum(aa[i] for i in S)
            acc += ((-1) ** r) * np.maximum(shifted, 0.0) ** k

    return np.clip(acc / (math.factorial(k) * np.prod(aa)), 0.0, 1.0)


def kappa_field(nhat, p0, shape, lo_index):
    """Exact fluid volume fraction for every cell of the padded array."""
    ii, jj, kk = np.meshgrid(
        np.arange(shape[0]) + lo_index,
        np.arange(shape[1]) + lo_index,
        np.arange(shape[2]) + lo_index,
        indexing="ij",
    )
    corner = np.stack([ii, jj, kk], axis=-1).astype(float) * DX   # lower corner
    nhat = np.asarray(nhat, dtype=float)

    # Decide entirely-solid / entirely-fluid cells from the corner extrema first.
    # The inclusion-exclusion below sums terms of order d^3 to produce an answer
    # of order 1, so for cells far from the plane it is pure cancellation noise
    # and returns kappa ~ 1e-13 where the answer is exactly 0. A real EB hands
    # back exact covered/regular flags; this reproduces that.
    base = (corner - p0) @ nhat
    lo_val = base + DX * np.sum(np.minimum(nhat, 0.0))
    hi_val = base + DX * np.sum(np.maximum(nhat, 0.0))

    d_local = ((p0 - corner) @ nhat) / DX
    kap = 1.0 - unit_cube_solid_fraction(nhat, d_local)
    kap = np.where(hi_val <= 0.0, 0.0, kap)     # every corner in the solid
    kap = np.where(lo_val >= 0.0, 1.0, kap)     # every corner in the fluid
    return kap


def box_fluid_volume(nhat, p0, side):
    """Fluid volume of the cube [0,side]^3 cut by the plane."""
    d_local = float(np.dot(p0, nhat)) / side
    return (side ** 3) * (1.0 - float(unit_cube_solid_fraction(nhat, d_local)))


# --------------------------------------------------------------------------- #
# Particles
# --------------------------------------------------------------------------- #
def sample_random(nhat, p0, n_target, side, rng):
    """Uniform random points in the fluid part of [0,side]^3."""
    out = []
    have = 0
    while have < n_target:
        batch = rng.uniform(0.0, side, size=(max(n_target, 1024), 3))
        keep = batch[(batch - p0) @ nhat > 0.0]
        out.append(keep)
        have += keep.shape[0]
    return np.concatenate(out, axis=0)[:n_target]


def sample_lattice(nhat, p0, side, sub):
    """
    Regular lattice, `sub` points per cell per direction, offset by half a
    spacing. Deterministic, so the bias measurement is not buried in shot noise.
    """
    h = DX / sub
    axis = (np.arange(int(side / h)) + 0.5) * h
    gx, gy, gz = np.meshgrid(axis, axis, axis, indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=-1)
    return pts[(pts - p0) @ nhat > 0.0], h


# --------------------------------------------------------------------------- #
# Deposition
# --------------------------------------------------------------------------- #
def cic_stencil(pos):
    """
    Returns (idx, wt) with idx of shape (Np, 8, 3) and wt of shape (Np, 8).
    idx values are absolute cell indices (may be negative).
    """
    g = pos / DX - 0.5
    base = np.floor(g).astype(np.int64)
    frac = g - base

    idx = np.empty((pos.shape[0], 8, 3), dtype=np.int64)
    wt = np.ones((pos.shape[0], 8), dtype=float)
    for c, (ox, oy, oz) in enumerate(np.ndindex(2, 2, 2)):
        for dim, off in enumerate((ox, oy, oz)):
            idx[:, c, dim] = base[:, dim] + off
            wt[:, c] *= frac[:, dim] if off else (1.0 - frac[:, dim])
    return idx, wt


def particle_cell(pos):
    return np.floor(pos / DX).astype(np.int64)


@dataclass
class Report:
    field: np.ndarray
    max_reach: int              # max |index offset| particle cell -> written cell
    max_ghost_layer: int        # deepest ghost layer receiving a kappa>0 write
    escaped_fluid: int          # kappa>0 writes outside the nominal ghost region
    escaped_covered: int        # kappa==0 write attempts outside it (harmless)
    deposited_mass: float
    weight_by_reach: dict       # reach in cells -> deposited weight fraction
    escaped_weight: float       # weight written outside the nominal ghost region


def deposit(pos, m, kappa, scheme, nhat, p0, lo_index, shape):
    """
    scheme in {"ngp", "plain", "mirror"}.
    Covered cells are skipped for "plain" and "mirror" alike -- writing into them
    and ignoring the result is the same thing, and skipping is what the mirror
    requires. "ngp" forces the whole weight into an irregular cell.
    """
    field = np.zeros(shape, dtype=float)
    pcell = particle_cell(pos)

    def accumulate(idx, wt, track, ref=None):
        nonlocal max_reach, max_ghost, esc_fluid, esc_cov, esc_wt
        ref = pcell if ref is None else ref
        shifted = idx - lo_index
        inside = np.all((shifted >= 0) & (shifted < np.array(shape)), axis=-1)
        live = (wt > 0.0) & inside

        kap = np.zeros(wt.shape, dtype=float)
        s = np.clip(shifted, 0, np.array(shape) - 1)
        kap[live] = kappa[s[..., 0][live], s[..., 1][live], s[..., 2][live]]

        write = live & (kap > 0.0)

        # containment instrumentation
        if track and np.any(write):
            reach_all = np.max(np.abs(idx - ref[:, None, :]), axis=-1)
            for r in np.unique(reach_all[write]):
                sel = write & (reach_all == r)
                by_reach[int(r)] = by_reach.get(int(r), 0.0) + float(
                    np.sum(wt[sel]) * m)
            out_all = np.any((idx < -N_GHOST) | (idx >= N_VALID + N_GHOST),
                             axis=-1)
            esc_wt += float(np.sum(wt[write & out_all]) * m)
            reach = np.max(reach_all[write])
            max_reach = max(max_reach, int(reach))
            wi = idx[write]
            low = int(np.min(wi)) if wi.size else 0
            high = int(np.max(wi)) if wi.size else 0
            max_ghost = max(max_ghost, max(0, -low), max(0, high - (N_VALID - 1)))
            out = (wi < -N_GHOST) | (wi >= N_VALID + N_GHOST)
            esc_fluid += int(np.count_nonzero(np.any(out, axis=-1)))
        if track:
            attempted_out = (wt > 0.0) & np.any(
                (idx < -N_GHOST) | (idx >= N_VALID + N_GHOST), axis=-1
            )
            esc_cov += int(np.count_nonzero(attempted_out & ~write))

        flat = np.ravel_multi_index(
            (s[..., 0], s[..., 1], s[..., 2]), shape, mode="clip"
        )
        contrib = np.where(write, wt * m, 0.0)
        np.add.at(field.reshape(-1), flat.ravel(), contrib.ravel())

    max_reach, max_ghost, esc_fluid, esc_cov = 0, 0, 0, 0
    by_reach, esc_wt = {}, 0.0

    if scheme == "ngp":
        s = pcell - lo_index
        kap = kappa[s[:, 0], s[:, 1], s[:, 2]]
        irreg = (kap > 0.0) & (kap < 1.0)
        idx, wt = cic_stencil(pos)
        # regular cells use CIC; irregular ones take the whole weight (NGP)
        wt = np.where(irreg[:, None], 0.0, wt)
        accumulate(idx, wt, track=True)
        ngp_idx = pcell[irreg][:, None, :]
        ngp_wt = np.ones((int(np.count_nonzero(irreg)), 1))
        if ngp_idx.size:
            accumulate(ngp_idx, ngp_wt, track=True, ref=pcell[irreg])
    else:
        idx, wt = cic_stencil(pos)
        accumulate(idx, wt, track=True)

        if scheme == "mirror":
            d = (pos - p0) @ nhat                       # > 0 in the fluid
            image = pos - 2.0 * d[:, None] * nhat[None, :]
            iidx, iwt = cic_stencil(image)
            iwt = np.where((d > 0.0)[:, None], iwt, 0.0)
            # "reflect iff the image's support holds a kappa>0 centre" is
            # enforced implicitly: cells with kappa==0 are skipped, so an image
            # buried in the solid contributes nothing.
            accumulate(iidx, iwt, track=True)

    field /= DX ** 3
    mass = float(np.sum(field * kappa) * DX ** 3)
    total = float(np.sum(pos.shape[0] * m))
    return Report(field, max_reach, max_ghost, esc_fluid, esc_cov,
                  mass, {k: v / total for k, v in sorted(by_reach.items())},
                  esc_wt / total)


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #
def scoreable_mask(kappa, lo_index, shape, side, margin=3):
    """
    Cells with kappa > 0 that are far enough from the edge of the sampled region
    for phi to be fully supplied.

    margin = 3, not 1. A cell centre on the solid side is reconstructed almost
    entirely from images, i.e. from particles near R(c), which sits up to
    sqrt(3)*dx from c. Requiring only c's own CIC support to be inside the box
    leaves such cells under-supplied and fakes a ~1.7% deficit at small kappa
    that is an artifact of the harness, not of the scheme. The requirement is
    1 + sqrt(3) = 2.73 -> 3.
    """
    ii, jj, kk = np.meshgrid(
        np.arange(shape[0]) + lo_index,
        np.arange(shape[1]) + lo_index,
        np.arange(shape[2]) + lo_index,
        indexing="ij",
    )
    interior = np.ones(shape, dtype=bool)
    for a in (ii, jj, kk):
        interior &= (a >= margin) & (a <= int(side / DX) - 1 - margin)
    return interior & (kappa > 0.0)


KAPPA_BINS = [(0.0, 0.05), (0.05, 0.15), (0.15, 0.30), (0.30, 0.50),
              (0.50, 0.75), (0.75, 0.999), (0.999, 1.001)]


def bin_report(field, kappa, mask):
    rows = []
    for lo, hi in KAPPA_BINS:
        sel = mask & (kappa > lo) & (kappa <= hi)
        n = int(np.count_nonzero(sel))
        if n == 0:
            rows.append((lo, hi, 0, float("nan"), float("nan"), float("nan")))
            continue
        dev = field[sel] - 1.0
        sem = float(np.std(dev, ddof=1) / math.sqrt(n)) if n > 1 else float("nan")
        rows.append((lo, hi, n, float(np.mean(dev)), sem,
                     float(np.max(np.abs(dev)))))
    return rows


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
CASES = [
    ("axis-aligned  x", np.array([1.0, 0.0, 0.0]), np.array([8.37, 8.0, 8.0])),
    ("edge  45 in xy",  np.array([1.0, 1.0, 0.0]), np.array([8.30, 8.0, 8.0])),
    ("body diagonal",   np.array([1.0, 1.0, 1.0]), np.array([8.20, 8.0, 8.0])),
    ("oblique",         np.array([0.2, 0.3, 0.93]), np.array([8.11, 8.0, 8.0])),
    ("near lo edge",    np.array([1.0, 0.0, 0.0]), np.array([2.40, 8.0, 8.0])),
    ("near lo, tilted", np.array([1.0, 0.7, 0.3]), np.array([2.35, 8.0, 8.0])),
    # centre-on-solid-side cases: s > 0, which is what makes the mirror reach
    ("axis, s=0.12",    np.array([1.0, 0.0, 0.0]), np.array([8.62, 8.0, 8.0])),
    ("axis, s=0.45",    np.array([1.0, 0.0, 0.0]), np.array([8.95, 8.0, 8.0])),
    ("45xy, s max",     np.array([1.0, 1.0, 0.0]), np.array([8.99, 8.0, 8.0])),
    ("diag, s max",     np.array([1.0, 1.0, 1.0]), np.array([8.99, 8.0, 8.0])),
]


def run(sampling="random", n_random=50_000, sub=4, seed=20260817):
    side = N_VALID * DX
    lo_index = -PAD
    shape = (N_VALID + 2 * PAD,) * 3
    rng = np.random.default_rng(seed)

    print(f"grid {N_VALID}^3, nominal ghost {N_GHOST}, dx {DX}, "
          f"sampling={sampling}\n")

    worst_reach = 0
    any_escape = 0

    for name, raw_n, p0 in CASES:
        nhat = raw_n / np.linalg.norm(raw_n)
        kappa = kappa_field(nhat, p0, shape, lo_index)

        if sampling == "lattice":
            pos, h = sample_lattice(nhat, p0, side, sub)
            m = h ** 3                        # unit density exactly
        else:
            pos = sample_random(nhat, p0, n_random, side, rng)
            m = box_fluid_volume(nhat, p0, side) / pos.shape[0]

        mask = scoreable_mask(kappa, lo_index, shape, side)
        expected = float(np.sum(pos.shape[0] * m))

        print(f"=== {name}   n = ({nhat[0]:.3f}, {nhat[1]:.3f}, {nhat[2]:.3f})"
              f"   particles = {pos.shape[0]}   scored cells = "
              f"{int(np.count_nonzero(mask))}")

        reports = {}
        for scheme in ("ngp", "plain", "mirror"):
            reports[scheme] = deposit(pos, m, kappa, scheme, nhat, p0,
                                      lo_index, shape)

        print(f"    {'kappa bin':>14} {'cells':>6} "
              f"{'ngp mean':>10} {'plain mean':>11} "
              f"{'mirror mean +- stderr':>24}")
        rows = {s: bin_report(reports[s].field, kappa, mask) for s in reports}
        for r in range(len(KAPPA_BINS)):
            lo, hi, n = rows["ngp"][r][:3]
            if n == 0:
                continue
            print(f"    {lo:5.3f}-{hi:5.3f} {n:6d} "
                  f"{rows['ngp'][r][3]:+10.4f} {rows['plain'][r][3]:+11.4f} "
                  f"{rows['mirror'][r][3]:+14.5f} +-{rows['mirror'][r][4]:8.5f}")

        rep = reports["mirror"]
        worst_reach = max(worst_reach, rep.max_reach)
        any_escape += rep.escaped_fluid
        tail = "  ".join(f"{k}:{v:.3e}" for k, v in rep.weight_by_reach.items())
        print(f"    weight by reach (fraction of total): {tail}")
        print(f"    weight written outside nominal ghosts: "
              f"{rep.escaped_weight:.3e}")
        print(f"    containment: reach {rep.max_reach} cells, "
              f"deepest ghost layer {rep.max_ghost_layer}/{N_GHOST}, "
              f"kappa>0 writes outside nominal ghosts: {rep.escaped_fluid}, "
              f"covered attempts outside: {rep.escaped_covered}")
        print(f"    mass  sum(kappa phi dV) / sum(m):  "
              f"ngp {reports['ngp'].deposited_mass/expected:.5f}   "
              f"plain {reports['plain'].deposited_mass/expected:.5f}   "
              f"mirror {rep.deposited_mass/expected:.5f}\n")

    print(f"OVERALL  worst mirror reach = {worst_reach} cells "
          f"(nominal ghost = {N_GHOST})   escapes = {any_escape}")
    print("PASS" if (worst_reach <= N_GHOST and any_escape == 0)
          else "FAIL: mirror writes outside the nominal ghost region")


if __name__ == "__main__":
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else "random"
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 50_000
    run(sampling=mode, n_random=count)
