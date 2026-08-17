#!/usr/bin/env python3
"""
Extended planar harness: the reflect criterion, the band width, the scoring
margin, and a demonstration that mirror_test.py's `lattice` mode is unsound.

mirror_test.py measures one reflect rule (reflect everything with d > 0) and
asserts a band width. This harness measures the alternatives and the band.

What it establishes:

1. `criteria` -- restricting WHICH particles reflect is the wrong lever.
   "Reflect only if my own stencil touches a covered cell" reflects 1.2% of
   particles (0.00% for axis-aligned normals, which is where a smaller figure
   comes from) and reproduces plain CIC, because a sliver cut cell shields the
   covered cells behind it. The mirror is a boundary condition on the
   reconstruction, not a refund of leaked weight.

2. `band` -- the reflect band is (3/2)*dx*sum|n_i| = 3*s_max, and that bound is
   tight. The image's CIC support is a cube of side 2 dx whose reach ALONG THE
   NORMAL is dx*sum|n_i| = 2*s_max, and the cut-cell centre it must reach can
   sit a further s_max inside the solid. A band of dx*(1 + sum|n_i|/2) counts
   the support reach as dx instead and is short by up to 0.71 dx on tilted
   normals. The two agree exactly at 1.5 dx for an axis-aligned normal, which is
   why the narrow form survives axis-aligned validation.

3. `margin` -- margin 3 is a defensible conservative choice, not a measured
   floor. The margin-1 artifact is real for the two `s max` cases; there is no
   plateau at 3, and `near lo edge` scores ZERO cut cells at margin >= 3, so one
   of the two near-domain-edge cases tests nothing.

4. `lattice` -- mirror_test.py's lattice sampling is documented as the
   deterministic, shot-noise-free mode. At its default sub=4 it reports the
   mirror as +20.9% biased at kappa ~ 0.63 with a standard error of exactly
   0.00000. It is quantization of the fluid boundary by the sampling lattice and
   it oscillates with lattice phase. This is a third instance of the failure
   class the README documents two of, and the zero standard error makes it the
   most convincing of the three. Do not use that mode.

Usage:
    python3 mirror_planar_ext.py criteria [nparticles]
    python3 mirror_planar_ext.py band     [nparticles]
    python3 mirror_planar_ext.py margin   [nparticles]
    python3 mirror_planar_ext.py lattice
"""

import sys

import numpy as np

from mirror_test import (
    CASES,
    DX,
    KAPPA_BINS,
    N_VALID,
    PAD,
    box_fluid_volume,
    cic_stencil,
    deposit,
    kappa_field,
    particle_cell,
    sample_lattice,
    sample_random,
    scoreable_mask,
)

SHAPE = (N_VALID + 2 * PAD,) * 3
LO = -PAD
SIDE = N_VALID * DX


def stencil_kappa(kappa, idx):
    sh = idx - LO
    ins = np.all((sh >= 0) & (sh < np.array(SHAPE)), axis=-1)
    s = np.clip(sh, 0, np.array(SHAPE) - 1)
    kap = np.zeros(idx.shape[:-1])
    kap[ins] = kappa[s[..., 0][ins], s[..., 1][ins], s[..., 2][ins]]
    return kap, ins, s


def deposit_sets(kappa, sets):
    field = np.zeros(SHAPE)
    for pos, mass in sets:
        if pos.shape[0] == 0:
            continue
        idx, wt = cic_stencil(pos)
        kap, ins, s = stencil_kappa(kappa, idx)
        live = (wt > 0) & ins & (kap > 0)
        flat = np.ravel_multi_index((s[..., 0], s[..., 1], s[..., 2]), SHAPE, mode="clip")
        np.add.at(
            field.reshape(-1), flat.ravel(), np.where(live, wt * mass[:, None], 0.0).ravel()
        )
    return field / DX**3


def bin_means(field, kappa, mask):
    out = {}
    for lo, hi in KAPPA_BINS:
        sel = mask & (kappa > lo) & (kappa <= hi)
        n = int(sel.sum())
        out[(lo, hi)] = (n, float(np.mean(field[sel] - 1.0)) if n else float("nan"))
    return out


def setup(case, n_random, rng):
    name, raw, p0 = case
    nhat = raw / np.linalg.norm(raw)
    kappa = kappa_field(nhat, p0, SHAPE, LO)
    pos = sample_random(nhat, p0, n_random, SIDE, rng)
    mass = box_fluid_volume(nhat, p0, SIDE) / pos.shape[0]
    return name, nhat, p0, kappa, pos, mass


def all_criteria(pos, nhat, p0, kappa):
    """Every reflect rule considered, including the two band formulas."""
    d = (pos - p0) @ nhat
    img = pos - 2.0 * d[:, None] * nhat[None, :]
    own_kap, own_ins, _ = stencil_kappa(kappa, cic_stencil(pos)[0])
    img_kap, img_ins, _ = stencil_kappa(kappa, cic_stencil(img)[0])
    sum_n = float(np.abs(nhat).sum())
    return img, d, {
        "reflect all d > 0": d > 0,
        "own stencil hits covered": (d > 0) & np.any(own_ins & (own_kap == 0.0), axis=1),
        "per-coord 2d|n_i| <= dx": (d > 0)
        & np.all(2.0 * d[:, None] * np.abs(nhat)[None, :] <= DX, axis=1),
        "d <= dx": (d > 0) & (d <= DX),
        "d <= dx(1+sum|n_i|/2)  OLD": (d > 0) & (d <= DX * (1.0 + sum_n / 2.0)),
        "d <= (3/2)dx sum|n_i|  NEW": (d > 0) & (d <= 1.5 * DX * sum_n),
        "image support has kappa>0": (d > 0) & np.any(img_ins & (img_kap > 0.0), axis=1),
    }


def main_criteria(n_random):
    print("### Reflect-criterion table (source-patch deposition)\n")
    rng = np.random.default_rng(20260817)
    agg = {}
    for case in CASES:
        name, nhat, p0, kappa, pos, mass = setup(case, n_random, rng)
        mask = scoreable_mask(kappa, LO, SHAPE, SIDE)
        img, _, crit = all_criteria(pos, nhat, p0, kappa)
        base = (pos, np.full(pos.shape[0], mass))
        print(f"--- {name}   n = ({nhat[0]:.3f}, {nhat[1]:.3f}, {nhat[2]:.3f})")
        for label, sel in crit.items():
            field = deposit_sets(kappa, [base, (img[sel], np.full(int(sel.sum()), mass))])
            b = bin_means(field, kappa, mask)
            low = next((v for (lo, hi), (n, v) in b.items() if n and hi <= 0.05), np.nan)
            mid = next((v for (lo, hi), (n, v) in b.items() if n and lo >= 0.05 and hi <= 0.30), np.nan)
            print(
                f"      {label:<28} reflected {100 * sel.mean():6.2f}%"
                f"   kappa<=0.05 {low:+8.4f}   0.05-0.30 {mid:+8.4f}"
            )
            agg.setdefault(label, []).append((sel.mean(), low, mid))
        print()
    print("### aggregate over all cases")
    for label, rows in agg.items():
        arr = np.array(rows, dtype=float)
        print(
            f"    {label:<28} reflected {100 * np.nanmean(arr[:, 0]):6.2f}%"
            f"   kappa<=0.05 {np.nanmean(arr[:, 1]):+8.4f}"
            f"   0.05-0.30 {np.nanmean(arr[:, 2]):+8.4f}"
        )


def main_band(n_random):
    print(
        "### How far out does the exact criterion actually fire?\n"
        "### The band must be a superset of it, or it silently overrides it.\n"
    )
    rng = np.random.default_rng(7)
    print(
        f"    {'case':<18} {'sum|n_i|':>9} {'old band':>9} {'NEW band':>9} "
        f"{'measured max d':>15}  verdict"
    )
    worst = -float("inf")  # NOT 0.0 -- that would clamp away the measured slack
    for case in CASES:
        name, nhat, p0, kappa, pos, _ = setup(case, n_random, rng)
        d = (pos - p0) @ nhat
        img = pos - 2.0 * d[:, None] * nhat[None, :]
        kap, ins, _ = stencil_kappa(kappa, cic_stencil(img)[0])
        fires = (d > 0) & np.any(ins & (kap > 0.0), axis=1)
        sum_n = float(np.abs(nhat).sum())
        old, new = DX * (1.0 + sum_n / 2.0), 1.5 * DX * sum_n
        measured = float(d[fires].max())
        worst = max(worst, measured - new)
        note = "OK" if measured <= old + 1e-9 else f"old band short by {measured - old:.3f} dx"
        print(
            f"    {name:<18} {sum_n:9.3f} {old:9.3f} {new:9.3f} {measured:15.3f}  {note}"
        )
    print(
        f"\n    New band exceeded by at most {worst:+.3f} dx across all cases"
        f" -- {'PASS, and tight' if worst <= 0 else 'FAIL: the new band is short too'}"
    )


def main_margin(n_random):
    print(
        "### Scoring-margin sweep. If margin 3 is a floor, 4 and 5 must agree with it.\n"
    )
    rng = np.random.default_rng(20260817)
    for case in CASES:
        name, nhat, p0, kappa, pos, mass = setup(case, n_random, rng)
        d = (pos - p0) @ nhat
        img = pos - 2.0 * d[:, None] * nhat[None, :]
        sel = d > 0
        field = deposit_sets(
            kappa,
            [(pos, np.full(pos.shape[0], mass)), (img[sel], np.full(int(sel.sum()), mass))],
        )
        print(f"--- {name}")
        for margin in (1, 2, 3, 4, 5):
            mask = scoreable_mask(kappa, LO, SHAPE, SIDE, margin=margin)
            b = bin_means(field, kappa, mask)
            cut = [(n, v) for (lo, hi), (n, v) in b.items() if n and hi <= 0.999]
            total = sum(n for n, _ in cut)
            worst = max((abs(v) for _, v in cut), default=float("nan"))
            low = next(
                (v for (lo, hi), (n, v) in b.items() if n and hi <= 0.05), float("nan")
            )
            print(
                f"      margin {margin}: cut cells scored {total:5d}"
                f"   kappa<=0.05 {low:+9.5f}   worst |dev| {worst:.5f}"
            )
        print()


def main_lattice():
    print(
        "### mirror_test.py's `lattice` mode, swept over sub.\n"
        "### A real bias would not change sign with the lattice spacing.\n"
    )
    name, raw, p0 = CASES[0]
    nhat = raw / np.linalg.norm(raw)
    kappa = kappa_field(nhat, p0, SHAPE, LO)
    mask = scoreable_mask(kappa, LO, SHAPE, SIDE)
    sel = mask & (kappa > 0.5) & (kappa <= 0.75)
    print(f"    plane at x = {p0[0]}, so the cut cells hold kappa = {kappa[8 - LO, 8 - LO, 8 - LO]:.3f}")
    print(f"\n    {'sub':>4} {'spacing':>9} {'particles':>10} {'mirror dev':>12}")
    for sub in (2, 4, 8, 16):
        pos, h = sample_lattice(nhat, p0, SIDE, sub)
        report = deposit(pos, h**3, kappa, "mirror", nhat, p0, LO, SHAPE)
        print(
            f"    {sub:>4} {h:9.4f} {pos.shape[0]:10d}"
            f" {float(np.mean(report.field[sel] - 1.0)):+12.5f}"
        )
        del pos, report
    print(
        "\n    Sign flips between sub=2 and sub=4 and collapses once the lattice\n"
        "    resolves the cut. It is boundary quantization, not a property of the\n"
        "    mirror. Use `random`; the standard error there is honest."
    )


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "criteria"
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 600_000
    if what == "lattice":
        main_lattice()
    else:
        {"criteria": main_criteria, "band": main_band, "margin": main_margin}[what](count)
