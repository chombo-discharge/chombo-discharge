"""
The reflect band, stated in CELLS rather than in distance, for both deposition kernels.

Uses the ACTUAL chombo-discharge overlap-integral weights (cloud of width L*dx against
the cell volume), not point-sampled hat functions. L = 1 is CIC, L = 2 is TSC.

Three tables per kernel, over the ten planar cases of mirror_test.CASES:

  reach      For every particle the EXACT reflect criterion accepts, the Chebyshev
             distance from the particle's own cell to the nearest CUT cell. This is
             what sets the band mask width, because a particle always sits in a valid
             cell and reads that cell's stored (x_c, n_c, S_c).

  retention  If the band mask is "cells within Chebyshev distance N of a cut cell",
             what fraction of the mirrored mass survives? Mirrored mass is counted as
             it is actually deposited: only the part of an image's cloud landing in a
             kappa>0 cell, since covered cells are never written.

Run:  python3 reach_cells.py
"""
import numpy as np

from mirror_test import CASES, DX, N_VALID, PAD, kappa_field, particle_cell

NP = 300_000
SEED = 20260817
SHAPE = (N_VALID + 2 * PAD,) * 3
LO = -PAD
KERNELS = (("CIC", 1.0), ("TSC", 2.0))


def cloud_weights(pos, kappa, L):
    """Overlap weight of a width-L*dx cloud at pos against every nearby cell.

    The two kernels are NOT the same cloud shape, and mass-weighted results depend on it:
    CIC is a TOP-HAT cloud, weight = min(b,L/2) - max(a,-L/2)      (CD_EBParticleMesh.H:732)
    TSC is a TRIANGULAR cloud, the (be-al) - (be|be|-al|al|)/L form (CD_EBParticleMesh.H:794)
    Only the support is shared. Revisions 3-4 of the plan used the TSC form for both and
    reported CIC retention of 95.5% at Chebyshev 1 where the true figure is 91.6%.

    Returns (useful, ) where useful[p] is the summed weight landing in kappa>0 cells.
    """
    cic = L == 1.0
    base = np.floor(pos / DX).astype(np.int64)
    rad = int(np.ceil(L / 2.0)) + 1
    useful = np.zeros(pos.shape[0])

    for ox in range(-rad, rad + 1):
        for oy in range(-rad, rad + 1):
            for oz in range(-rad, rad + 1):
                iv = base + np.array([ox, oy, oz])
                w = np.ones(pos.shape[0])

                for d in range(3):
                    a = iv[:, d] - pos[:, d] / DX
                    b = a + 1.0
                    al = np.maximum(a, -0.5 * L)
                    be = np.minimum(b, 0.5 * L)
                    ok = al < be
                    ww = (be - al) if cic else ((be - al) - (be * np.abs(be) - al * np.abs(al)) / L)
                    w *= np.where(ok, ww, 0.0)

                sh = iv - LO
                ins = np.all((sh >= 0) & (sh < np.array(SHAPE)), axis=-1)
                s = np.clip(sh, 0, np.array(SHAPE) - 1)
                kap = np.zeros(pos.shape[0])
                kap[ins] = kappa[s[ins, 0], s[ins, 1], s[ins, 2]]
                useful += np.where((w > 1e-14) & ins & (kap > 0.0), w, 0.0)

    return useful


def measure(raw, p0, L, rng):
    nhat = np.asarray(raw, float)
    nhat = nhat / np.linalg.norm(nhat)
    kap = kappa_field(nhat, p0, SHAPE, LO)
    cut_idx = np.argwhere((kap > 0.0) & (kap < 1.0)) + LO

    pos = rng.uniform(LO * DX, (LO + SHAPE[0]) * DX, size=(NP, 3))
    d = (pos - p0) @ nhat
    pos, d = pos[d > 0.0], d[d > 0.0]

    img = pos - 2.0 * d[:, None] * nhat[None, :]
    useful = cloud_weights(img, kap, L)

    sel = useful > 0.0
    cells = particle_cell(pos[sel])
    cheb = np.abs(cells[:, None, :] - cut_idx[None, :, :]).max(axis=2).min(axis=1)
    return d[sel], cheb, useful[sel]


def main():
    for label, L in KERNELS:
        rng = np.random.default_rng(SEED)
        rows = [(name, *measure(raw, p0, L, rng)) for name, raw, p0 in CASES]

        print(f"### {label}  (cloud width L = {L:.0f}*dx)")
        print(f"    {'case':18s} {'max d':>7} {'max cheb':>9} "
              f"{'N=0':>9} {'N=1':>9} {'N=2':>9} {'N=3':>9}   (% mass retained)")
        for name, d, cheb, useful in rows:
            tot = useful.sum()
            keep = [100.0 * useful[cheb <= n].sum() / tot for n in (0, 1, 2, 3)]
            print(f"    {name:18s} {d.max():7.3f} {cheb.max():9d} "
                  + " ".join(f"{k:9.3f}" for k in keep))
        print()


if __name__ == "__main__":
    main()
