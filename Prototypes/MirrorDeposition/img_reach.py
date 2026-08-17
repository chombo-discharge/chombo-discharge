"""Where the IMAGE's cell sits, in cells, for both deposition kernels.

Two columns matter, and they answer two different questions in the plan:

  img -> cut   Chebyshev distance from the image's cell to the nearest CUT cell. The plan
               argues from a 1-D construction that this is always <= 1 ("the image is never
               two cells away -- the particle is"). This measures it in 3-D on tilted normals.

  img -> src   Chebyshev distance from the image's cell to the SOURCE particle's cell. Plus
               the kernel reach of 1, this is the ghost width a design would need if it
               skipped the remap and let each source patch deposit its own images.

Uses the real kernels: CIC is a top-hat cloud, TSC is a triangular cloud. Only the support
enters the selection here, but the two are kept distinct so nothing is copied out wrong.

Run:  python3 img_reach.py
"""
import numpy as np

from mirror_test import CASES, DX, N_VALID, PAD, kappa_field, particle_cell

NP = 300_000
SEED = 20260817
SHAPE = (N_VALID + 2 * PAD,) * 3
LO = -PAD
KERNELS = (("CIC", 1.0), ("TSC", 2.0))


def useful_weight(pos, kappa, L):
    """Summed cloud weight landing in kappa>0 cells, per particle."""
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


def cheb(a_cells, a_ref):
    """Min over reference cells of the Chebyshev distance, per row of a_cells."""
    return np.abs(a_cells[:, None, :] - a_ref[None, :, :]).max(axis=2).min(axis=1)


def main():
    for label, L in KERNELS:
        rng = np.random.default_rng(SEED)

        print(f"### {label}  (cloud width L = {L:.0f}*dx)")
        print(f"    {'case':18s} {'img->cut':>9} {'src->cut':>9} {'img->src':>9} {'ghost if no remap':>19}")

        for name, raw, p0 in CASES:
            nhat = np.asarray(raw, float)
            nhat = nhat / np.linalg.norm(nhat)
            kap = kappa_field(nhat, p0, SHAPE, LO)
            cut_idx = np.argwhere((kap > 0.0) & (kap < 1.0)) + LO

            pos = rng.uniform(LO * DX, (LO + SHAPE[0]) * DX, size=(NP, 3))
            d = (pos - p0) @ nhat
            pos = pos[d > 0.0]
            d = d[d > 0.0]

            img = pos - 2.0 * d[:, None] * nhat[None, :]
            sel = useful_weight(img, kap, L) > 0.0

            src_cells = particle_cell(pos[sel])
            img_cells = particle_cell(img[sel])

            i2c = int(cheb(img_cells, cut_idx).max())
            s2c = int(cheb(src_cells, cut_idx).max())
            i2s = int(np.abs(img_cells - src_cells).max(axis=1).max())

            print(f"    {name:18s} {i2c:9d} {s2c:9d} {i2s:9d} {i2s + 1:19d}")

        print()


if __name__ == "__main__":
    main()
