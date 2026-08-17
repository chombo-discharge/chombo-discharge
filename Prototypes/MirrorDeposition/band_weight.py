"""
How much of the mirrored mass actually lands in kappa>0 cells, as a function of the
source particle's distance d to the EB?

The image deposits into its 8 CIC cells; only those with kappa>0 are written, the rest
is discarded. So the useful weight of one image is sum_{kappa_j>0} W(x_j - R(p)).
This prices any proposed cutoff "do not reflect beyond d_cut".
"""
import numpy as np
from mirror_test import CASES, DX, N_VALID, PAD, kappa_field, cic_stencil

rng = np.random.default_rng(20260817)
lo_index = -PAD
shape = (N_VALID + 2 * PAD,) * 3
NP = 2_000_000
CUTS = [1.0, 1.5, 2.0, 2.5]

print(f"{'case':18s} {'3*s_max':>8} {'max d':>7} {'mean w':>7}   "
      + "  ".join(f"{'>':>1}{c:<4.1f}dx" for c in CUTS))
print(f"{'':18s} {'(dx)':>8} {'(dx)':>7} {'/img':>7}   "
      + "  ".join(f"{'(% of mirrored mass)':>7}" if i == 0 else f"{'':>7}"
                 for i, _ in enumerate(CUTS)))

for name, raw, p0 in CASES:
    nhat = np.asarray(raw, float)
    nhat = nhat / np.linalg.norm(nhat)
    kap = kappa_field(nhat, p0, shape, lo_index)
    s_max = 0.5 * DX * np.abs(nhat).sum()

    pos = rng.uniform(lo_index * DX, (lo_index + shape[0]) * DX, size=(NP, 3))
    d = (pos - p0) @ nhat
    keep = d > 0.0
    pos, d = pos[keep], d[keep]

    img = pos - 2.0 * d[:, None] * nhat[None, :]
    idx, wt = cic_stencil(img)
    loc = idx - lo_index
    inb = np.all((loc >= 0) & (loc < shape[0]), axis=-1)
    g = np.clip(loc, 0, shape[0] - 1)
    kv = kap[g[..., 0], g[..., 1], g[..., 2]]
    useful = (wt * (kv > 0.0) * inb).sum(axis=1)      # weight landing in fluid cells

    tot = useful.sum()
    row = [f"{100.0 * useful[d > c].sum() / tot:7.3f}" for c in CUTS]
    print(f"{name:18s} {3 * s_max:8.3f} {d[useful > 0].max():7.3f} "
          f"{useful[useful > 0].mean():7.4f}   " + "  ".join(row))
