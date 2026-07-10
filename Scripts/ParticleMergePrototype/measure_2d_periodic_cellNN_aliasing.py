"""
Third reference point in the aliasing ladder: "NN merging per cell" -- restrict the
nearest-neighbor SEARCH to within a single cell (never look at neighbors outside it, no cross-
patch coordination needed at all -- a cell belongs entirely to one patch), but position the merged
particle at the TRUE weighted centroid of the pair, not snapped to the cell center. This isolates
whether restricting search SCOPE alone (independent of the earlier, much cruder cell-center-snap
pathology) produces patch- or cell-frequency aliasing -- i.e. does the mere fact of never looking
outside a cell (as opposed to our real algorithm's full patch+ghost search) impose any spatial
structure, even when the position rule itself is "honest" (true centroid, real data)?

Same periodic domain/injection/generation-loop setup as the other 2D periodic measurements, for
direct comparability. No MPI coordination is needed for the merge step itself (every merge is
local to one cell, hence one patch) -- only for injection/gather bookkeeping, mirroring the
existing driver structure.
"""

import argparse
import json
import random
import sys

import numpy as np
from mpi4py import MPI

import merge_lib_2d_periodic as ml
from run_scenario_2d_periodic import patch_to_rank


def radial_power_spectrum(hist, bin_width, n_bins_radial=40):
    n = hist.shape[0]
    f = np.fft.fft2(hist - hist.mean())
    power2d = np.abs(f) ** 2
    freqs = np.fft.fftfreq(n, d=bin_width)
    kx, ky = np.meshgrid(freqs, freqs, indexing="ij")
    kmag = np.sqrt(kx ** 2 + ky ** 2)
    k_max = freqs.max()
    k_edges = np.linspace(0, k_max, n_bins_radial + 1)
    k_centers = 0.5 * (k_edges[:-1] + k_edges[1:])
    power_radial = np.zeros(n_bins_radial)
    flat_k = kmag.ravel()
    flat_p = power2d.ravel()
    idx = np.digitize(flat_k, k_edges) - 1
    for b in range(n_bins_radial):
        mask = idx == b
        if mask.any():
            power_radial[b] = flat_p[mask].mean()
    return k_centers, power_radial


def cell_nn_merge_round(local, my_patches, cell_dx, Lx, Ly, thresh, id_counter):
    """Within each cell, independently: repeatedly find the nearest-neighbor pair and combine
    them at their TRUE weighted centroid, until the cell's count drops to threshold. Purely
    local -- no ghosts, no cross-patch/cross-rank interaction, matching a cell's guarantee of
    belonging entirely to one patch."""
    n_merged = 0
    for pid in my_patches:
        cells = {}
        for p in local[pid]:
            k = ml.cell_key(p["pos"], cell_dx, Lx, Ly)
            cells.setdefault(k, []).append(p)
        new_list = []
        for cell_key, plist in cells.items():
            while len(plist) > thresh:
                best_i, best_j, best_d2 = None, None, float("inf")
                for a in range(len(plist)):
                    for b in range(a + 1, len(plist)):
                        d2 = ml.dist2(plist[a]["pos"], plist[b]["pos"])
                        if d2 < best_d2:
                            best_d2, best_i, best_j = d2, a, b
                p1, p2 = plist[best_i], plist[best_j]
                for idx in sorted((best_i, best_j), reverse=True):
                    plist.pop(idx)
                id_counter[0] += 1
                merged = ml.combine(p1, p2, id_counter[0], p1["rank"], pid)
                plist.append(merged)
                n_merged += 1
            new_list.extend(plist)
        local[pid] = new_list
    return local, n_merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-patches-x", type=int, default=5)
    ap.add_argument("--n-patches-y", type=int, default=5)
    ap.add_argument("--patch-size", type=float, default=1.0)
    ap.add_argument("--cell-dx", type=float, default=0.2)
    ap.add_argument("--thresh", type=int, default=2)
    ap.add_argument("--n-inject-per-patch", type=int, default=40)
    ap.add_argument("--n-generations", type=int, default=25)
    ap.add_argument("--hist-bins-per-cell", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    n_patches = args.n_patches_x * args.n_patches_y
    if size > args.n_patches_y:
        if rank == 0:
            print(f"size ({size}) must not exceed n_patches_y ({args.n_patches_y})", file=sys.stderr)
        comm.Abort(1)

    my_patches = [pid for pid in range(n_patches)
                  if patch_to_rank(pid, args.n_patches_x, args.n_patches_y, size) == rank]
    Lx = args.patch_size * args.n_patches_x
    Ly = args.patch_size * args.n_patches_y

    gen_rng = random.Random(args.seed * 7919 + rank)
    id_counter = [100_000_000 + rank * 1_000_000]

    local = {pid: [] for pid in my_patches}
    n_total_trace = []
    for gen in range(args.n_generations):
        for pid in my_patches:
            i, j = pid % args.n_patches_x, pid // args.n_patches_x
            new = ml.gen_particles_uniform(pid, i, j, args.patch_size, gen_rng, args.n_inject_per_patch, id_counter)
            local[pid].extend(new)

        local, n_merged = cell_nn_merge_round(local, my_patches, args.cell_dx, Lx, Ly, args.thresh, id_counter)

        n_local = sum(len(local[pid]) for pid in my_patches)
        n_global = comm.allreduce(n_local, op=MPI.SUM)
        if rank == 0:
            n_total_trace.append(n_global)

    my_final = [p for pid in local for p in local[pid]]
    all_final = comm.gather(my_final, root=0)

    if rank == 0:
        final_flat = [p for sub in all_final for p in sub]
        n_final = len(final_flat)
        positions = [p["pos"] for p in final_flat]

        null_rng = random.Random(args.seed + 424242)
        null_positions = [(null_rng.uniform(0, Lx), null_rng.uniform(0, Ly)) for _ in range(n_final)]

        n_hist_bins = int(round(Lx / args.cell_dx * args.hist_bins_per_cell))
        bin_width = Lx / n_hist_bins

        hist_final, _, _ = np.histogram2d([p[0] for p in positions], [p[1] for p in positions],
                                           bins=n_hist_bins, range=[[0, Lx], [0, Ly]])
        hist_null, _, _ = np.histogram2d([p[0] for p in null_positions], [p[1] for p in null_positions],
                                          bins=n_hist_bins, range=[[0, Lx], [0, Ly]])

        k_final, p_final = radial_power_spectrum(hist_final, bin_width)
        k_null, p_null = radial_power_spectrum(hist_null, bin_width)

        report = {
            "n_patches_x": args.n_patches_x, "n_patches_y": args.n_patches_y,
            "n_generations": args.n_generations, "n_final": n_final, "n_total_trace": n_total_trace,
            "k_patch": 1.0 / args.patch_size, "k_cell": 1.0 / args.cell_dx,
            "k_values": k_final.tolist(), "power_final": p_final.tolist(), "power_null": p_null.tolist(),
        }
        with open(args.out, "w") as f:
            json.dump(report, f, indent=2)

        print(f"n_final={n_final} n_total_trace(last 8)={n_total_trace[-8:]}")


if __name__ == "__main__":
    main()
