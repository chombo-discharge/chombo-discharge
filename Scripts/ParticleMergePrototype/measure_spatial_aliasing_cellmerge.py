"""
Control/reference case for the aliasing comparison: the classic naive "merge by cell" scheme --
whenever a cell holds more than N particles, repeatedly combine pairs within that cell, with the
merged particle's position SNAPPED TO THE CELL CENTER (not the true pairwise centroid). This is
the textbook pathological case the whole aliasing investigation was checking our nearest-neighbor
merger against: many merged particles across the domain land at the exact same relative position
within their cell, which is expected to produce a strong, sharp density comb at the cell spatial
frequency (and, since it never looks across a patch boundary either, a separate comb at the patch
frequency too).

Deliberately structural: a pure cell-based merge never needs to look outside its own cell, so
(unlike our real algorithm) it needs no ghost-fill and no cross-rank coordination at all -- each
patch/rank merges its own cells in complete isolation. Same domain, same injection protocol, same
generation count, same seeds as the no-motion multi-patch measurement of the real algorithm, so
the two results are directly, quantitatively comparable.
"""

import argparse
import json
import random
import sys

import numpy as np
from mpi4py import MPI

import merge_lib as ml


def make_uniform_grid(n_side):
    patches = []
    rank = 0
    for j in range(n_side):
        for i in range(n_side):
            patches.append(((float(i), float(j)), (float(i + 1), float(j + 1)), rank))
            rank += 1
    return patches, rank


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


def cell_merge_round(local, my_patches, cell_dx, thresh, id_counter):
    """Naive cell-based merge: within each cell independently, repeatedly find the NEAREST-NEIGHBOR
    pair of particles within that cell and combine them into one, positioned at the CELL CENTER
    (never the true pairwise centroid), until the cell's count drops to the threshold. No
    cross-patch/cross-rank interaction at all -- purely local per cell. This isolates the position
    rule (cell-center snap vs. true centroid) as the only difference from our real algorithm --
    partner SELECTION is the same nearest-neighbor rule in both cases."""
    n_merged = 0
    for pid in my_patches:
        cells = {}
        for p in local[pid]:
            k = ml.cell_key(p["pos"], cell_dx)
            cells.setdefault(k, []).append(p)
        new_list = []
        for cell_key, plist in cells.items():
            while len(plist) > thresh:
                best_i, best_j, best_d2 = None, None, float("inf")
                for i in range(len(plist)):
                    for j in range(i + 1, len(plist)):
                        d2 = ml.dist2(plist[i]["pos"], plist[j]["pos"])
                        if d2 < best_d2:
                            best_d2, best_i, best_j = d2, i, j
                a, b = plist[best_i], plist[best_j]
                for idx in sorted((best_i, best_j), reverse=True):
                    plist.pop(idx)
                id_counter[0] += 1
                cell_center = tuple((c + 0.5) * cell_dx for c in cell_key)
                merged = {
                    "pos": cell_center,
                    "weight": a["weight"] + b["weight"],
                    "id": id_counter[0],
                    "energy": (a["weight"] * a["energy"] + b["weight"] * b["energy"]) / (a["weight"] + b["weight"]),
                    "rank": a["rank"],
                    "patch": pid,
                    "lineage": a["lineage"] | b["lineage"],
                }
                plist.append(merged)
                n_merged += 1
            new_list.extend(plist)
        local[pid] = new_list
    return local, n_merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-side", type=int, default=5)
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

    patches, nranks = make_uniform_grid(args.n_side)
    if size != nranks:
        if rank == 0:
            print(f"need {nranks} ranks for a {args.n_side}x{args.n_side} grid, got {size}", file=sys.stderr)
        comm.Abort(1)

    my_patches = [pid for pid, box in enumerate(patches) if box[2] == rank]
    local = {pid: [] for pid in my_patches}
    domain_size = float(args.n_side)

    gen_rng = random.Random(args.seed * 7919 + rank)
    id_counter = [100_000_000 + rank * 1_000_000]

    n_total_trace = []
    for gen in range(args.n_generations):
        for pid in my_patches:
            box = patches[pid]
            new = ml.gen_particles_for_patch(pid, box, gen_rng, n_interior=args.n_inject_per_patch,
                                              n_boundary_cluster=0)
            local[pid].extend(new)

        local, n_merged = cell_merge_round(local, my_patches, args.cell_dx, args.thresh, id_counter)

        n_local = sum(len(local[pid]) for pid in my_patches)
        n_global = comm.allreduce(n_local, op=MPI.SUM)
        if rank == 0:
            n_total_trace.append(n_global)

    final_local = [p for pid in local for p in local[pid]]
    all_final = comm.gather(final_local, root=0)

    if rank == 0:
        final_flat = [p for sub in all_final for p in sub]
        n_final = len(final_flat)
        final_positions = [p["pos"] for p in final_flat]

        null_rng = random.Random(args.seed + 424242)
        null_positions = [(null_rng.uniform(0, domain_size), null_rng.uniform(0, domain_size))
                           for _ in range(n_final)]

        n_hist_bins = int(round(domain_size / args.cell_dx * args.hist_bins_per_cell))
        bin_width = domain_size / n_hist_bins

        hist_final, _, _ = np.histogram2d(
            [p[0] for p in final_positions], [p[1] for p in final_positions],
            bins=n_hist_bins, range=[[0, domain_size], [0, domain_size]])
        hist_null, _, _ = np.histogram2d(
            [p[0] for p in null_positions], [p[1] for p in null_positions],
            bins=n_hist_bins, range=[[0, domain_size], [0, domain_size]])

        k_final, p_final = radial_power_spectrum(hist_final, bin_width)
        k_null, p_null = radial_power_spectrum(hist_null, bin_width)

        report = {
            "n_side": args.n_side, "n_generations": args.n_generations,
            "n_inject_per_patch": args.n_inject_per_patch, "n_final": n_final,
            "n_total_trace": n_total_trace, "k_patch": 1.0, "k_cell": 1.0 / args.cell_dx,
            "k_values": k_final.tolist(), "power_final": p_final.tolist(),
            "power_null": p_null.tolist(),
        }
        with open(args.out, "w") as f:
            json.dump(report, f, indent=2)

        print(f"n_final={n_final} n_total_trace(last 8)={n_total_trace[-8:]}")


if __name__ == "__main__":
    main()
