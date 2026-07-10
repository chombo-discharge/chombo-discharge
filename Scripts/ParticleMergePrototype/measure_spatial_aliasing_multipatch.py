"""
Same spatial-aliasing diagnostic as measure_spatial_aliasing.py, but multi-patch: checks whether
the PATCH TILING itself (not just the merge cell grid) imprints a periodicity onto the surviving
particle distribution after repeated inject-then-merge generations. This is the natural follow-up
to the single-patch cell-frequency check (clean, no spike found) and to the earlier boundary-
retention-bias measurement (near-boundary regions demonstrably merge less per round, converging to
a stable spatial pattern) -- the open question those two results leave is whether THAT retention
bias shows up as excess power at the patch spatial frequency specifically, the way a cell-snapping
artifact would show up at the cell frequency.

Setup: an n_side x n_side grid of unit patches, ONE PATCH PER RANK (worst case -- every internal
boundary is cross-rank, matching the earlier boundary-bias measurement's domain exactly). Each
generation, every rank injects a uniformly-random batch of new particles into its own patch (so
injection is uniform over the whole domain by construction, not just per-patch), then the REAL
multi-rank run_one_round() (ghost-fill, propose, judge, verdict, placement -- unmodified, the
actual validated protocol) runs once. Repeated for many generations. After reaching a stable
population, bin the final particle positions into a fine domain-wide histogram, take its radial
power spectrum, and check specifically for a spike at k_patch = 1/patch_size and harmonics,
relative to both a null (never-merged, same-count uniform) reference and the local smooth trend.
"""

import argparse
import json
import random
import sys

import numpy as np
from mpi4py import MPI

import merge_lib as ml
from run_scenario import run_one_round


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-side", type=int, default=5)
    ap.add_argument("--cell-dx", type=float, default=0.2)
    ap.add_argument("--thresh", type=int, default=2)
    ap.add_argument("--ghost-width", type=float, default=0.15)
    ap.add_argument("--n-inject-per-patch", type=int, default=40)
    ap.add_argument("--n-generations", type=int, default=30)
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

    gen_rng = random.Random(args.seed * 7919 + rank)
    shuffle_rng = random.Random(args.seed + 999 + rank)
    id_counter = [100_000_000 + rank * 1_000_000]

    n_total_trace = []
    for gen in range(args.n_generations):
        for pid in my_patches:
            box = patches[pid]
            lo, hi = box[0], box[1]
            new = ml.gen_particles_for_patch(pid, box, gen_rng, n_interior=args.n_inject_per_patch,
                                              n_boundary_cluster=0)
            local[pid].extend(new)

        local, n_trivial, n_judge = run_one_round(
            comm, rank, size, patches, my_patches, local, args.thresh, args.ghost_width, args.cell_dx,
            shuffle_rng, id_counter)

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
        domain_size = float(args.n_side)
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

        k_patch = 1.0 / 1.0  # patch size is 1 unit by construction
        k_cell = 1.0 / args.cell_dx

        report = {
            "n_side": args.n_side, "n_generations": args.n_generations,
            "n_inject_per_patch": args.n_inject_per_patch, "n_final": n_final,
            "n_total_trace": n_total_trace, "k_patch": k_patch, "k_cell": k_cell,
            "k_values": k_final.tolist(), "power_final": p_final.tolist(),
            "power_null": p_null.tolist(),
        }
        with open(args.out, "w") as f:
            json.dump(report, f, indent=2)

        print(f"n_final={n_final} n_total_trace(last 8)={n_total_trace[-8:]}")
        print(f"k_patch={k_patch:.3f} k_cell={k_cell:.3f}")
        print(f"{'k':>8s} {'power_final':>14s} {'power_null':>12s} {'ratio':>8s}")
        for i in range(len(k_final)):
            marker = ""
            for h in (1, 2, 3):
                if abs(k_final[i] - h * k_patch) < (k_final[1] - k_final[0]) / 2:
                    marker = f" <== patch harmonic {h}"
            ratio = p_final[i] / p_null[i] if p_null[i] > 0 else float("nan")
            print(f"{k_final[i]:8.3f} {p_final[i]:14.2f} {p_null[i]:12.2f} {ratio:8.3f}{marker}")


if __name__ == "__main__":
    main()
