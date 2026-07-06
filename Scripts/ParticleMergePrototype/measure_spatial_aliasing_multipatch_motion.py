"""
Same setup as measure_spatial_aliasing_multipatch.py (n_side x n_side unit patches, one rank per
patch, repeated inject-then-merge generations), but now with a small per-round random walk applied
to every SURVIVING particle before each generation's merge -- a minimal stand-in for real particle
motion between timesteps. Tests whether the strong patch-frequency spike found in the static
(no-motion) version survives, shrinks, or disappears once particles don't sit at a fixed position
relative to a patch boundary generation after generation.

A particle that drifts across a patch boundary is migrated to its new owning rank (point-location
+ one alltoall, mirroring the placement mechanism already used elsewhere in the protocol for newly
merged particles) BEFORE that generation's merge round runs, so ownership is always correct by the
time run_one_round() does its own ghost-fill.
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


def reflect(v, lo, hi):
    # Simple elastic reflection off the domain's true outer edge, so a random walk never leaks
    # particles out of the finite domain (which would just be population loss, not a meaningful
    # part of what we're testing).
    span = hi - lo
    if span <= 0:
        return lo
    v = v - lo
    v = v % (2 * span)
    if v > span:
        v = 2 * span - v
    v = v + lo
    return min(max(v, lo), hi - 1e-9)


def migrate(comm, rank, size, patches, my_patches, local, domain_size, drift_sigma, rng):
    """Displace every local particle by a small random walk, reflecting off the domain edge, then
    ship anything that crossed into a different rank's patch to its new owner."""
    outgoing = [[] for _ in range(size)]
    new_local = {pid: [] for pid in my_patches}

    for pid in my_patches:
        for p in local[pid]:
            x, y = p["pos"]
            x = reflect(x + rng.gauss(0, drift_sigma), 0.0, domain_size)
            y = reflect(y + rng.gauss(0, drift_sigma), 0.0, domain_size)
            new_pid = ml.locate_patch((x, y), patches)
            new_rank = patches[new_pid][2]
            p = dict(p)
            p["pos"] = (x, y)
            p["rank"] = new_rank
            p["patch"] = new_pid
            if new_rank == rank:
                new_local[new_pid].append(p)
            else:
                outgoing[new_rank].append(p)

    incoming = comm.alltoall(outgoing)
    for batch in incoming:
        for p in batch:
            new_local[p["patch"]].append(p)

    return new_local


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-side", type=int, default=5)
    ap.add_argument("--cell-dx", type=float, default=0.2)
    ap.add_argument("--thresh", type=int, default=2)
    ap.add_argument("--ghost-width", type=float, default=0.15)
    ap.add_argument("--n-inject-per-patch", type=int, default=40)
    ap.add_argument("--n-generations", type=int, default=25)
    ap.add_argument("--drift-sigma", type=float, default=0.05)
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
    shuffle_rng = random.Random(args.seed + 999 + rank)
    drift_rng = random.Random(args.seed * 31337 + rank)
    id_counter = [100_000_000 + rank * 1_000_000]

    n_total_trace = []
    for gen in range(args.n_generations):
        if args.drift_sigma > 0 and gen > 0:
            local = migrate(comm, rank, size, patches, my_patches, local, domain_size,
                             args.drift_sigma, drift_rng)

        for pid in my_patches:
            box = patches[pid]
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

        k_patch = 1.0
        k_cell = 1.0 / args.cell_dx

        report = {
            "n_side": args.n_side, "n_generations": args.n_generations,
            "drift_sigma": args.drift_sigma,
            "n_inject_per_patch": args.n_inject_per_patch, "n_final": n_final,
            "n_total_trace": n_total_trace, "k_patch": k_patch, "k_cell": k_cell,
            "k_values": k_final.tolist(), "power_final": p_final.tolist(),
            "power_null": p_null.tolist(),
        }
        with open(args.out, "w") as f:
            json.dump(report, f, indent=2)

        print(f"drift_sigma={args.drift_sigma} n_final={n_final} n_total_trace(last 8)={n_total_trace[-8:]}")


if __name__ == "__main__":
    main()
