"""
2D periodic testbed for the patch-frequency aliasing question -- 3x3 patches, periodic BCs (no
domain-edge confound), corners preserved (unlike the 1D reduction, which turned out to
structurally lack them and showed no spike as a result). Same inject-then-merge-repeatedly
protocol, using run_one_round_2d (an unmodified periodic port of the validated algorithm).

Includes the same basic correctness check (weight + lineage ledger) as the 1D port before
trusting the spectrum -- this is fresh code, not a re-run of already-validated code.
"""

import argparse
import json
import random
import sys

import numpy as np
from mpi4py import MPI

import merge_lib_2d_periodic as ml
from run_scenario_2d_periodic import run_one_round_2d, patch_to_rank


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
    ap.add_argument("--n-patches-x", type=int, default=3)
    ap.add_argument("--n-patches-y", type=int, default=3)
    ap.add_argument("--patch-size", type=float, default=1.0)
    ap.add_argument("--cell-dx", type=float, default=0.2)
    ap.add_argument("--thresh", type=int, default=2)
    ap.add_argument("--ghost-width", type=float, default=0.15)
    ap.add_argument("--n-inject-per-patch", type=int, default=40)
    ap.add_argument("--n-generations", type=int, default=25)
    ap.add_argument("--hist-bins-per-cell", type=int, default=4)
    ap.add_argument("--max-cell-distance", type=int, default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    n_patches = args.n_patches_x * args.n_patches_y
    # Block assignment: n_patches_y rows split into `size` contiguous blocks, so the grid can
    # scale arbitrarily while the rank count stays fixed at whatever the machine actually has.
    if size > args.n_patches_y:
        if rank == 0:
            print(f"size ({size}) must not exceed n_patches_y ({args.n_patches_y})", file=sys.stderr)
        comm.Abort(1)

    my_patches = [pid for pid in range(n_patches)
                  if patch_to_rank(pid, args.n_patches_x, args.n_patches_y, size) == rank]
    Lx = args.patch_size * args.n_patches_x
    Ly = args.patch_size * args.n_patches_y

    gen_rng = random.Random(args.seed * 7919 + rank)
    shuffle_rng = random.Random(args.seed + 999 + rank)
    id_counter = [100_000_000 + rank * 1_000_000]

    local = {pid: [] for pid in my_patches}
    initial_local = []
    n_total_trace = []
    for gen in range(args.n_generations):
        for pid in my_patches:
            i, j = pid % args.n_patches_x, pid // args.n_patches_x
            new = ml.gen_particles_uniform(pid, i, j, args.patch_size, gen_rng, args.n_inject_per_patch, id_counter)
            local[pid].extend(new)
            initial_local.extend(dict(p) for p in new)

        local, n_trivial, n_judge = run_one_round_2d(
            comm, rank, size, args.n_patches_x, args.n_patches_y, args.patch_size, my_patches, local,
            args.thresh, args.ghost_width, args.cell_dx, shuffle_rng, id_counter,
            max_cell_distance=args.max_cell_distance)

        n_local = sum(len(local[pid]) for pid in my_patches)
        n_global = comm.allreduce(n_local, op=MPI.SUM)
        if rank == 0:
            n_total_trace.append(n_global)

    my_final = [p for pid in local for p in local[pid]]
    all_final = comm.gather(my_final, root=0)
    all_initial = comm.gather(initial_local, root=0)

    if rank == 0:
        final_flat = [p for sub in all_final for p in sub]
        initial_flat = [p for sub in all_initial for p in sub]
        n_final = len(final_flat)

        orig_weight = {p["id"]: p["weight"] for p in initial_flat}
        all_orig_ids = set(orig_weight.keys())
        seen_ids = set()
        ok_lineage_disjoint = True
        ok_lineage_weight = True
        for p in final_flat:
            lineage = p["lineage"]
            if lineage & seen_ids:
                ok_lineage_disjoint = False
            seen_ids |= lineage
            lw = sum(orig_weight[k] for k in lineage)
            if abs(lw - p["weight"]) > 1e-9:
                ok_lineage_weight = False
        ok_lineage_complete = (seen_ids == all_orig_ids)
        ids = [p["id"] for p in final_flat]
        ok_ids_unique = len(ids) == len(set(ids))
        total_w_initial = sum(p["weight"] for p in initial_flat)
        total_w_final = sum(p["weight"] for p in final_flat)
        ok_weight = abs(total_w_initial - total_w_final) < 1e-9

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

        k_patch = 1.0 / args.patch_size
        k_cell = 1.0 / args.cell_dx

        report = {
            "n_patches_x": args.n_patches_x, "n_patches_y": args.n_patches_y,
            "n_generations": args.n_generations, "n_inject_per_patch": args.n_inject_per_patch,
            "n_final": n_final, "n_total_trace": n_total_trace,
            "ok_weight": ok_weight, "ok_ids_unique": ok_ids_unique,
            "ok_lineage_disjoint": ok_lineage_disjoint, "ok_lineage_weight": ok_lineage_weight,
            "ok_lineage_complete": ok_lineage_complete,
            "k_patch": k_patch, "k_cell": k_cell,
            "k_values": k_final.tolist(), "power_final": p_final.tolist(), "power_null": p_null.tolist(),
        }
        with open(args.out, "w") as f:
            json.dump(report, f, indent=2)

        status = "PASS" if (ok_weight and ok_ids_unique and ok_lineage_disjoint and ok_lineage_weight
                            and ok_lineage_complete) else "FAIL"
        print(f"[{status}] n_final={n_final} n_total_trace(last 8)={n_total_trace[-8:]} "
              f"weight_ok={ok_weight} ids_unique={ok_ids_unique} lineage_disjoint={ok_lineage_disjoint} "
              f"lineage_weight={ok_lineage_weight} lineage_complete={ok_lineage_complete}")


if __name__ == "__main__":
    main()
