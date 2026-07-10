"""
1D periodic testbed for the patch-frequency aliasing question -- a faster, simpler harness for
iterating on fixes, per the request to reduce to 1D (10 patches, 16 cells/patch, periodic BCs so
there's no domain-edge confound). Same inject-then-merge-repeatedly protocol as the 2D multipatch
aliasing measurement, using the REAL validated algorithm (run_one_round_1d, an unmodified 1D port
of run_scenario.py's protocol -- see merge_lib_1d.py / run_scenario_1d.py).

Includes a basic correctness check (weight + lineage conservation, via the same per-particle
ledger audit used in the 2D validation) before trusting the spectrum -- this is a fresh port, not
just a re-run of already-validated code, so it needs its own sanity check.
"""

import argparse
import json
import random
import sys

import numpy as np
from mpi4py import MPI

import merge_lib_1d as ml
from run_scenario_1d import run_one_round_1d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-patches", type=int, default=10)
    ap.add_argument("--cells-per-patch", type=int, default=16)
    ap.add_argument("--patch-size", type=float, default=1.0)
    ap.add_argument("--thresh", type=int, default=2)
    ap.add_argument("--ghost-width", type=float, default=0.15)
    ap.add_argument("--n-inject-per-patch", type=int, default=40)
    ap.add_argument("--n-generations", type=int, default=25)
    ap.add_argument("--hist-bins-per-cell", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if size != args.n_patches:
        if rank == 0:
            print(f"need {args.n_patches} ranks, got {size}", file=sys.stderr)
        comm.Abort(1)

    cell_dx = args.patch_size / args.cells_per_patch
    L = args.patch_size * args.n_patches

    gen_rng = random.Random(args.seed * 7919 + rank)
    shuffle_rng = random.Random(args.seed + 999 + rank)
    id_counter = [100_000_000 + rank * 1_000_000]

    local = []
    initial_local = []
    n_total_trace = []
    for gen in range(args.n_generations):
        new = ml.gen_particles_uniform(rank, args.patch_size, gen_rng, args.n_inject_per_patch, id_counter)
        local.extend(new)
        initial_local.extend(dict(p) for p in new)

        local, n_trivial, n_judge = run_one_round_1d(
            comm, rank, size, args.n_patches, args.patch_size, local, args.thresh, args.ghost_width,
            cell_dx, shuffle_rng, id_counter)

        n_global = comm.allreduce(len(local), op=MPI.SUM)
        if rank == 0:
            n_total_trace.append(n_global)

    all_final = comm.gather(local, root=0)
    all_initial = comm.gather(initial_local, root=0)

    if rank == 0:
        final_flat = [p for sub in all_final for p in sub]
        initial_flat = [p for sub in all_initial for p in sub]
        n_final = len(final_flat)

        # ---- Basic correctness check: weight + lineage ledger (this is a fresh port) ----
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
            lw = sum(orig_weight[i] for i in lineage)
            if abs(lw - p["weight"]) > 1e-9:
                ok_lineage_weight = False
        ok_lineage_complete = (seen_ids == all_orig_ids)
        ids = [p["id"] for p in final_flat]
        ok_ids_unique = len(ids) == len(set(ids))
        total_w_initial = sum(p["weight"] for p in initial_flat)
        total_w_final = sum(p["weight"] for p in final_flat)
        ok_weight = abs(total_w_initial - total_w_final) < 1e-9

        # ---- Spectral analysis ----
        positions = [p["pos"] for p in final_flat]
        null_rng = random.Random(args.seed + 424242)
        null_positions = [null_rng.uniform(0, L) for _ in range(n_final)]

        n_hist_bins = int(round(L / cell_dx * args.hist_bins_per_cell))
        bin_width = L / n_hist_bins

        hist_final, _ = np.histogram(positions, bins=n_hist_bins, range=(0, L))
        hist_null, _ = np.histogram(null_positions, bins=n_hist_bins, range=(0, L))

        def power_spectrum_1d(hist, bin_width):
            f = np.fft.fft(hist - hist.mean())
            power = np.abs(f) ** 2
            freqs = np.fft.fftfreq(len(hist), d=bin_width)
            pos_mask = freqs > 0
            return freqs[pos_mask], power[pos_mask]

        k_final, p_final = power_spectrum_1d(hist_final, bin_width)
        k_null, p_null = power_spectrum_1d(hist_null, bin_width)

        k_patch = 1.0 / args.patch_size
        k_cell = 1.0 / cell_dx

        report = {
            "n_patches": args.n_patches, "cells_per_patch": args.cells_per_patch,
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
        print(f"k_patch={k_patch} k_cell={k_cell}")


if __name__ == "__main__":
    main()
