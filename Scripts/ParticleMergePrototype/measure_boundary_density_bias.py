"""
Physical diagnostic, not a correctness test: does repeated merging of an initially UNIFORM
particle distribution introduce a systematic bias near patch/rank boundaries, of the kind that
would corrupt a PIC code's charge density?

Mass conservation is already proven exactly (lineage ledger, 360 runs) -- so charge in any given
region is not at risk. What COULD still be biased is the *coarse-graining*: if boundary-adjacent
cells complete fewer merges per round than interior cells (which the completion-rate measurement
already found -- ~93% vs ~43% per round), boundary regions are left with MORE, SMALLER-weight
particles than the interior, at exactly the same total mass. Same charge, different macroparticle
texture -- more shot noise near boundaries than in the interior. That's the artifact worth
measuring directly, not inferring from the completion-rate proxy.

Method:
  1. Tile a uniform grid of unit patches (default 5x5 = 25 patches), ONE PATCH PER RANK -- the
     worst case (every internal boundary is cross-rank).
  2. Seed particles uniformly at random *within each patch's own box* (gen_particles_for_patch's
     existing per-patch uniform sampling) at equal density everywhere -- so the initial state is
     uniform by construction, not because we hope so.
  3. Run several rounds of the REAL, validated run_one_round() (imported unmodified from
     run_scenario.py) with NO particle motion between rounds.
  4. For both the initial and final particle sets, classify every particle by its distance to the
     nearest INTERNAL patch boundary (excluding a margin near the true outer domain edge, which
     has no neighbor on that side and isn't a fair comparison point).
  5. Bin by that distance and report, per bin: particle count, total mass, and mean weight, for
     both initial and final -- so "final count / initial count" and "final mean weight" as a
     function of distance-to-boundary directly show whether coarse-graining is uneven.
"""

import argparse
import json
import math
import random
import sys

from mpi4py import MPI

import merge_lib as ml
from run_scenario import run_one_round


def make_uniform_grid(n_side):
    """n_side x n_side unit patches tiling [0, n_side] x [0, n_side], one rank per patch."""
    patches = []
    rank = 0
    for j in range(n_side):
        for i in range(n_side):
            patches.append(((float(i), float(j)), (float(i + 1), float(j + 1)), rank))
            rank += 1
    return patches, rank


def dist_to_nearest_internal_boundary(pos, n_side):
    """Minimum distance to any patch grid line, considering ONLY internal lines (1..n_side-1 in
    each axis) -- the domain's own outer edge (0 and n_side) is excluded since it has no neighbor
    on the far side and isn't a fair comparison point for a cross-patch merge-bias question."""
    best = math.inf
    for c in pos:
        for line in range(1, n_side):
            best = min(best, abs(c - line))
    return best


def margin_ok(pos, n_side, margin):
    """True iff pos is far enough from the domain's true outer edge that its "nearest boundary"
    distance isn't contaminated by an edge that has no neighbor at all."""
    return all(margin <= c <= n_side - margin for c in pos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-side", type=int, default=5)
    ap.add_argument("--particles-per-patch", type=int, default=150)
    ap.add_argument("--thresh", type=int, default=2)
    ap.add_argument("--ghost-width", type=float, default=0.15)
    ap.add_argument("--cell-dx", type=float, default=0.2)
    ap.add_argument("--rounds", type=int, default=8)
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

    local = {}
    for pid in my_patches:
        gen_rng = random.Random(args.seed * 1000 + pid)
        box = patches[pid]
        local[pid] = ml.gen_particles_for_patch(pid, box, gen_rng, n_interior=args.particles_per_patch,
                                                  n_boundary_cluster=0)

    initial_local = [dict(p) for pid in my_patches for p in local[pid]]

    shuffle_rng = random.Random(12345 + rank)
    id_counter = [100_000_000 + rank * 1_000_000]

    for r in range(args.rounds):
        local, n_trivial, n_judge = run_one_round(
            comm, rank, size, patches, my_patches, local, args.thresh, args.ghost_width, args.cell_dx,
            shuffle_rng, id_counter)

    final_local = [p for pid in local for p in local[pid]]

    all_initial = comm.gather(initial_local, root=0)
    all_final = comm.gather(final_local, root=0)

    if rank == 0:
        initial_flat = [p for sub in all_initial for p in sub]
        final_flat = [p for sub in all_final for p in sub]

        margin = args.ghost_width  # exclude the outer ring so every classified particle has a
                                    # genuine internal boundary as its nearest one

        def bin_particles(particles, edges):
            bins = [{"count": 0, "mass": 0.0} for _ in range(len(edges) - 1)]
            excluded = 0
            for p in particles:
                if not margin_ok(p["pos"], args.n_side, margin):
                    excluded += 1
                    continue
                d = dist_to_nearest_internal_boundary(p["pos"], args.n_side)
                for i in range(len(edges) - 1):
                    if edges[i] <= d < edges[i + 1]:
                        bins[i]["count"] += 1
                        bins[i]["mass"] += p["weight"]
                        break
            return bins, excluded

        edges = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
        initial_bins, init_excluded = bin_particles(initial_flat, edges)
        final_bins, final_excluded = bin_particles(final_flat, edges)

        report = {"n_side": args.n_side, "seed": args.seed, "rounds": args.rounds,
                  "n_initial": len(initial_flat), "n_final": len(final_flat),
                  "edges": edges, "initial_bins": initial_bins, "final_bins": final_bins,
                  "excluded_initial": init_excluded, "excluded_final": final_excluded}
        with open(args.out, "w") as f:
            json.dump(report, f, indent=2)

        print(f"n_initial={len(initial_flat)} n_final={len(final_flat)} "
              f"overall_retained_count_frac={len(final_flat)/len(initial_flat):.4f}")
        print(f"{'dist_band':>14s} {'init_count':>11s} {'final_count':>12s} {'retained_frac':>14s} "
              f"{'init_mean_w':>12s} {'final_mean_w':>13s}")
        for i in range(len(edges) - 1):
            ib, fb = initial_bins[i], final_bins[i]
            band = f"[{edges[i]:.2f},{edges[i+1]:.2f})"
            retained = fb["count"] / ib["count"] if ib["count"] else float("nan")
            init_mw = ib["mass"] / ib["count"] if ib["count"] else float("nan")
            final_mw = fb["mass"] / fb["count"] if fb["count"] else float("nan")
            print(f"{band:>14s} {ib['count']:11d} {fb['count']:12d} {retained:14.4f} "
                  f"{init_mw:12.4f} {final_mw:13.4f}")


if __name__ == "__main__":
    main()
