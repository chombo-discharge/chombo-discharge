"""
Single-round drain test: seed EXACTLY 64 particles per cell (not a randomly-injected steady
state), set thresh=16 (eligible while count>16, so the "ideal" target is draining down to
exactly 16), run ONE round of the real protocol (run_one_round_2d -- trivial tier, then
propose/judge/verdict, then STOP -- no second round, no repeated generations), and report what
actually remains per cell. Split by interior (non-exposed) vs. boundary-exposed cells, since the
earlier completion-rate measurement found those behave very differently per round (~93% vs ~43%
completion) -- this test checks whether that gap holds up under a much larger, one-shot backlog.
"""

import argparse
import json
import random
import sys

from mpi4py import MPI

import merge_lib_2d_periodic as ml
from run_scenario_2d_periodic import run_one_round_2d, patch_to_rank


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-patches-x", type=int, default=5)
    ap.add_argument("--n-patches-y", type=int, default=5)
    ap.add_argument("--patch-size", type=float, default=1.0)
    ap.add_argument("--cell-dx", type=float, default=0.2)
    ap.add_argument("--n-per-cell", type=int, default=64)
    ap.add_argument("--thresh", type=int, default=16)
    ap.add_argument("--ghost-width", type=float, default=0.15)
    ap.add_argument("--iterate-local-to-convergence", action="store_true",
                     help="loop the local trivial-tier search/resolve step to convergence within "
                          "this one round, before the single cross-patch propose/judge/verdict "
                          "exchange -- see RULES.md, 'Optional: iterate the local tier...'")
    ap.add_argument("--max-fallback-candidates", type=int, default=0,
                     help="when a particle's chosen candidate is blocked, retry up to this many "
                          "next-nearest still-available candidates before giving up")
    ap.add_argument("--max-cell-distance", type=int, default=None,
                     help="physical merge-eligibility cap in whole grid cells (Chebyshev cell-index "
                          "distance); None (default) = unrestricted, matching original behavior")
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
    cells_per_patch_side = int(round(args.patch_size / args.cell_dx))

    gen_rng = random.Random(args.seed * 7919 + rank)
    shuffle_rng = random.Random(args.seed + 999 + rank)
    id_counter = [100_000_000 + rank * 1_000_000]

    # ---- Seed EXACTLY n_per_cell particles per cell, uniformly within each cell's bounds ----
    local = {pid: [] for pid in my_patches}
    exposed_by_cell = {}  # (pid, cell_i, cell_j) -> is this cell's CENTER exposed
    for pid in my_patches:
        i, j = pid % args.n_patches_x, pid // args.n_patches_x
        patch_lo_x, patch_lo_y = i * args.patch_size, j * args.patch_size
        for ci in range(cells_per_patch_side):
            for cj in range(cells_per_patch_side):
                cell_lo_x = patch_lo_x + ci * args.cell_dx
                cell_lo_y = patch_lo_y + cj * args.cell_dx
                cell_center = (cell_lo_x + args.cell_dx / 2, cell_lo_y + args.cell_dx / 2)
                is_exposed = ml.is_boundary_exposed(cell_center, i, j, args.ghost_width,
                                                     args.patch_size, args.n_patches_x, args.n_patches_y)
                exposed_by_cell[(pid, ci, cj)] = is_exposed
                for _ in range(args.n_per_cell):
                    pos = (gen_rng.uniform(cell_lo_x + 1e-9, cell_lo_x + args.cell_dx - 1e-9),
                           gen_rng.uniform(cell_lo_y + 1e-9, cell_lo_y + args.cell_dx - 1e-9))
                    id_counter[0] += 1
                    local[pid].append(ml.make_particle(pos, 1.0, pid, pid, gen_rng, id_counter[0]))

    n_cells_total = len(exposed_by_cell)
    n_exposed_cells = sum(1 for v in exposed_by_cell.values() if v)
    n_local_initial = sum(len(local[pid]) for pid in my_patches)
    n_initial_global = comm.allreduce(n_local_initial, op=MPI.SUM)

    # ---- Run the protocol EXACTLY ONCE -- no second round ----
    local, n_trivial, n_judge = run_one_round_2d(
        comm, rank, size, args.n_patches_x, args.n_patches_y, args.patch_size, my_patches, local,
        args.thresh, args.ghost_width, args.cell_dx, shuffle_rng, id_counter,
        iterate_local_to_convergence=args.iterate_local_to_convergence,
        max_fallback_candidates=args.max_fallback_candidates,
        max_cell_distance=args.max_cell_distance)

    n_trivial_g = comm.reduce(n_trivial, op=MPI.SUM, root=0)
    n_judge_g = comm.reduce(n_judge, op=MPI.SUM, root=0)

    # ---- Classify each FINAL particle's cell as interior or exposed (by cell center, matching
    # the classification used at seeding time) and tally counts per cell ----
    my_final = [p for pid in local for p in local[pid]]
    final_counts_interior = {}
    final_counts_exposed = {}
    for pid in my_patches:
        i, j = pid % args.n_patches_x, pid // args.n_patches_x
        for p in local[pid]:
            ci = int((p["pos"][0] - i * args.patch_size) / args.cell_dx)
            cj = int((p["pos"][1] - j * args.patch_size) / args.cell_dx)
            ci = min(max(ci, 0), cells_per_patch_side - 1)
            cj = min(max(cj, 0), cells_per_patch_side - 1)
            key = (pid, ci, cj)
            is_exposed = exposed_by_cell.get(key, False)
            d = final_counts_exposed if is_exposed else final_counts_interior
            d[key] = d.get(key, 0) + 1

    all_interior = comm.gather(final_counts_interior, root=0)
    all_exposed = comm.gather(final_counts_exposed, root=0)
    all_exposed_map = comm.gather(exposed_by_cell, root=0)

    if rank == 0:
        n_local_final = len(my_final)
        n_final_global = comm.reduce(n_local_final, op=MPI.SUM, root=0)

        interior_counts = []
        for d in all_interior:
            interior_counts.extend(d.values())
        exposed_counts = []
        for d in all_exposed:
            exposed_counts.extend(d.values())

        # Any interior/exposed cell with NO surviving particles at all (fully drained, or --
        # for interior cells targeting 16, shouldn't happen -- would show as a gap) needs to
        # count as 0, not be silently absent from the list.
        all_keys = {}
        for m in all_exposed_map:
            all_keys.update(m)
        interior_keys = [k for k, v in all_keys.items() if not v]
        exposed_keys = [k for k, v in all_keys.items() if v]

        interior_full = []
        exposed_full = []
        for d in all_interior:
            pass
        merged_interior = {}
        for d in all_interior:
            merged_interior.update(d)
        merged_exposed = {}
        for d in all_exposed:
            merged_exposed.update(d)
        for k in interior_keys:
            interior_full.append(merged_interior.get(k, 0))
        for k in exposed_keys:
            exposed_full.append(merged_exposed.get(k, 0))

        report = {
            "n_per_cell_initial": args.n_per_cell, "thresh": args.thresh,
            "n_cells_total": n_cells_total, "n_exposed_cells": n_exposed_cells,
            "n_interior_cells": n_cells_total - n_exposed_cells,
            "n_initial": n_initial_global, "n_final": n_final_global,
            "n_trivial": n_trivial_g, "n_judge": n_judge_g,
            "interior_final_counts": interior_full,
            "exposed_final_counts": exposed_full,
            "interior_mean": sum(interior_full) / len(interior_full) if interior_full else None,
            "exposed_mean": sum(exposed_full) / len(exposed_full) if exposed_full else None,
        }
        with open(args.out, "w") as f:
            json.dump(report, f, indent=2)

        print(f"n_initial={n_initial_global} n_final={n_final_global} n_trivial={n_trivial_g} n_judge={n_judge_g}")
        if interior_full:
            print(f"interior cells: {len(interior_full)}  mean_final={report['interior_mean']:.3f}  "
                  f"min={min(interior_full)}  max={max(interior_full)}")
        else:
            print("interior cells: none")
        if exposed_full:
            print(f"exposed cells:  {len(exposed_full)}  mean_final={report['exposed_mean']:.3f}  "
                  f"min={min(exposed_full)}  max={max(exposed_full)}")
        else:
            print("exposed cells: none")


if __name__ == "__main__":
    main()
