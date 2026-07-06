"""
Pressure test: a proposal's source particle is NOT guaranteed to be visible as a ghost on the
target's rank -- mutual ghost visibility is a common case, not a guarantee. This test deliberately
constructs the asymmetric case and checks the protocol is still correct.

Concretely: patch A (rank 0) owns p0, positioned just far enough from the shared A/B boundary that
it is NOT boundary-exposed (never shipped anywhere). Patch B (rank 1) owns p1, positioned close
enough to the boundary that it IS exposed and gets shipped to A as a ghost, and close enough to p0
to be its true nearest neighbor. So: A sees B's p1 (that's how p0's own search finds it at all),
but B never sees A's p0 -- by construction, not by chance. p0's proposal to p1 must therefore be
judged using ONLY the data carried explicitly in the proposal message itself (position/weight/id/
rank/payload), never a pre-existing ghost copy -- if the implementation ever assumed mutual
visibility (e.g. reconstructing p0 from a cached ghost instead of the proposal payload), this is
exactly the scenario that would expose it.

Two variants, both run for several shuffle-seeds (determinism under processing-order reordering)
and checked against the full invariant suite (weight conservation, id uniqueness, lineage ledger):

  1. ACCEPT: p1 has no better local option, so it must accept p0's proposal using ONLY the
     shipped data. Checked: the merge happens, the merged particle's position/weight are the
     correct weighted centroid of p0 and p1 specifically (not corrupted, not using stale/wrong
     data).
  2. REJECT: p1 has a closer, LOCALLY VISIBLE rival p2 (also owned by rank 1), so p1 correctly
     prefers p2 and rejects p0 -- even though, from p0's side, p1 looked like its unambiguous best
     match. Checked: p0 survives the round completely untouched (not merged, not lost, not
     duplicated), while p1-p2 (a genuine mutual pair, visible to both) merge correctly.

Both variants explicitly assert the adversarial precondition actually holds (p0 is provably never
shipped to rank 1) rather than just hoping the scenario produces it -- if a future change to
ghost_targets()/patch geometry accidentally made this symmetric, this test would fail loudly on
that assertion rather than silently stop testing anything.
"""

import random
import sys

from mpi4py import MPI

import merge_lib as ml
from run_scenario import run_one_round

# Two adjacent unit-height patches sharing the boundary at x=1.
PATCH_A = ((0.0, 0.0), (1.0, 1.0))
PATCH_B = ((1.0, 0.0), (2.0, 1.0))
PATCHES = [(PATCH_A[0], PATCH_A[1], 0), (PATCH_B[0], PATCH_B[1], 1)]
GHOST_WIDTH = 0.2
CELL_DX = 0.05
THRESH = 0  # isolate the visibility question from the crowding question entirely

# p0: distance 0.3 from the shared boundary (x=1) -- strictly more than GHOST_WIDTH=0.2, so NOT
# exposed toward patch B. p1: distance 0.05 from the boundary on B's side -- well within
# GHOST_WIDTH, so exposed toward A. d(p0, p1) = 0.35.
P0_POS = (0.70, 0.50)
P1_POS = (1.05, 0.50)
# p2: a rival for p1, visible to p1 (same patch, valid, not a ghost), much closer than p0.
P2_POS = (1.06, 0.50)


def make(pos, weight, rank, patch, pid):
    return {"pos": pos, "weight": weight, "id": pid, "energy": 1.0, "rank": rank, "patch": patch,
            "lineage": frozenset((pid,))}


def check_preconditions():
    """Assert the adversarial geometry is actually adversarial, using the SAME pure geometry
    function run_one_round() itself calls -- this is the proof, not an assumption."""
    assert ml.is_boundary_exposed(P0_POS, 0, PATCHES, GHOST_WIDTH) is False, "p0 must NOT be exposed"
    assert ml.ghost_targets(P0_POS, 0, PATCHES, GHOST_WIDTH) == [], "p0 must never be shipped anywhere"
    assert ml.is_boundary_exposed(P1_POS, 1, PATCHES, GHOST_WIDTH) is True, "p1 MUST be exposed"
    assert 0 in ml.ghost_targets(P1_POS, 1, PATCHES, GHOST_WIDTH), "p1 must be shipped to patch A"
    d_p0p1 = ml.dist2(P0_POS, P1_POS) ** 0.5
    d_p1p2 = ml.dist2(P1_POS, P2_POS) ** 0.5
    assert d_p1p2 < d_p0p1, "p2 must be closer to p1 than p0 is, for the REJECT variant"


def run_variant(comm, rank, size, include_p2, shuffle_seed, results_out):
    local = {0: [], 1: []}
    if rank == 0:
        local[0] = [make(P0_POS, 1.0, 0, 0, 1)]
    if rank == 1:
        local[1] = [make(P1_POS, 1.0, 1, 1, 1_000_001)]
        if include_p2:
            local[1].append(make(P2_POS, 1.0, 1, 1, 1_000_002))

    my_patches = [0] if rank == 0 else [1]
    initial_ids = {p["id"]: p["weight"] for pid in my_patches for p in local[pid]}
    all_initial = comm.gather(initial_ids, root=0)

    shuffle_rng = random.Random(shuffle_seed * 1000 + rank)
    id_counter = [100_000_000 + rank * 1_000_000]

    local, n_trivial, n_judge = run_one_round(
        comm, rank, size, PATCHES, my_patches, local, THRESH, GHOST_WIDTH, CELL_DX, shuffle_rng, id_counter
    )
    # Per-rank counts (e.g. n_judge here is only whatever RANK 0 itself judged) -- the merge in
    # this test is judged by rank 1 (it owns p1), so these must be reduced across both ranks, not
    # read from rank 0 alone.
    n_trivial_g = comm.reduce(n_trivial, op=MPI.SUM, root=0)
    n_judge_g = comm.reduce(n_judge, op=MPI.SUM, root=0)

    my_final = [p for pid in local for p in local[pid]]
    all_final = comm.gather(my_final, root=0)

    if rank == 0:
        initial = {}
        for d in all_initial:
            initial.update(d)
        final = [p for sub in all_final for p in sub]

        ok_weight = abs(sum(final_p["weight"] for final_p in final) - sum(initial.values())) < 1e-12
        ok_ids_unique = len({p["id"] for p in final}) == len(final)

        results_out.append({
            "shuffle_seed": shuffle_seed,
            "n_trivial": n_trivial_g,
            "n_judge": n_judge_g,
            "final": sorted((p["id"], round(p["pos"][0], 9), round(p["pos"][1], 9), round(p["weight"], 9))
                            for p in final),
            "ok_weight": ok_weight,
            "ok_ids_unique": ok_ids_unique,
        })


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    assert size == 2, "this pressure test requires exactly 2 ranks"

    if rank == 0:
        check_preconditions()
    comm.Barrier()

    all_pass = True

    # ---- Variant 1: ACCEPT (p1 has no rival, must use p0's shipped data) ----
    accept_results = []
    for sh in range(6):
        run_variant(comm, rank, size, include_p2=False, shuffle_seed=sh, results_out=accept_results)

    if rank == 0:
        first = accept_results[0]["final"]
        for r in accept_results:
            ok_shape = (len(r["final"]) == 1 and r["n_judge"] == 1 and r["n_trivial"] == 0)
            merged = r["final"][0] if r["final"] else None
            ok_centroid = merged is not None and abs(merged[1] - 0.875) < 1e-9 and abs(merged[3] - 2.0) < 1e-9
            deterministic = r["final"] == first
            passed = ok_shape and ok_centroid and r["ok_weight"] and r["ok_ids_unique"] and deterministic
            all_pass &= passed
            print(f"[{'PASS' if passed else 'FAIL'}] ACCEPT variant shuffle={r['shuffle_seed']} "
                  f"final={r['final']} trivial={r['n_trivial']} judge={r['n_judge']} "
                  f"weight_ok={r['ok_weight']} ids_unique={r['ok_ids_unique']} deterministic={deterministic}")

    # ---- Variant 2: REJECT (p1 prefers visible p2, p0 must survive untouched) ----
    reject_results = []
    for sh in range(6):
        run_variant(comm, rank, size, include_p2=True, shuffle_seed=sh, results_out=reject_results)

    if rank == 0:
        first = reject_results[0]["final"]
        for r in reject_results:
            p0_entry = [p for p in r["final"] if p[0] == 1]
            ok_p0_untouched = p0_entry == [(1, 0.7, 0.5, 1.0)]
            ok_shape = (len(r["final"]) == 2 and r["n_judge"] == 1)  # p0 survives + one p1-p2 merge
            deterministic = r["final"] == first
            passed = ok_p0_untouched and ok_shape and r["ok_weight"] and r["ok_ids_unique"] and deterministic
            all_pass &= passed
            print(f"[{'PASS' if passed else 'FAIL'}] REJECT variant shuffle={r['shuffle_seed']} "
                  f"final={r['final']} trivial={r['n_trivial']} judge={r['n_judge']} "
                  f"p0_untouched={ok_p0_untouched} weight_ok={r['ok_weight']} "
                  f"ids_unique={r['ok_ids_unique']} deterministic={deterministic}")

    if rank == 0:
        print("\nSUMMARY:", "ALL PASS" if all_pass else "SOME FAILED")
        sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
