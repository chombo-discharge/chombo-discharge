"""
One MPI-rank's worth of the merge protocol, run under `mpirun -np K`, for --rounds consecutive
rounds (a stand-in for consecutive simulation timesteps).

Per-round structure (all ranks execute the same code, SPMD):
  1. Ghost-fill (alltoall #1) -- ALWAYS fresh, every round, before ANY merge decision. Never
     reused across rounds: round 2 must see round 1's merged particles (different weight, a
     different birth patch) as ordinary ghosts, not stale round-1 data.
  2. Per patch: trivial-tier merges (both endpoints non-exposed) resolved immediately, in a
     shuffled-order-robust way; everything else becomes an outgoing proposal.
  3. Proposals routed to the target's owner (alltoall #2) -- always by partner ownership.
  4. Judge processing, locally, per rank, over all incoming proposals for particles it owns.
  5. Verdicts routed back to proposers (alltoall #3).
  6. Proposers apply verdicts (remove particle if won).
  7. Placement: newly merged particles point-located against the global patch registry;
     anything not in one of my own patches goes out via one more alltoall (#4).

The fresh-id counter is a single monotonic counter that persists across ALL rounds within a run
(never reset), so a round-2 merged particle can never collide with a round-1 one.

--shuffle-seed perturbs local processing order (patch order, particle order, and the order the
received proposal list is processed in) WITHOUT changing the physical particle configuration
(--seed controls that). Running the same --seed with many different --shuffle-seed values and
diffing the final states is the race-freedom check.
"""

import argparse
import json
import random
import sys

from mpi4py import MPI

import merge_lib as ml


def run_one_round(comm, rank, size, patches, my_patches, local, thresh, ghost_width, cell_dx, shuffle_rng, id_counter):
    # ---- Ghost-fill (alltoall #1), fresh every round, before any decision ----
    outgoing_ghosts = [[] for _ in range(size)]
    for pid in my_patches:
        for p in local[pid]:
            for dest_patch in ml.ghost_targets(p["pos"], pid, patches, ghost_width):
                dest_rank = patches[dest_patch][2]
                outgoing_ghosts[dest_rank].append((p, dest_patch))
    incoming_ghosts = comm.alltoall(outgoing_ghosts)

    ghosts = {pid: [] for pid in my_patches}
    for items in incoming_ghosts:
        for p, dest_patch in items:
            ghosts[dest_patch].append(p)

    # ---- Per-patch: trivial tier + candidate generation ----
    alive_valid = {}
    exposed = {}
    live_count = {}
    has_outgoing = {}
    outgoing_target = {}   # particle id -> id of the particle it proposed to (for mutual-match detection)
    trivial_merges = []
    outgoing_proposals = [[] for _ in range(size)]

    patch_order = my_patches[:]
    shuffle_rng.shuffle(patch_order)

    for pid in patch_order:
        valid = local[pid]
        alive_valid[pid] = {p["id"] for p in valid}
        exposed_here = {p["id"]: ml.is_boundary_exposed(p["pos"], pid, patches, ghost_width)
                         for p in valid}
        exposed.update(exposed_here)

        lc = {}
        for p in valid:
            lc[ml.cell_key(p["pos"], cell_dx)] = lc.get(ml.cell_key(p["pos"], cell_dx), 0) + 1
        for g in ghosts[pid]:
            lc[ml.cell_key(g["pos"], cell_dx)] = lc.get(ml.cell_key(g["pos"], cell_dx), 0) + 1
        live_count[pid] = lc

        ghost_ids = {g["id"] for g in ghosts[pid]}
        pool = valid + ghosts[pid]

        order = list(range(len(valid)))
        shuffle_rng.shuffle(order)

        edges = []
        for idx in order:
            p = valid[idx]
            if lc.get(ml.cell_key(p["pos"], cell_dx), 0) <= thresh:
                continue
            nbr, d2 = ml.nearest_neighbor(p, pool, alive_valid[pid] | ghost_ids)
            if nbr is None:
                continue
            edges.append((d2, p, nbr))

        edges.sort(key=lambda e: (e[0], e[1]["id"], e[2]["id"]))

        for d2, p, nbr in edges:
            if p["id"] not in alive_valid[pid]:
                continue
            if nbr["patch"] == pid:
                if nbr["id"] not in alive_valid[pid]:
                    continue
                if has_outgoing.get(nbr["id"], False):
                    continue
                trivial_ok = (not exposed_here[p["id"]]) and (not exposed.get(nbr["id"], False))
            else:
                trivial_ok = False
            ka, kb = ml.cell_key(p["pos"], cell_dx), ml.cell_key(nbr["pos"], cell_dx)
            if lc.get(ka, 0) <= thresh or lc.get(kb, 0) <= thresh:
                continue
            if trivial_ok:
                alive_valid[pid].discard(p["id"])
                alive_valid[pid].discard(nbr["id"])
                lc[ka] = lc.get(ka, 0) - 1
                lc[kb] = lc.get(kb, 0) - 1
                trivial_merges.append((p, nbr, pid))
            else:
                has_outgoing[p["id"]] = True
                outgoing_target[p["id"]] = nbr["id"]
                dest_rank = nbr["rank"]
                outgoing_proposals[dest_rank].append({"source": p, "target_id": nbr["id"], "d2": d2})

    # ---- Proposals routed by target ownership (alltoall #2) ----
    incoming_proposals = comm.alltoall(outgoing_proposals)
    flat_incoming = [prop for bucket in incoming_proposals for prop in bucket]
    shuffle_rng.shuffle(flat_incoming)

    # ---- Judge processing (local, per rank) ----
    by_target = {}
    for prop in flat_incoming:
        by_target.setdefault(prop["target_id"], []).append(prop)

    judge_merges = []
    verdicts_out = [[] for _ in range(size)]

    target_lookup = {}
    for pid in my_patches:
        for p in local[pid]:
            if p["id"] in alive_valid[pid]:
                target_lookup[p["id"]] = (pid, p)

    # Phase 1: for each target, determine its candidate winner (mutual-match or plain argmin)
    # and the distance to sort on -- WITHOUT committing anything yet. Targets with no viable
    # candidate at all (not ours, or has_outgoing with no reciprocal match) are rejected
    # immediately, since nothing about them can ever change from further sorting.
    candidates = []  # (d2, target_id, winner_prop, tpid, tparticle)
    deferred_props = {}  # target_id -> props, only for targets actually added to candidates
    for target_id, props in by_target.items():
        if target_id not in target_lookup:
            for prop in props:
                verdicts_out[prop["source"]["rank"]].append({"proposer_id": prop["source"]["id"], "accepted": False})
            continue
        tpid, tparticle = target_lookup[target_id]
        if has_outgoing.get(target_id, False):
            # target_id has its own outgoing commitment elsewhere -- normally reject everyone
            # (it's "spoken for"). EXCEPT: if target_id's own outgoing proposal was TO one of
            # its current incoming proposers, this is a genuinely mutual best-match pair (both
            # sides' own independent, full-information searches agree), which is unambiguous
            # and safe to accept -- rejecting it unconditionally would starve it forever, since
            # nothing about a stable mutual pair changes from one round to the next. Both sides
            # detect the same mutual match independently, so a deterministic tiebreak (lower id
            # performs the accept) is required to avoid BOTH sides accepting it and double
            # merging the same two particles.
            my_target = outgoing_target.get(target_id)
            mutual = [pr for pr in props if pr["source"]["id"] == my_target] if my_target is not None else []
            if mutual and target_id < my_target:
                candidates.append((mutual[0]["d2"], target_id, mutual[0], tpid, tparticle))
                deferred_props[target_id] = props
            else:
                for prop in props:
                    verdicts_out[prop["source"]["rank"]].append({"proposer_id": prop["source"]["id"], "accepted": False})
            continue
        winner = min(props, key=lambda pr: (pr["d2"], pr["source"]["id"]))
        candidates.append((winner["d2"], target_id, winner, tpid, tparticle))
        deferred_props[target_id] = props

    # Phase 2: process ALL targets (mutual-match and plain alike) in one globally sorted,
    # closest-first order -- exactly the same discipline as the local edge processing above.
    # This matters: accepting one target decrements its cell's live count, which can push a
    # DIFFERENT target sharing that cell below threshold. Iterating by_target in whatever order
    # the (shuffled) incoming-proposal list happened to produce made that outcome shuffle-order
    # dependent; sorting by distance with an id tiebreak makes it not.
    candidates.sort(key=lambda c: (c[0], c[1]))
    accepted_targets = {}  # target_id -> winner_prop, for phase 3's verdict emission
    for d2, target_id, winner, tpid, tparticle in candidates:
        if target_id not in alive_valid[tpid]:
            continue
        ka = ml.cell_key(tparticle["pos"], cell_dx)
        if live_count[tpid].get(ka, 0) <= thresh:
            continue
        alive_valid[tpid].discard(target_id)
        live_count[tpid][ka] = live_count[tpid].get(ka, 0) - 1
        judge_merges.append((tparticle, winner["source"], tpid))
        accepted_targets[target_id] = winner

    # Phase 3: emit verdicts only for targets that were actually deferred to phase 2 (targets
    # rejected outright in phase 1 already got their verdicts there -- iterating by_target again
    # here would emit harmless but wasteful duplicate reject verdicts for them).
    for target_id, props in deferred_props.items():
        winner = accepted_targets.get(target_id)
        for prop in props:
            won = winner is not None and prop is winner
            verdicts_out[prop["source"]["rank"]].append({"proposer_id": prop["source"]["id"], "accepted": won})

    # ---- Verdicts routed back to proposers (alltoall #3) ----
    incoming_verdicts = comm.alltoall(verdicts_out)
    flat_verdicts = [v for bucket in incoming_verdicts for v in bucket]
    accepted_source_ids = {v["proposer_id"] for v in flat_verdicts if v["accepted"]}

    for pid in my_patches:
        alive_valid[pid] -= accepted_source_ids

    # ---- Assemble final valid lists, place newly-created merged particles ----
    def fresh_id():
        id_counter[0] += 1
        return id_counter[0]

    final_local = {pid: [p for p in local[pid] if p["id"] in alive_valid[pid]] for pid in my_patches}

    pending_placement = []
    for a, b, pid in trivial_merges:
        pending_placement.append(ml.combine(a, b, fresh_id(), rank, None))
    for target_p, winner_p, pid in judge_merges:
        pending_placement.append(ml.combine(target_p, winner_p, fresh_id(), rank, None))

    outgoing_placement = [[] for _ in range(size)]
    for m in pending_placement:
        dest_patch = ml.locate_patch(m["pos"], patches)
        if dest_patch is None:
            dest_patch = my_patches[0]
        dest_rank = patches[dest_patch][2]
        m["patch"] = dest_patch
        m["rank"] = dest_rank
        if dest_rank == rank:
            final_local.setdefault(dest_patch, [])
            final_local[dest_patch].append(m)
        else:
            outgoing_placement[dest_rank].append(m)

    incoming_placement = comm.alltoall(outgoing_placement)
    for items in incoming_placement:
        for m in items:
            final_local.setdefault(m["patch"], [])
            final_local[m["patch"]].append(m)

    return final_local, len(trivial_merges), len(judge_merges)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--shuffle-seed", type=int, default=0)
    ap.add_argument("--thresh", type=int, default=2)
    ap.add_argument("--ghost-width", type=float, default=0.15)
    ap.add_argument("--cell-dx", type=float, default=0.2)
    ap.add_argument("--rounds", type=int, default=1)
    ap.add_argument("--n-interior", type=int, default=5)
    ap.add_argument("--n-cluster", type=int, default=3)
    ap.add_argument("--cluster-radius", type=float, default=0.08)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    patches, nranks, stress_points = ml.SCENARIOS[args.scenario]()
    if size != nranks:
        if rank == 0:
            print(f"scenario {args.scenario} requires {nranks} ranks, got {size}", file=sys.stderr)
        comm.Abort(1)

    my_patches = [pid for pid, box in enumerate(patches) if box[2] == rank]

    # ---- Particle generation (deterministic given --seed; independent of shuffle-seed) ----
    local = {}
    for pid in my_patches:
        gen_rng = random.Random(args.seed * 1000 + pid)
        box = patches[pid]
        parts = ml.gen_particles_for_patch(pid, box, gen_rng, n_interior=args.n_interior, n_boundary_cluster=0)
        for spt in stress_points:
            if ml.point_in_box(spt, box, grow=1e-9):
                parts += ml.gen_particles_for_patch(pid, box, gen_rng, n_interior=0,
                                                      n_boundary_cluster=args.n_cluster,
                                                      boundary_point=spt, cluster_radius=args.cluster_radius)
        local[pid] = parts

    initial_total_weight = sum(p["weight"] for pid in my_patches for p in local[pid])
    initial_total_weight = comm.allreduce(initial_total_weight, op=MPI.SUM)
    n_initial_local = sum(len(local[pid]) for pid in my_patches)
    n_initial_global = comm.allreduce(n_initial_local, op=MPI.SUM)
    initial_flat_local = [dict(p) for pid in my_patches for p in local[pid]]

    shuffle_rng = random.Random(args.shuffle_seed * 1000 + rank)
    id_counter = [100_000_000 + rank * 1_000_000]

    round_stats = []
    for r in range(args.rounds):
        local, n_trivial, n_judge = run_one_round(
            comm, rank, size, patches, my_patches, local, args.thresh, args.ghost_width, args.cell_dx,
            shuffle_rng, id_counter)

        n_trivial_g = comm.reduce(n_trivial, op=MPI.SUM, root=0)
        n_judge_g = comm.reduce(n_judge, op=MPI.SUM, root=0)
        n_local = sum(len(local[pid]) for pid in my_patches)
        n_global = comm.allreduce(n_local, op=MPI.SUM)

        # Global count of cells still over threshold, to show monotonic convergence across rounds.
        lc = {}
        for pid in my_patches:
            for p in local[pid]:
                lc[ml.cell_key(p["pos"], args.cell_dx)] = lc.get(ml.cell_key(p["pos"], args.cell_dx), 0) + 1
        n_over = sum(1 for c in lc.values() if c > args.thresh)
        n_over_g = comm.reduce(n_over, op=MPI.SUM, root=0)

        if rank == 0:
            round_stats.append({"round": r, "n_particles": n_global, "n_trivial": n_trivial_g,
                                 "n_judge": n_judge_g, "n_cells_over_thresh": n_over_g})

    # ---- Gather final state, check invariants across ALL rounds combined, write result ----
    my_final = [p for pid in local for p in local[pid]]
    all_final = comm.gather(my_final, root=0)
    all_initial = comm.gather(initial_flat_local, root=0)

    if rank == 0:
        final_flat = [p for sub in all_final for p in sub]
        initial_flat = [p for sub in all_initial for p in sub]

        final_weight = sum(p["weight"] for p in final_flat)
        ok_weight = abs(final_weight - initial_total_weight) < 1e-9

        ids = [p["id"] for p in final_flat]
        ok_unique = len(ids) == len(set(ids))

        # Per-particle lineage ledger audit -- strictly stronger than the global-sum check
        # above: verifies mass conservation at the level of individual particle ancestry, not
        # just the aggregate total, so it cannot be fooled by two independent bugs (one
        # duplicating mass, one silently dropping mass) that happen to cancel in the sum.
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
            lineage_weight = sum(orig_weight[i] for i in lineage)
            if abs(lineage_weight - p["weight"]) > 1e-9:
                ok_lineage_weight = False
        ok_lineage_complete = (seen_ids == all_orig_ids)
        ok_lineage = ok_lineage_disjoint and ok_lineage_weight and ok_lineage_complete

        oracle_pairs = ml.oracle_matching(initial_flat, args.thresh, args.cell_dx)

        # Two invariants, both universally valid (unlike "must reach zero" -- a cell whose
        # residents' true nearest neighbors ALL live in non-crowded cells is a legitimate,
        # by-design non-merge under the "both cells over threshold" rule, and this static,
        # no-particle-motion test harness gives it no way to ever resolve on its own; a real
        # timestep would drift particles between rounds and very likely break such a standoff):
        #  1. Never regress -- a cell that drops out of "over threshold" must never become over
        #     threshold again.
        #  2. Once merge activity stops (0 trivial + 0 judge for 2 rounds straight), the
        #     over-threshold count must stay perfectly flat -- confirms a genuine, stable
        #     standoff rather than some new oscillation/instability.
        over_seq = [rs["n_cells_over_thresh"] for rs in round_stats]
        ok_nonregress = all(b <= a for a, b in zip(over_seq, over_seq[1:]))
        ok_plateau_stable = True
        for i in range(2, len(round_stats)):
            quiet = (round_stats[i - 1]["n_trivial"] == 0 and round_stats[i - 1]["n_judge"] == 0
                     and round_stats[i - 2]["n_trivial"] == 0 and round_stats[i - 2]["n_judge"] == 0)
            if quiet and over_seq[i] != over_seq[i - 1]:
                ok_plateau_stable = False
        ok_converging = ok_nonregress and ok_plateau_stable

        result = {
            "scenario": args.scenario,
            "seed": args.seed,
            "shuffle_seed": args.shuffle_seed,
            "rounds": args.rounds,
            "n_initial": n_initial_global,
            "n_final": len(final_flat),
            "initial_weight": initial_total_weight,
            "final_weight": final_weight,
            "ok_weight_conserved": ok_weight,
            "ok_ids_unique": ok_unique,
            "ok_converging": ok_converging,
            "ok_lineage_disjoint": ok_lineage_disjoint,
            "ok_lineage_weight": ok_lineage_weight,
            "ok_lineage_complete": ok_lineage_complete,
            "n_oracle_merges": len(oracle_pairs),
            "round_stats": round_stats,
            "dim": len(final_flat[0]["pos"]) if final_flat else len(patches[0][0]),
            "final_particles": sorted(
                [(p["id"], *[round(c, 9) for c in p["pos"]], round(p["weight"], 9),
                  round(p["energy"], 9), p["rank"], p["patch"]) for p in final_flat]
            ),
        }
        with open(args.out, "w") as f:
            json.dump(result, f, indent=2)

        status = "PASS" if (ok_weight and ok_unique and ok_converging and ok_lineage) else "FAIL"
        totals = f"trivial={sum(rs['n_trivial'] for rs in round_stats)} judge={sum(rs['n_judge'] for rs in round_stats)}"
        over_trace = "->".join(str(rs["n_cells_over_thresh"]) for rs in round_stats)
        print(f"[{status}] scenario={args.scenario} seed={args.seed} shuffle={args.shuffle_seed} "
              f"rounds={args.rounds} n_initial={n_initial_global} n_final={len(final_flat)} "
              f"{totals} oracle={len(oracle_pairs)} over_thresh_trace={over_trace} "
              f"weight_ok={ok_weight} ids_unique={ok_unique} converging={ok_converging} "
              f"lineage_disjoint={ok_lineage_disjoint} lineage_weight={ok_lineage_weight} "
              f"lineage_complete={ok_lineage_complete}")


if __name__ == "__main__":
    main()
