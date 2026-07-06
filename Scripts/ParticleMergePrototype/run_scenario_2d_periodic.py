"""
2D periodic mirror of run_scenario.py's run_one_round() -- identical algorithm/protocol
(RULES.md); geometry differs (8-neighbor periodic wraparound), and patches are assigned to ranks
many-to-one: n_patches_y rows are split into n_ranks contiguous blocks (each rank owning one or
more whole rows), so the grid can be scaled arbitrarily (toward the real target of ~10k patches
per direction) while the actual MPI rank count stays fixed at whatever the machine actually has
(e.g. 12) -- not tied 1:1 to patch or even row count.
"""

import merge_lib_2d_periodic as ml


def patch_to_rank(pid, n_patches_x, n_patches_y, n_ranks):
    row = pid // n_patches_x
    return (row * n_ranks) // n_patches_y


def run_one_round_2d(comm, rank, size, n_patches_x, n_patches_y, patch_size, my_patches, local,
                      thresh, ghost_width, cell_dx, shuffle_rng, id_counter,
                      iterate_local_to_convergence=False, max_fallback_candidates=0,
                      max_cell_distance=None):
    """iterate_local_to_convergence: see CD_NearestNeighborParticleMerge.H's
    a_iterateLocalTierToConvergence and RULES.md, "Optional: iterate the local (trivial) tier to
    convergence within a round". When False (default): candidate search + trivial-tier resolve
    runs exactly once, matching the originally documented single-pass behavior. When True: that
    step repeats -- reusing the SAME ghost-fill and SAME exposed flags (neither depends on which
    local particles are still alive) -- until a full local iteration commits zero further merges.
    The cross-patch propose/judge/verdict exchange below is unaffected either way -- exactly once,
    always.

    max_fallback_candidates: when a particle's chosen candidate turns out to be blocked (already
    consumed, or already has its own outgoing commitment) -- the single largest loss mechanism
    measured empirically -- immediately re-query for the next-nearest still-available, non-
    excluded candidate, up to this many extra tries, INLINE at the point of failure (not deferred
    to preserve strict global distance ordering across every particle's retries -- a deliberate,
    pragmatic simplification: retries only ever fire for edges whose first choice already failed,
    so they're not competing on equal footing with anyone's first attempt anyway). 0 (default)
    disables this and matches the original single-candidate-per-particle behavior. The proposal
    routed to generateProposals()/judgeProposals() for a particle that exhausts every fallback
    still uses its TRUE (rank-1) nearest neighbor, matching existing proposal semantics exactly --
    only local trivial-tier eligibility is affected by the fallback search.

    max_cell_distance: physical merge-eligibility cap, in whole grid cells (Chebyshev distance
    between cell indices -- same cell = 0, any Moore-neighborhood-adjacent cell including
    diagonals = 1). None (default) disables this and matches the original unrestricted-search
    behavior (any candidate, however far, is eligible). A particle whose true nearest neighbor
    exceeds this cap is NOT proposed to (unlike the busy/exposed/ghost cases -- this is a
    physical-validity rule, not a coordination one; a distant pair is exactly as invalid via
    propose/judge as it would be trivially) -- if a_maxFallbackCandidates > 0 it gets a fallback
    retry for a CLOSER candidate instead, otherwise it is simply left unmatched this round, with
    no side effects on anyone else's processing (same treatment as a stale candidate)."""
    Lx, Ly = patch_size * n_patches_x, patch_size * n_patches_y
    n_cells_x = round(Lx / cell_dx)
    n_cells_y = round(Ly / cell_dx)

    # ---- Ghost-fill (alltoall #1), fresh once, reused across every local iteration below ----
    outgoing_ghosts = [[] for _ in range(size)]
    for pid in my_patches:
        i, j = pid % n_patches_x, pid // n_patches_x
        for p in local[pid]:
            for dest_pid, shipped_pos in ml.ghost_targets(p["pos"], i, j, ghost_width, patch_size,
                                                            n_patches_x, n_patches_y):
                dest_rank = patch_to_rank(dest_pid, n_patches_x, n_patches_y, size)
                g = dict(p)
                g["pos"] = shipped_pos
                outgoing_ghosts[dest_rank].append((g, dest_pid))
    incoming_ghosts = comm.alltoall(outgoing_ghosts)

    ghosts = {pid: [] for pid in my_patches}
    for items in incoming_ghosts:
        for g, dest_pid in items:
            ghosts[dest_pid].append(g)

    # ---- Per-patch state that's fixed for the whole round, computed once regardless of how
    # many local iterations follow: which ids are alive, exposure (a geometric property of a
    # cell, never changes), and each patch's own valid-particle-list/ghost-pool/ghost-ids. ----
    alive_valid = {}
    exposed = {}
    live_count = {}
    valid_by_patch = {}
    pool_by_patch = {}
    ghost_ids_by_patch = {}

    for pid in my_patches:
        i, j = pid % n_patches_x, pid // n_patches_x
        valid = local[pid]
        valid_by_patch[pid] = valid
        alive_valid[pid] = {p["id"] for p in valid}
        exposed.update({p["id"]: ml.is_boundary_exposed(p["pos"], i, j, ghost_width, patch_size,
                                                          n_patches_x, n_patches_y)
                         for p in valid})

        lc = {}
        for p in valid:
            k = ml.cell_key(p["pos"], cell_dx, Lx, Ly)
            lc[k] = lc.get(k, 0) + 1
        for g in ghosts[pid]:
            k = ml.cell_key(g["pos"], cell_dx, Lx, Ly)
            lc[k] = lc.get(k, 0) + 1
        live_count[pid] = lc

        ghost_ids_by_patch[pid] = {g["id"] for g in ghosts[pid]}
        pool_by_patch[pid] = valid + ghosts[pid]

    # ---- Local (trivial-tier) search/resolve -- one pass, or looped to convergence ----
    trivial_merges = []
    has_outgoing = {}
    outgoing_target = {}
    outgoing_proposals = [[] for _ in range(size)]

    while True:
        # Reset per-iteration bookkeeping: a particle's fate this iteration is decided fresh,
        # not carried over from an earlier iteration (its situation may have changed).
        has_outgoing = {}
        outgoing_target = {}
        outgoing_proposals = [[] for _ in range(size)]
        n_committed_this_iteration = 0

        patch_order = my_patches[:]
        shuffle_rng.shuffle(patch_order)

        for pid in patch_order:
            valid = valid_by_patch[pid]
            lc = live_count[pid]
            pool = pool_by_patch[pid]
            ghost_ids = ghost_ids_by_patch[pid]
            exposed_here = exposed

            order = list(range(len(valid)))
            shuffle_rng.shuffle(order)

            edges = []
            for idx in order:
                p = valid[idx]
                if p["id"] not in alive_valid[pid]:
                    continue
                if lc.get(ml.cell_key(p["pos"], cell_dx, Lx, Ly), 0) <= thresh:
                    continue
                nbr, d2 = ml.nearest_neighbor(p, pool, alive_valid[pid] | ghost_ids)
                if nbr is None:
                    continue
                edges.append((d2, p, nbr))

            edges.sort(key=lambda e: (e[0], e[1]["id"], e[2]["id"]))

            for d2, p, nbr in edges:
                if p["id"] not in alive_valid[pid]:
                    continue

                # Try the original (true nearest) candidate first, then -- if it's blocked (not
                # if it's genuinely ineligible due to crowding, and not if p itself is exposed,
                # since neither of those can be fixed by trying a different candidate -- up to
                # max_fallback_candidates progressively-farther still-available candidates before
                # giving up and proposing with the ORIGINAL true nearest neighbor, unchanged.
                cand, cd2 = nbr, d2
                tried_bad = set()
                committed = False
                attempts = 0
                p_exposed = exposed_here[p["id"]]
                p_cell = ml.cell_key(p["pos"], cell_dx, Lx, Ly)

                def _too_far(other):
                    return max_cell_distance is not None and ml.cell_chebyshev_distance(
                        p_cell, ml.cell_key(other["pos"], cell_dx, Lx, Ly), n_cells_x, n_cells_y
                    ) > max_cell_distance

                # Whether the ORIGINAL (un-retried) candidate was already gone, OR too far away,
                # by the time this edge was processed -- these are the cases that must NOT
                # trigger a proposal (see resolveTrivialTier()'s docs: ghost / exposed /
                # candidate-already-busy all trigger a proposal, but a candidate that no longer
                # exists, or one that's too physically far away to be a valid merge regardless of
                # coordination, does not -- there is nothing meaningful/valid to propose,
                # and marking p "has an outgoing commitment" here would needlessly block OTHER
                # particles from using p, which hasn't actually committed to anything).
                orig_candidate_stale = (nbr["patch"] == pid) and (nbr["id"] not in alive_valid[pid])
                orig_too_far = _too_far(nbr)
                while True:
                    if cand["patch"] == pid:
                        cand_alive = cand["id"] in alive_valid[pid]
                        cand_busy = cand_alive and has_outgoing.get(cand["id"], False)
                        trivial_ok = cand_alive and (not cand_busy) and (not p_exposed) and \
                            (not exposed.get(cand["id"], False)) and (not _too_far(cand))
                    else:
                        trivial_ok = False

                    if trivial_ok:
                        ka = ml.cell_key(p["pos"], cell_dx, Lx, Ly)
                        kb = ml.cell_key(cand["pos"], cell_dx, Lx, Ly)
                        if lc.get(ka, 0) <= thresh or lc.get(kb, 0) <= thresh:
                            break  # genuine crowding limit -- not fixable by a different candidate
                        alive_valid[pid].discard(p["id"])
                        alive_valid[pid].discard(cand["id"])
                        lc[ka] = lc.get(ka, 0) - 1
                        lc[kb] = lc.get(kb, 0) - 1
                        trivial_merges.append((p, cand, pid))
                        n_committed_this_iteration += 1
                        committed = True
                        break

                    if p_exposed or attempts >= max_fallback_candidates:
                        break
                    tried_bad.add(cand["id"])
                    attempts += 1
                    new_cand, new_d2 = ml.nearest_neighbor(p, pool, alive_valid[pid] | ghost_ids,
                                                           exclude_ids=tried_bad)
                    if new_cand is None:
                        break
                    cand, cd2 = new_cand, new_d2

                if not committed and not orig_candidate_stale and not orig_too_far:
                    # Proposal semantics unchanged: always the ORIGINAL true nearest neighbor,
                    # never a fallback candidate -- fallbacks only affect local trivial-tier
                    # eligibility, never what gets proposed cross-patch.
                    has_outgoing[p["id"]] = True
                    outgoing_target[p["id"]] = nbr["id"]
                    dest_rank = patch_to_rank(nbr["patch"], n_patches_x, n_patches_y, size)
                    outgoing_proposals[dest_rank].append({"source": p, "target_id": nbr["id"], "d2": d2})
                # else (orig_candidate_stale, or orig_too_far, not committed): this edge is simply
                # invalid or physically ineligible -- p is left exactly as it was, fully available
                # to others, deferred to next round's fresh candidate search. No proposal, no
                # has_outgoing marking.

        if not iterate_local_to_convergence or n_committed_this_iteration == 0:
            break

    # ---- Proposals routed by target ownership (alltoall #2) ----
    incoming_proposals = comm.alltoall(outgoing_proposals)
    flat_incoming = [prop for bucket in incoming_proposals for prop in bucket]
    shuffle_rng.shuffle(flat_incoming)

    # ---- Judge processing (local, per rank, across all of this rank's own patches) ----
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

    candidates = []
    deferred_props = {}
    for target_id, props in by_target.items():
        if target_id not in target_lookup:
            for prop in props:
                verdicts_out[patch_to_rank(prop["source"]["patch"], n_patches_x, n_patches_y, size)].append(
                    {"proposer_id": prop["source"]["id"], "accepted": False})
            continue
        tpid, tparticle = target_lookup[target_id]
        if has_outgoing.get(target_id, False):
            my_target = outgoing_target.get(target_id)
            mutual = [pr for pr in props if pr["source"]["id"] == my_target] if my_target is not None else []
            if mutual and target_id < my_target:
                candidates.append((mutual[0]["d2"], target_id, mutual[0], tpid, tparticle))
                deferred_props[target_id] = props
            else:
                for prop in props:
                    verdicts_out[patch_to_rank(prop["source"]["patch"], n_patches_x, n_patches_y, size)].append(
                        {"proposer_id": prop["source"]["id"], "accepted": False})
            continue
        winner = min(props, key=lambda pr: (pr["d2"], pr["source"]["id"]))
        candidates.append((winner["d2"], target_id, winner, tpid, tparticle))
        deferred_props[target_id] = props

    candidates.sort(key=lambda c: (c[0], c[1]))
    accepted_targets = {}
    for d2, target_id, winner, tpid, tparticle in candidates:
        if target_id not in alive_valid[tpid]:
            continue
        ka = ml.cell_key(tparticle["pos"], cell_dx, Lx, Ly)
        if live_count[tpid].get(ka, 0) <= thresh:
            continue
        alive_valid[tpid].discard(target_id)
        live_count[tpid][ka] = live_count[tpid].get(ka, 0) - 1
        judge_merges.append((tparticle, winner["source"], tpid))
        accepted_targets[target_id] = winner

    for target_id, props in deferred_props.items():
        winner = accepted_targets.get(target_id)
        for prop in props:
            won = winner is not None and prop is winner
            verdicts_out[patch_to_rank(prop["source"]["patch"], n_patches_x, n_patches_y, size)].append(
                {"proposer_id": prop["source"]["id"], "accepted": won})

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
        dest_pid = ml.locate_patch(m["pos"], patch_size, n_patches_x, n_patches_y)
        dest_rank = patch_to_rank(dest_pid, n_patches_x, n_patches_y, size)
        m["patch"] = dest_pid
        m["rank"] = dest_rank
        m["pos"] = (ml.wrap(m["pos"][0], Lx), ml.wrap(m["pos"][1], Ly))
        if dest_rank == rank:
            final_local.setdefault(dest_pid, [])
            final_local[dest_pid].append(m)
        else:
            outgoing_placement[dest_rank].append(m)

    incoming_placement = comm.alltoall(outgoing_placement)
    for items in incoming_placement:
        for m in items:
            final_local.setdefault(m["patch"], [])
            final_local[m["patch"]].append(m)

    return final_local, len(trivial_merges), len(judge_merges)
