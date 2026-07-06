"""
1D periodic mirror of run_scenario.py's run_one_round() -- identical algorithm/protocol
(RULES.md), only the geometry helpers differ (periodic wraparound instead of a bounded 2D box).
One patch per rank, rank index == patch index, patches tile [0, n_patches*patch_size) with patch i
spanning [i*patch_size, (i+1)*patch_size), wrapping at the domain length.
"""

import merge_lib_1d as ml


def run_one_round_1d(comm, rank, size, n_patches, patch_size, local, thresh, ghost_width, cell_dx,
                      shuffle_rng, id_counter):
    pid = rank  # one patch per rank, 1:1

    # ---- Ghost-fill (alltoall #1), fresh every round, before any decision ----
    outgoing_ghosts = [[] for _ in range(size)]
    for p in local:
        for dest_patch, shipped_pos in ml.ghost_targets(p["pos"], pid, ghost_width, patch_size, n_patches):
            g = dict(p)
            g["pos"] = shipped_pos
            outgoing_ghosts[dest_patch].append(g)  # dest_patch == dest_rank, 1:1
    incoming_ghosts = comm.alltoall(outgoing_ghosts)
    ghosts = [g for batch in incoming_ghosts for g in batch]

    # ---- Trivial tier + candidate generation (single patch per rank, so no per-patch loop) ----
    valid = local
    alive_valid = {p["id"] for p in valid}
    exposed = {p["id"]: ml.is_boundary_exposed(p["pos"], pid, ghost_width, patch_size, n_patches)
               for p in valid}

    lc = {}
    for p in valid:
        k = ml.cell_key(p["pos"], cell_dx, patch_size * n_patches)
        lc[k] = lc.get(k, 0) + 1
    for g in ghosts:
        k = ml.cell_key(g["pos"], cell_dx, patch_size * n_patches)
        lc[k] = lc.get(k, 0) + 1

    ghost_ids = {g["id"] for g in ghosts}
    pool = valid + ghosts

    order = list(range(len(valid)))
    shuffle_rng.shuffle(order)

    has_outgoing = {}
    outgoing_target = {}
    trivial_merges = []
    outgoing_proposals = [[] for _ in range(size)]
    L = patch_size * n_patches

    edges = []
    for idx in order:
        p = valid[idx]
        if lc.get(ml.cell_key(p["pos"], cell_dx, L), 0) <= thresh:
            continue
        nbr, d2 = ml.nearest_neighbor(p, pool, alive_valid | ghost_ids)
        if nbr is None:
            continue
        edges.append((d2, p, nbr))

    edges.sort(key=lambda e: (e[0], e[1]["id"], e[2]["id"]))

    for d2, p, nbr in edges:
        if p["id"] not in alive_valid:
            continue
        if nbr["patch"] == pid:
            if nbr["id"] not in alive_valid:
                continue
            if has_outgoing.get(nbr["id"], False):
                continue
            trivial_ok = (not exposed[p["id"]]) and (not exposed.get(nbr["id"], False))
        else:
            trivial_ok = False
        ka = ml.cell_key(p["pos"], cell_dx, L)
        kb = ml.cell_key(nbr["pos"], cell_dx, L)
        if lc.get(ka, 0) <= thresh or lc.get(kb, 0) <= thresh:
            continue
        if trivial_ok:
            alive_valid.discard(p["id"])
            alive_valid.discard(nbr["id"])
            lc[ka] = lc.get(ka, 0) - 1
            lc[kb] = lc.get(kb, 0) - 1
            trivial_merges.append((p, nbr))
        else:
            has_outgoing[p["id"]] = True
            outgoing_target[p["id"]] = nbr["id"]
            dest_rank = nbr["rank"]
            outgoing_proposals[dest_rank].append({"source": p, "target_id": nbr["id"], "d2": d2})

    # ---- Proposals routed by target ownership (alltoall #2) ----
    incoming_proposals = comm.alltoall(outgoing_proposals)
    flat_incoming = [prop for bucket in incoming_proposals for prop in bucket]
    shuffle_rng.shuffle(flat_incoming)

    # ---- Judge processing ----
    by_target = {}
    for prop in flat_incoming:
        by_target.setdefault(prop["target_id"], []).append(prop)

    judge_merges = []
    verdicts_out = [[] for _ in range(size)]

    target_lookup = {p["id"]: p for p in valid if p["id"] in alive_valid}

    candidates = []
    deferred_props = {}
    for target_id, props in by_target.items():
        if target_id not in target_lookup:
            for prop in props:
                verdicts_out[prop["source"]["rank"]].append({"proposer_id": prop["source"]["id"], "accepted": False})
            continue
        tparticle = target_lookup[target_id]
        if has_outgoing.get(target_id, False):
            my_target = outgoing_target.get(target_id)
            mutual = [pr for pr in props if pr["source"]["id"] == my_target] if my_target is not None else []
            if mutual and target_id < my_target:
                candidates.append((mutual[0]["d2"], target_id, mutual[0], tparticle))
                deferred_props[target_id] = props
            else:
                for prop in props:
                    verdicts_out[prop["source"]["rank"]].append({"proposer_id": prop["source"]["id"], "accepted": False})
            continue
        winner = min(props, key=lambda pr: (pr["d2"], pr["source"]["id"]))
        candidates.append((winner["d2"], target_id, winner, tparticle))
        deferred_props[target_id] = props

    candidates.sort(key=lambda c: (c[0], c[1]))
    accepted_targets = {}
    for d2, target_id, winner, tparticle in candidates:
        if target_id not in alive_valid:
            continue
        ka = ml.cell_key(tparticle["pos"], cell_dx, L)
        if lc.get(ka, 0) <= thresh:
            continue
        alive_valid.discard(target_id)
        lc[ka] = lc.get(ka, 0) - 1
        judge_merges.append((tparticle, winner["source"]))
        accepted_targets[target_id] = winner

    for target_id, props in deferred_props.items():
        winner = accepted_targets.get(target_id)
        for prop in props:
            won = winner is not None and prop is winner
            verdicts_out[prop["source"]["rank"]].append({"proposer_id": prop["source"]["id"], "accepted": won})

    # ---- Verdicts routed back to proposers (alltoall #3) ----
    incoming_verdicts = comm.alltoall(verdicts_out)
    flat_verdicts = [v for bucket in incoming_verdicts for v in bucket]
    accepted_source_ids = {v["proposer_id"] for v in flat_verdicts if v["accepted"]}
    alive_valid -= accepted_source_ids

    # ---- Assemble final valid list, place newly-created merged particles ----
    def fresh_id():
        id_counter[0] += 1
        return id_counter[0]

    final_local = [p for p in valid if p["id"] in alive_valid]

    pending_placement = []
    for a, b in trivial_merges:
        pending_placement.append(ml.combine(a, b, fresh_id(), rank, None))
    for target_p, winner_p in judge_merges:
        pending_placement.append(ml.combine(target_p, winner_p, fresh_id(), rank, None))

    outgoing_placement = [[] for _ in range(size)]
    for m in pending_placement:
        dest_patch = ml.locate_patch(m["pos"], patch_size, n_patches)
        m["patch"] = dest_patch
        m["rank"] = dest_patch
        m["pos"] = ml.wrap(m["pos"], L)
        if dest_patch == rank:
            final_local.append(m)
        else:
            outgoing_placement[dest_patch].append(m)

    incoming_placement = comm.alltoall(outgoing_placement)
    for items in incoming_placement:
        for m in items:
            final_local.append(m)

    return final_local, len(trivial_merges), len(judge_merges)
