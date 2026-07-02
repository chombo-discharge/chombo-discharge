"""
Core (MPI-agnostic) logic for the nearest-neighbor particle-merge prototype. Dimension-generic
(2D and 3D) -- all geometry operates on tuples of whatever length the scenario uses, matching how
the real C++ algorithm is templated on SpaceDim rather than hardcoded to 2D.

Data model
----------
A particle is a plain dict:
    {'pos': (x, y[, z]), 'weight': w, 'id': pid, 'energy': e, 'rank': owning_rank, 'patch': owning_patch_id}

'id' is a globally unique int64-ish integer. 'rank'/'patch' are the CURRENT owner -- always trusted
as the single source of truth (never inferred from message transport).

A patch is (lo, hi, owning_rank) -- lo/hi are position tuples of the scenario's dimensionality,
a half-open box [lo, hi). Every rank knows the full patch list (mirrors how a DisjointBoxLayout is
globally replicated in the real AMR code; only particle data is distributed).
"""

import math
import random


# ---------------------------------------------------------------------------
# Domain / scenario definitions
# ---------------------------------------------------------------------------

def scenario_corner4():
    """2x2 patches (2D), 4 ranks, one patch per rank -- all four patches meet at a single
    interior corner point (1,1). Every adjacent pair is cross-rank. The classic N-way corner
    stress case: a particle right at the corner can be ghosted to THREE different neighbors at
    once (the two edge-adjacent patches plus the diagonal one -- ghost_targets includes
    diagonal/corner neighbors, not just face neighbors)."""
    patches = [
        ((0, 0), (1, 1), 0),
        ((1, 0), (2, 1), 1),
        ((0, 1), (1, 2), 2),
        ((1, 1), (2, 2), 3),
    ]
    return patches, 4, [(1.0, 1.0)]


def scenario_mixed():
    """2x2 patches (2D), 2 ranks. Rank 0 owns the bottom row, rank 1 the top row. Two of the
    four edges are same-rank, two are cross-rank, both kinds converging at the same corner --
    exactly the case that broke the naive "same rank = safe" shortcut."""
    patches = [
        ((0, 0), (1, 1), 0),
        ((1, 0), (2, 1), 0),
        ((0, 1), (1, 2), 1),
        ((1, 1), (2, 2), 1),
    ]
    return patches, 2, [(1.0, 1.0)]


def scenario_grid3x3():
    """3x3 patches (9 total, 2D) over a 3x3 domain, 5 ranks with a deliberately irregular mix of
    single- and multi-patch ownership (including same-rank patches that are only DIAGONALLY
    adjacent, e.g. patch(0,0) and patch(2,1) both on rank 0, far apart). Four interior corners
    are all stressed AT ONCE in the same round -- the center patch touches all four of them
    simultaneously, so its rank is juggling four concurrent N-way judge relationships in a
    single pass. Checks that spatially-separated judge computations don't interfere."""
    def pid(i, j):
        return i * 3 + j

    ranks = {
        (0, 0): 0, (1, 0): 1, (2, 0): 2,
        (0, 1): 1, (1, 1): 4, (2, 1): 0,
        (0, 2): 2, (1, 2): 3, (2, 2): 4,
    }
    patches = [None] * 9
    for (i, j), r in ranks.items():
        patches[pid(i, j)] = ((float(i), float(j)), (float(i + 1), float(j + 1)), r)
    stress_points = [(1.0, 1.0), (2.0, 1.0), (1.0, 2.0), (2.0, 2.0)]
    return patches, 5, stress_points


def scenario_tjunction():
    """A genuine T-junction (3-way, not 4-way, 2D): patch A is DOUBLE-WIDE ([0,2] x [0,1], where
    a uniform grid would have had two separate unit patches), sitting below patch B ([0,1]x[1,2])
    and patch C ([1,2]x[1,2]). At (1,1), A meets it via the MIDDLE of its top edge (not a
    corner), while B and C each meet it at a corner -- exactly three patches converge, never
    four. Uneven patch sizes are common in real AMR domain decompositions; this exercises N=3
    arising from geometry, not from an accidental tie."""
    patches = [
        ((0.0, 0.0), (2.0, 1.0), 0),  # patch 0 (A), rank 0, double-wide
        ((0.0, 1.0), (1.0, 2.0), 1),  # patch 1 (B), rank 1
        ((1.0, 1.0), (2.0, 2.0), 2),  # patch 2 (C), rank 2
    ]
    return patches, 3, [(1.0, 1.0)]


def scenario_corner8():
    """2x2x2 patches (8 total, 3D), 8 ranks, one patch per rank -- all eight meet at a single
    interior corner point (1,1,1). The direct 3D analogue of corner4, but harder: a particle
    right at the corner can be ghosted to up to SEVEN different neighbors at once (3 face, 3
    edge, 1 opposite-corner), all cross-rank, all needing to route through and be reconciled by
    a single judge in one pass."""
    patches = []
    rank = 0
    coords = {}
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                coords[(i, j, k)] = rank
                patches.append(((float(i), float(j), float(k)), (float(i + 1), float(j + 1), float(k + 1)), rank))
                rank += 1
    return patches, 8, [(1.0, 1.0, 1.0)]


def scenario_mixed3d():
    """2x2x2 patches (3D), 4 ranks, each rank owning a diagonal PAIR of patches (opposite
    corners of the small cube) -- so every rank has one patch touching the central corner point
    directly and one far away on the opposite side, mixing same-rank non-adjacent ownership
    with the full 3D cross-rank corner conflict, analogous to `mixed` but one dimension up."""
    coords = [(i, j, k) for i in (0, 1) for j in (0, 1) for k in (0, 1)]
    # Pair each corner with its point-reflection through the cube center (0.5,0.5,0.5) -> same rank.
    rank_of = {}
    next_rank = 0
    for c in coords:
        if c in rank_of:
            continue
        opposite = tuple(1 - x for x in c)
        rank_of[c] = next_rank
        rank_of[opposite] = next_rank
        next_rank += 1
    patches = []
    for (i, j, k) in coords:
        r = rank_of[(i, j, k)]
        patches.append(((float(i), float(j), float(k)), (float(i + 1), float(j + 1), float(k + 1)), r))
    return patches, next_rank, [(1.0, 1.0, 1.0)]


SCENARIOS = {
    "corner4": scenario_corner4,
    "mixed": scenario_mixed,
    "grid3x3": scenario_grid3x3,
    "tjunction": scenario_tjunction,
    "corner8": scenario_corner8,
    "mixed3d": scenario_mixed3d,
}


# ---------------------------------------------------------------------------
# Geometry helpers (dimension-generic)
# ---------------------------------------------------------------------------

def point_in_box(pos, box, grow=0.0):
    lo, hi, _rank = box
    return all((lo[d] - grow) <= pos[d] < (hi[d] + grow) for d in range(len(pos)))


def locate_patch(pos, patches):
    """Which patch owns this position. Returns patch_id or None if outside the whole domain."""
    for pid, box in enumerate(patches):
        if point_in_box(pos, box, grow=0.0):
            return pid
    return None


def dist2(a, b):
    return sum((a[d] - b[d]) ** 2 for d in range(len(a)))


def cell_key(pos, dx=1.0):
    """Unclamped integer cell coordinate -- used for the per-cell crowding-threshold count.
    Deliberately NOT clamped to any one patch's box, so ghost particles bucket into their true
    origin cell rather than a boundary cell of the receiving patch."""
    return tuple(math.floor(pos[d] / dx) for d in range(len(pos)))


# ---------------------------------------------------------------------------
# Particle generation
# ---------------------------------------------------------------------------

_id_counter = [0]


def make_particle(pos, weight, rank, patch, rng, id_hint=None):
    if id_hint is not None:
        pid = id_hint
    else:
        # Poor-man's globally-unique id allocator for a prototype: rank-namespaced counter.
        # Mirrors what a real implementation might do without a central id service.
        _id_counter[0] += 1
        pid = rank * 1_000_000 + _id_counter[0]
    return {
        "pos": tuple(pos),
        "weight": weight,
        "id": pid,
        "energy": rng.uniform(0.0, 10.0),
        "rank": rank,
        "patch": patch,
        # Lineage: the set of ORIGINAL particle ids whose mass this particle's weight
        # ultimately derives from. For a freshly generated particle, that's just itself.
        # Propagated (unioned) through combine() -- a pure audit trail, never read by any
        # decision logic, so it can't influence the protocol, only verify it after the fact.
        "lineage": frozenset((pid,)),
    }


def gen_particles_for_patch(patch_id, box, rng, n_interior=6, n_boundary_cluster=0,
                             boundary_point=None, cluster_radius=0.08):
    """Generate a patch's own valid particles: some scattered through the interior, optionally a
    dense cluster right against a named boundary point (to force crowding exactly where the
    ghost/corner logic needs to be exercised)."""
    lo, hi, rank = box
    dim = len(lo)
    out = []
    for _ in range(n_interior):
        pos = tuple(rng.uniform(lo[d] + 0.05, hi[d] - 0.05) for d in range(dim))
        out.append(make_particle(pos, 1.0, rank, patch_id, rng))
    if n_boundary_cluster and boundary_point is not None:
        for _ in range(n_boundary_cluster):
            # Clip into this patch's own box so gen never accidentally seeds a particle that
            # geometrically belongs to a neighboring patch. The margin is itself randomized
            # (not a fixed 1e-6) so that clamping never collapses two distinct particles onto
            # the exact same floating-point position -- an exact-position tie is a meaningless
            # artifact of the generator, not something the merge protocol should be judged on.
            pos = []
            for d in range(dim):
                margin = rng.uniform(1e-6, 1e-3)
                v = boundary_point[d] + rng.uniform(-cluster_radius, cluster_radius)
                v = min(max(v, lo[d] + margin), hi[d] - margin)
                pos.append(v)
            out.append(make_particle(tuple(pos), 1.0, rank, patch_id, rng))
    return out


# ---------------------------------------------------------------------------
# Ghost targeting
# ---------------------------------------------------------------------------

def ghost_targets(pos, owner_patch_id, patches, ghost_width):
    """Which OTHER patches should receive this particle as a ghost. Includes diagonal
    (edge/corner, in 3D) neighbors, not just face neighbors -- deliberately, since that's
    exactly the topology that creates multi-way judging conflicts at a corner."""
    targets = []
    for pid, box in enumerate(patches):
        if pid == owner_patch_id:
            continue
        if point_in_box(pos, box, grow=ghost_width):
            targets.append(pid)
    return targets


def is_boundary_exposed(pos, owner_patch_id, patches, ghost_width):
    return len(ghost_targets(pos, owner_patch_id, patches, ghost_width)) > 0


# ---------------------------------------------------------------------------
# Nearest-neighbor candidate generation (brute force -- fine at this scale)
# ---------------------------------------------------------------------------

def nearest_neighbor(query_particle, pool, alive_ids):
    """pool: list of particles (locals + ghosts) query_particle may pair with (excludes itself).
    alive_ids: set of ids currently eligible. Returns (partner_particle, d2) or (None, None)."""
    best = None
    best_d2 = math.inf
    for p in pool:
        if p["id"] == query_particle["id"] or p["id"] not in alive_ids:
            continue
        d2 = dist2(query_particle["pos"], p["pos"])
        if d2 < best_d2:
            best_d2 = d2
            best = p
    return best, best_d2


# ---------------------------------------------------------------------------
# Merge arithmetic
# ---------------------------------------------------------------------------

def combine(p1, p2, new_id, new_rank, new_patch):
    w = p1["weight"] + p2["weight"]
    dim = len(p1["pos"])
    pos = tuple((p1["weight"] * p1["pos"][d] + p2["weight"] * p2["pos"][d]) / w for d in range(dim))
    e = (p1["weight"] * p1["energy"] + p2["weight"] * p2["energy"]) / w
    lineage = p1["lineage"] | p2["lineage"]
    return {"pos": pos, "weight": w, "id": new_id, "energy": e, "rank": new_rank, "patch": new_patch,
            "lineage": lineage}


# ---------------------------------------------------------------------------
# Oracle: global, boundary-free greedy nearest-pair matching (reference only)
# ---------------------------------------------------------------------------

def oracle_matching(all_particles, thresh, dx=1.0):
    """Same greedy nearest-pair-first algorithm, applied to the WHOLE domain as if it were one
    giant patch with no boundaries at all. Used only to sanity-check that every merge the
    distributed protocol actually performs is a genuine, locally-justified nearest pair -- NOT
    expected to match the distributed result exactly (the distributed protocol is deliberately
    conservative and can leave some boundary particles unmerged that this oracle would match)."""
    particles = {p["id"]: dict(p) for p in all_particles}
    alive = set(particles.keys())
    live_count = {}
    for p in particles.values():
        k = cell_key(p["pos"], dx)
        live_count[k] = live_count.get(k, 0) + 1

    edges = []
    for pid, p in particles.items():
        nbr, d2 = nearest_neighbor(p, list(particles.values()), alive)
        if nbr is not None:
            edges.append((d2, pid, nbr["id"]))
    edges.sort(key=lambda e: (e[0], e[1], e[2]))

    pairs = []
    for d2, a, b in edges:
        if a not in alive or b not in alive:
            continue
        ka = cell_key(particles[a]["pos"], dx)
        kb = cell_key(particles[b]["pos"], dx)
        if live_count.get(ka, 0) <= thresh or live_count.get(kb, 0) <= thresh:
            continue
        alive.discard(a)
        alive.discard(b)
        live_count[ka] = live_count.get(ka, 0) - 1
        live_count[kb] = live_count.get(kb, 0) - 1
        pairs.append((a, b, math.sqrt(d2)))
    return pairs
