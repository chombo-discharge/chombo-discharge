"""
1D, periodic-boundary pure-logic module -- a faster testbed for iterating on fixes to the
patch-frequency aliasing found in the 2D prototype. Same algorithm (RULES.md), same validated
propose/judge/verdict protocol; the only thing that changes is the geometry: 1D, periodic, so
there is no true domain edge to confound the patch-periodicity measurement -- every patch boundary
is an ordinary interior boundary, and there are exactly n_patches of them, evenly spaced, which is
exactly the structure we want a clean spectrum against.

Periodic ghost-shipping: a ghost shipped across the wraparound seam (patch 0's left neighbor is
the last patch, and vice versa) is given a position shifted by +/- the domain length L, so the
RECEIVING patch sees it at its true local geometric proximity rather than its literal far-away
wrapped coordinate. cell_key() wraps internally, so this shifted position still maps back to the
ghost's correct true origin cell -- the shift only matters for the local distance computation.
"""

import math


def wrap(x, L):
    return x % L


def dist2(a, b):
    return (a - b) ** 2


def cell_key(pos, cell_dx, L):
    return int(math.floor(wrap(pos, L) / cell_dx))


def locate_patch(pos, patch_size, n_patches):
    L = patch_size * n_patches
    return int(wrap(pos, L) / patch_size) % n_patches


def ghost_targets(pos, owner_patch, ghost_width, patch_size, n_patches):
    """Returns [(dest_patch, shipped_pos), ...]."""
    L = patch_size * n_patches
    lo = owner_patch * patch_size
    hi = lo + patch_size
    out = []
    if (pos - lo) < ghost_width:
        dest = (owner_patch - 1) % n_patches
        shift = L if owner_patch == 0 else 0.0
        out.append((dest, pos + shift))
    if (hi - pos) < ghost_width:
        dest = (owner_patch + 1) % n_patches
        shift = -L if owner_patch == n_patches - 1 else 0.0
        out.append((dest, pos + shift))
    return out


def is_boundary_exposed(pos, owner_patch, ghost_width, patch_size, n_patches):
    return len(ghost_targets(pos, owner_patch, ghost_width, patch_size, n_patches)) > 0


def nearest_neighbor(query_particle, pool, alive_ids):
    """pool: list of particle dicts (locals use their own canonical 'pos'; ghosts carry their
    shifted 'pos', already adjusted at ghost-fill time). Returns (partner, d2) or (None, None)."""
    best = None
    best_d2 = math.inf
    qpos = query_particle["pos"]
    qid = query_particle["id"]
    for p in pool:
        if p["id"] == qid or p["id"] not in alive_ids:
            continue
        d2 = dist2(qpos, p["pos"])
        if d2 < best_d2:
            best_d2 = d2
            best = p
    return best, best_d2


def combine(p1, p2, new_id, new_rank, new_patch):
    w = p1["weight"] + p2["weight"]
    pos = (p1["pos"] * p1["weight"] + p2["pos"] * p2["weight"]) / w
    energy = (p1["weight"] * p1["energy"] + p2["weight"] * p2["energy"]) / w
    return {
        "pos": pos, "weight": w, "id": new_id, "energy": energy,
        "rank": new_rank, "patch": new_patch,
        "lineage": p1["lineage"] | p2["lineage"],
    }


def make_particle(pos, weight, rank, patch, rng, id_hint):
    return {
        "pos": pos, "weight": weight, "id": id_hint,
        "energy": rng.uniform(0.0, 10.0),
        "rank": rank, "patch": patch,
        "lineage": frozenset((id_hint,)),
    }


def gen_particles_uniform(patch_id, patch_size, rng, n, id_counter):
    lo = patch_id * patch_size
    out = []
    for _ in range(n):
        id_counter[0] += 1
        pos = rng.uniform(lo + 1e-9, lo + patch_size - 1e-9)
        out.append(make_particle(pos, 1.0, patch_id, patch_id, rng, id_counter[0]))
    return out
