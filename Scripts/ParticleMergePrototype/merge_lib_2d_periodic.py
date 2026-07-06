"""
2D, periodic-boundary pure-logic module -- the corner-preserving fix to the 1D reduction (which
turned out to structurally lack corners, and showed no patch-frequency spike as a result). Same
algorithm (RULES.md), periodic in BOTH x and y so there's no domain-edge confound, on an
n_patches_x x n_patches_y grid of unit patches, one rank per patch.

ghost_targets checks all 8 neighbor offsets (face + diagonal), matching the original 2D
merge_lib.py's "check every other patch" behavior -- for a 3x3 periodic grid this enumerates
exactly the other 8 patches. A ghost shipped across either periodic seam (or both, for a diagonal
wrap) gets its position shifted by the domain length in that axis, so the receiving patch sees it
at its true local geometric proximity.
"""

import math


def wrap(x, L):
    return x % L


def dist2(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


def cell_key(pos, cell_dx, Lx, Ly):
    return (int(math.floor(wrap(pos[0], Lx) / cell_dx)), int(math.floor(wrap(pos[1], Ly) / cell_dx)))


def cell_chebyshev_distance(key_a, key_b, n_cells_x, n_cells_y):
    """Distance between two cell-index keys (as returned by cell_key), in units of whole cells,
    using Chebyshev (max of per-axis) distance -- same cell = 0, any of the 8 Moore-neighborhood
    adjacent cells (including diagonals) = 1, etc. Periodic-wraparound-aware: the shortest path
    across either periodic seam counts, matching how ghost_targets already treats the domain."""
    dx = abs(key_a[0] - key_b[0])
    dx = min(dx, n_cells_x - dx)
    dy = abs(key_a[1] - key_b[1])
    dy = min(dy, n_cells_y - dy)
    return max(dx, dy)


def locate_patch(pos, patch_size, n_patches_x, n_patches_y):
    Lx, Ly = patch_size * n_patches_x, patch_size * n_patches_y
    i = int(wrap(pos[0], Lx) / patch_size) % n_patches_x
    j = int(wrap(pos[1], Ly) / patch_size) % n_patches_y
    return j * n_patches_x + i


def ghost_targets(pos, i, j, ghost_width, patch_size, n_patches_x, n_patches_y):
    """Returns [(dest_pid, shipped_pos), ...] over all 8 neighbor offsets."""
    Lx, Ly = patch_size * n_patches_x, patch_size * n_patches_y
    out = []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            ni, nj = i + di, j + dj
            lo_x, lo_y = ni * patch_size, nj * patch_size
            hi_x, hi_y = lo_x + patch_size, lo_y + patch_size
            if (lo_x - ghost_width) <= pos[0] < (hi_x + ghost_width) and \
               (lo_y - ghost_width) <= pos[1] < (hi_y + ghost_width):
                dest_i = ni % n_patches_x
                dest_j = nj % n_patches_y
                shift_x = -Lx if ni >= n_patches_x else (Lx if ni < 0 else 0.0)
                shift_y = -Ly if nj >= n_patches_y else (Ly if nj < 0 else 0.0)
                dest_pid = dest_j * n_patches_x + dest_i
                out.append((dest_pid, (pos[0] + shift_x, pos[1] + shift_y)))
    return out


def is_boundary_exposed(pos, i, j, ghost_width, patch_size, n_patches_x, n_patches_y):
    return len(ghost_targets(pos, i, j, ghost_width, patch_size, n_patches_x, n_patches_y)) > 0


def nearest_neighbor(query_particle, pool, alive_ids, exclude_ids=None):
    """exclude_ids: candidates to skip even if alive -- used for the fallback-candidate retry,
    where a particle's first choice turned out to be blocked and it needs its NEXT nearest,
    excluding candidates already tried and found unusable this pass."""
    best = None
    best_d2 = math.inf
    qpos = query_particle["pos"]
    qid = query_particle["id"]
    for p in pool:
        if p["id"] == qid or p["id"] not in alive_ids:
            continue
        if exclude_ids is not None and p["id"] in exclude_ids:
            continue
        d2 = dist2(qpos, p["pos"])
        if d2 < best_d2:
            best_d2 = d2
            best = p
    return best, best_d2


def combine(p1, p2, new_id, new_rank, new_patch):
    w = p1["weight"] + p2["weight"]
    pos = ((p1["pos"][0] * p1["weight"] + p2["pos"][0] * p2["weight"]) / w,
           (p1["pos"][1] * p1["weight"] + p2["pos"][1] * p2["weight"]) / w)
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


def gen_particles_uniform(pid, i, j, patch_size, rng, n, id_counter):
    lo_x, lo_y = i * patch_size, j * patch_size
    out = []
    for _ in range(n):
        id_counter[0] += 1
        pos = (rng.uniform(lo_x + 1e-9, lo_x + patch_size - 1e-9),
               rng.uniform(lo_y + 1e-9, lo_y + patch_size - 1e-9))
        out.append(make_particle(pos, 1.0, pid, pid, rng, id_counter[0]))
    return out
