"""
Outer (non-MPI) harness: launches many `mpirun -np K run_scenario.py ...` combinations and
cross-checks the results across all six scenarios (2D and 3D, corners/mixed-ownership/T-junction),
each run for multiple rounds (consecutive simulation timesteps). Checks beyond each run's own
weight/uniqueness/convergence invariants:

  Determinism / race-freedom: for a FIXED (scenario, seed) -- i.e. a fixed physical particle
  configuration -- running with many different --shuffle-seed values (which only perturb local
  processing order, never the physics) must produce the IDENTICAL canonicalized final state
  every time, after every round. Freshly-created particle ids are arbitrary labels (their numeric
  value depends on the order merges happened to be created in), so they're stripped before
  comparing -- what has to match is the multiset of (position, weight, energy) for newly created
  particles, and the exact id set for particles that simply survived untouched.
"""

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VENV_PY = os.path.join(HERE, "venv", "bin", "python")
RESULTS = os.path.join(HERE, "results")

SCENARIO_RANKS = {
    "corner4": 4,
    "mixed": 2,
    "grid3x3": 5,
    "tjunction": 3,
    "corner8": 8,
    "mixed3d": 4,
}


def run_one(scenario, seed, shuffle_seed, thresh, ghost_width, cell_dx, rounds, nranks):
    out_path = os.path.join(RESULTS, f"{scenario}_s{seed}_sh{shuffle_seed}.json")
    cmd = [
        "mpirun", "--oversubscribe", "-np", str(nranks), VENV_PY, "run_scenario.py",
        "--scenario", scenario, "--seed", str(seed), "--shuffle-seed", str(shuffle_seed),
        "--thresh", str(thresh), "--ghost-width", str(ghost_width), "--cell-dx", str(cell_dx),
        "--rounds", str(rounds), "--out", out_path,
    ]
    proc = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"[CRASH] {scenario} seed={seed} shuffle={shuffle_seed}")
        print(proc.stdout)
        print(proc.stderr)
        return None
    print(proc.stdout.strip())
    with open(out_path) as f:
        return json.load(f)


def canonicalize(final_particles, dim, n_original_max_id=99_999_999):
    """Split into (frozen set of surviving-original ids+data, multiset of new-particle physical
    tuples) so runs that assign different arbitrary fresh ids can still be compared for physical
    equivalence. Row shape: (id, x, y[, z], weight, energy, rank, patch)."""
    survivors = []
    created = []
    for row in final_particles:
        pid = row[0]
        rest = row[1:1 + dim + 2]  # position (dim) + weight + energy
        if pid <= n_original_max_id:
            survivors.append((pid,) + tuple(rest))
        else:
            created.append(tuple(rest))
    return frozenset(survivors), tuple(sorted(created))


def check_determinism(results):
    by_key = {}
    for r in results:
        if r is None:
            continue
        key = (r["scenario"], r["seed"])
        by_key.setdefault(key, []).append(r)

    all_ok = True
    for (scenario, seed), runs in by_key.items():
        dim = runs[0]["dim"]
        canon = [canonicalize(r["final_particles"], dim) for r in runs]
        first = canon[0]
        mismatches = [i for i, c in enumerate(canon) if c != first]
        if mismatches:
            all_ok = False
            print(f"[FAIL] determinism: scenario={scenario} seed={seed} -- "
                  f"{len(mismatches)}/{len(runs)} shuffle-seeds disagree with shuffle_seed={runs[0]['shuffle_seed']}")
            for i in mismatches[:2]:
                print(f"    shuffle_seed={runs[i]['shuffle_seed']} differs")
        else:
            print(f"[PASS] determinism: scenario={scenario} seed={seed} -- "
                  f"{len(runs)} shuffle-seeds agree (n_final={runs[0]['n_final']})")
    return all_ok


def main():
    os.makedirs(RESULTS, exist_ok=True)

    seeds = list(range(1, 11))
    shuffle_seeds = list(range(6))
    thresh = 2
    ghost_width = 0.15
    cell_dx = 0.2
    rounds = 6

    all_results = []
    any_run_failed = False

    for scenario, nranks in SCENARIO_RANKS.items():
        for seed in seeds:
            for shuffle_seed in shuffle_seeds:
                r = run_one(scenario, seed, shuffle_seed, thresh, ghost_width, cell_dx, rounds, nranks)
                all_results.append(r)
                if r is None or not (r["ok_weight_conserved"] and r["ok_ids_unique"] and r["ok_converging"]
                                      and r["ok_lineage_disjoint"] and r["ok_lineage_weight"] and r["ok_lineage_complete"]):
                    any_run_failed = True

    print()
    det_ok = check_determinism(all_results)

    print()
    print("Pressure test: asymmetric ghost visibility (proposal source never a pre-existing ghost)")
    pressure_proc = subprocess.run(
        ["mpirun", "--oversubscribe", "-np", "2", VENV_PY, "pressure_test_asymmetric_visibility.py"],
        cwd=HERE, capture_output=True, text=True,
    )
    print(pressure_proc.stdout.strip())
    pressure_ok = pressure_proc.returncode == 0
    if not pressure_ok:
        print(pressure_proc.stderr)

    print()
    if any_run_failed:
        print("SUMMARY: per-run invariant FAILURES present -- see above.")
    if not det_ok:
        print("SUMMARY: determinism FAILURES present -- see above.")
    if not pressure_ok:
        print("SUMMARY: asymmetric-visibility pressure test FAILED -- see above.")
    if not any_run_failed and det_ok and pressure_ok:
        print("SUMMARY: all runs passed, all determinism checks passed, pressure test passed.")
        sys.exit(0)
    sys.exit(1)


if __name__ == "__main__":
    main()
