"""
Spatial-aliasing diagnostic: the classic failure mode this is checking for is a cell-based merge
scheme that snaps/gravitates computational particles toward cell centers, which stamps a
Dirac-comb-like structure onto the particle distribution -- visible as spurious excess power at
the cell spatial frequency (and harmonics) in the Fourier spectrum of the particle density, even
though the scheme conserves mass exactly. Our scheme merges to the actual weighted centroid of a
specific nearest-neighbor PAIR (never snapped to any grid point), so there's no obvious mechanism
for this -- but the crowding TRIGGER is still cell-based (a pair is only eligible once its cell
exceeds a threshold), so it's worth checking empirically rather than assuming continuity of
merge position rules out any imprint.

Protocol (as specified): repeatedly (1) inject a fresh batch of uniformly-sampled particles into
whatever currently exists, (2) merge (old and new together), (3) repeat for many generations --
mimicking a real simulation where new particles are continuously born (e.g. by ionization) and
periodically thinned by merging. After many generations, check the SURVIVING particle
distribution's spatial power spectrum for excess power at the cell frequency (or any other
periodicity) relative to a null reference (a uniformly-random point set of the same count, never
merged at all).

Deliberately single-patch, single-rank: this isolates the pure "cell-based crowding trigger +
nearest-neighbor merge" effect from the cross-rank propose/judge protocol (already characterized
separately for boundary bias) -- with only one patch, nothing is ever boundary-exposed, so
everything resolves via the trivial tier alone. If aliasing exists at all, it must come from the
crowding-trigger/NN-merge mechanism itself, not from anything rank- or patch-boundary-specific.
"""

import argparse
import json
import random
import sys

import numpy as np
from mpi4py import MPI

import merge_lib as ml
from run_scenario import run_one_round


def radial_power_spectrum(hist, bin_width, n_bins_radial=40):
    """2D FFT power spectrum of a density histogram, radially averaged. Returns (k_values, power)
    where k is a spatial frequency in cycles/unit-length (not angular)."""
    n = hist.shape[0]
    f = np.fft.fft2(hist - hist.mean())
    power2d = np.abs(f) ** 2
    freqs = np.fft.fftfreq(n, d=bin_width)  # cycles per unit length
    kx, ky = np.meshgrid(freqs, freqs, indexing="ij")
    kmag = np.sqrt(kx ** 2 + ky ** 2)

    k_max = freqs.max()
    k_edges = np.linspace(0, k_max, n_bins_radial + 1)
    k_centers = 0.5 * (k_edges[:-1] + k_edges[1:])
    power_radial = np.zeros(n_bins_radial)
    counts = np.zeros(n_bins_radial)
    flat_k = kmag.ravel()
    flat_p = power2d.ravel()
    idx = np.digitize(flat_k, k_edges) - 1
    for b in range(n_bins_radial):
        mask = idx == b
        if mask.any():
            power_radial[b] = flat_p[mask].mean()
            counts[b] = mask.sum()
    return k_centers, power_radial


def make_density_hist(positions, domain_size, n_hist_bins):
    hist, _, _ = np.histogram2d(
        [p[0] for p in positions], [p[1] for p in positions],
        bins=n_hist_bins, range=[[0, domain_size], [0, domain_size]]
    )
    return hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain-size", type=float, default=4.0)
    ap.add_argument("--cell-dx", type=float, default=0.2)
    ap.add_argument("--thresh", type=int, default=2)
    ap.add_argument("--n-inject", type=int, default=250)
    ap.add_argument("--n-generations", type=int, default=60)
    ap.add_argument("--hist-bins-per-cell", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    assert size == 1, "single-patch aliasing test is deliberately serial"

    patches = [((0.0, 0.0), (args.domain_size, args.domain_size), 0)]
    my_patches = [0]
    local = {0: []}

    rng = random.Random(args.seed)
    shuffle_rng = random.Random(args.seed + 999)
    id_counter = [1]

    n_over_trace = []
    n_total_trace = []

    for gen in range(args.n_generations):
        box = patches[0]
        for _ in range(args.n_inject):
            pos = tuple(rng.uniform(0.0, args.domain_size) for _ in range(2))
            local[0].append(ml.make_particle(pos, 1.0, 0, 0, rng, id_hint=id_counter[0]))
            id_counter[0] += 1

        local, n_trivial, n_judge = run_one_round(
            comm, rank, size, patches, my_patches, local, args.thresh, 0.0, args.cell_dx,
            shuffle_rng, id_counter)

        n_total_trace.append(len(local[0]))
        lc = {}
        for p in local[0]:
            k = ml.cell_key(p["pos"], args.cell_dx)
            lc[k] = lc.get(k, 0) + 1
        n_over_trace.append(sum(1 for c in lc.values() if c > args.thresh))

    final = local[0]
    n_final = len(final)
    total_mass = sum(p["weight"] for p in final)

    # ---- Null reference: same COUNT of purely uniform-random points, never merged ----
    null_rng = random.Random(args.seed + 424242)
    null_positions = [(null_rng.uniform(0, args.domain_size), null_rng.uniform(0, args.domain_size))
                       for _ in range(n_final)]

    n_hist_bins = int(round(args.domain_size / args.cell_dx * args.hist_bins_per_cell))
    bin_width = args.domain_size / n_hist_bins

    final_positions = [p["pos"] for p in final]
    final_weights = [p["weight"] for p in final]

    hist_count_final = make_density_hist(final_positions, args.domain_size, n_hist_bins)
    hist_count_null = make_density_hist(null_positions, args.domain_size, n_hist_bins)

    # Weighted (charge) density histogram for the final (merged) set.
    hist_mass_final, _, _ = np.histogram2d(
        [p[0] for p in final_positions], [p[1] for p in final_positions],
        bins=n_hist_bins, range=[[0, args.domain_size], [0, args.domain_size]],
        weights=final_weights,
    )

    k_final_count, p_final_count = radial_power_spectrum(hist_count_final, bin_width)
    k_null, p_null = radial_power_spectrum(hist_count_null, bin_width)
    k_final_mass, p_final_mass = radial_power_spectrum(hist_mass_final, bin_width)

    k_cell = 1.0 / args.cell_dx  # cycles per unit length corresponding to the cell spacing

    report = {
        "n_generations": args.n_generations, "n_inject": args.n_inject,
        "n_final": n_final, "total_mass": total_mass,
        "n_total_trace": n_total_trace, "n_over_trace": n_over_trace,
        "cell_dx": args.cell_dx, "k_cell": k_cell,
        "k_values": k_final_count.tolist(),
        "power_final_count": p_final_count.tolist(),
        "power_null_count": p_null.tolist(),
        "power_final_mass": p_final_mass.tolist(),
    }
    with open(args.out, "w") as f:
        json.dump(report, f, indent=2)

    print(f"n_final={n_final} total_mass={total_mass:.1f} "
          f"n_total_trace(last 10)={n_total_trace[-10:]} n_over_trace(last 10)={n_over_trace[-10:]}")
    print(f"cell frequency k_cell = {k_cell:.3f} cycles/unit-length (cell_dx={args.cell_dx})")
    print()
    print(f"{'k':>8s} {'power_final_count':>18s} {'power_null_count':>17s} {'ratio':>8s} {'power_final_mass':>17s}")
    for i in range(len(k_final_count)):
        k = k_final_count[i]
        pf, pn, pm = p_final_count[i], p_null[i], p_final_mass[i]
        ratio = pf / pn if pn > 0 else float("nan")
        marker = " <== cell freq" if abs(k - k_cell) < (k_final_count[1] - k_final_count[0]) else ""
        print(f"{k:8.3f} {pf:18.2f} {pn:17.2f} {ratio:8.3f} {pm:17.2f}{marker}")


if __name__ == "__main__":
    main()
