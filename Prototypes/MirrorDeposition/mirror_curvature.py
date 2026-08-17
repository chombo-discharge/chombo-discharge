#!/usr/bin/env python3
"""
The curvature Jacobian, and how to get it without a Hessian.

The mirror needs J = prod_i (1 - k_i d)/(1 + k_i d) on a curved boundary, because
reflection is measure-preserving only about a plane. The design question is not
whether J is needed -- mirror_sphere_ext.py measures that -- but what has to be
stored to evaluate it.

Three routes are compared here, and they must agree:

  closed form    prod_i (1 - k_i d)/(1 + k_i d), from principal curvatures known
                 analytically for each test surface. The reference.
  determinant    |det grad R| of the reflection map built from the implicit
                 function. No curvature anywhere -- just a 3x3 finite difference
                 of the map itself.
  invariants     (1 - 2H d + K d^2)/(1 + 2H d + K d^2), where M = grad(nhat) at
                 the FOOT POINT, 2H = tr M, and K = 0.5[(tr M)^2 - tr(M^2)].

The third is what the plan stores: two Reals per band cell, computed once at
regrid. M has nhat in its kernel, so its two nonzero eigenvalues are exactly
k_1, k_2 and the second invariant is exactly the Gaussian curvature -- which is
why no eigen-decomposition and no principal directions are needed.

A torus is the load-bearing test. A sphere has k_1 = k_2 and a cylinder has
k_2 = 0, so both are satisfied by wrong implementations; a torus has two
unequal curvatures that vary over the surface and change relative sign between
its outer and inner regions.

SCOPE, and it is the important one: everything here is CONTINUUM. The implicit
functions are analytic and the finite differences use steps of 1e-4 to 1e-6.
K is a second derivative of the geometry, and nothing below measures what it
does when grad(nhat) is differenced on the actual mesh, at the actual dx, from
a discretely sampled IF whose own error is O(dx^2). That can only be measured
in the code, and it is the reason the exact form could still disappoint.

Usage:
    python3 mirror_curvature.py            # all checks
    python3 mirror_curvature.py exponent   # only the D-1 exponent trap
    python3 mirror_curvature.py singular   # only the concave denominator
"""

import sys

import numpy as np

FD_MAP = 1e-4  # step for differencing the reflection map
FD_GRAD = 1e-6  # step for differencing the implicit function


# --------------------------------------------------------------------------- #
# The three routes
# --------------------------------------------------------------------------- #
def reflect(implicit_fn, grad_fn, x):
    """Reflection of x through the level set implicit_fn == 0."""
    g = grad_fn(x)
    gn = np.linalg.norm(g)
    return x - 2.0 * (implicit_fn(x) / gn) * (g / gn)


def foot_point(implicit_fn, grad_fn, x):
    """Nearest point on the surface, to first order."""
    g = grad_fn(x)
    gn = np.linalg.norm(g)
    return x - (implicit_fn(x) / gn) * (g / gn)


def signed_distance(implicit_fn, grad_fn, x):
    return implicit_fn(x) / np.linalg.norm(grad_fn(x))


def jacobian_by_determinant(implicit_fn, grad_fn, x, h=FD_MAP):
    """|det grad R| by differencing the reflection map. No curvature needed."""
    cols = [
        (reflect(implicit_fn, grad_fn, x + h * e) - reflect(implicit_fn, grad_fn, x - h * e))
        / (2.0 * h)
        for e in np.eye(3)
    ]
    return abs(np.linalg.det(np.stack(cols, axis=1)))


def curvature_invariants(implicit_fn, grad_fn, x, h=FD_MAP):
    """
    (2H, K) at the foot point of x.

    Evaluating at the foot point rather than at x is not cosmetic. For a signed
    distance function div(nhat)(y) = sum_i k_i/(1 + k_i d(y)), so differencing at
    a band cell a distance d out returns the curvature at that offset, not at the
    surface -- roughly a 17% error at d ~ dx, R ~ 6 dx, in the very quantity the
    correction is built from.
    """
    p = foot_point(implicit_fn, grad_fn, x)

    def nhat(y):
        g = grad_fn(y)
        return g / np.linalg.norm(g)

    cols = [(nhat(p + h * e) - nhat(p - h * e)) / (2.0 * h) for e in np.eye(3)]
    m = np.stack(cols, axis=1)
    two_h = np.trace(m)
    k = 0.5 * (two_h**2 - np.trace(m @ m))
    return float(two_h), float(k)


def jacobian_from_invariants(two_h, k, d):
    """J = (1 - 2H d + K d^2)/(1 + 2H d + K d^2), the product form expanded."""
    num = 1.0 - two_h * d + k * d * d
    den = 1.0 + two_h * d + k * d * d
    return abs(num / den), den


def principal_curvatures(two_h, k):
    """
    Roots of t^2 - 2H t + K = 0. The discriminant is (k_1 - k_2)^2, so it is
    non-negative analytically and only goes slightly negative through finite
    differences -- clamp rather than admitting a complex pair.
    """
    disc = max(two_h * two_h - 4.0 * k, 0.0)
    root = np.sqrt(disc)
    return 0.5 * (two_h + root), 0.5 * (two_h - root)


# --------------------------------------------------------------------------- #
# Test surfaces. Each returns (name, implicit_fn, grad_fn, sampler, exact_J)
# --------------------------------------------------------------------------- #
def make_sphere(radius, fluid_outside):
    centre = np.zeros(3)
    sign = 1.0 if fluid_outside else -1.0

    def implicit_fn(y):
        return sign * (np.linalg.norm(y - centre) - radius)

    def grad_fn(y):
        v = y - centre
        return sign * v / np.linalg.norm(v)

    def sampler(rng, d):
        u = rng.normal(size=3)
        u /= np.linalg.norm(u)
        r = radius + sign * d
        return centre + r * u

    def exact_j(x, d):
        # both principal curvatures are sign/radius
        kap = sign / radius
        return abs(((1 - kap * d) / (1 + kap * d)) ** 2)

    side = "convex" if fluid_outside else "concave"
    return (f"sphere R={radius} ({side})", implicit_fn, grad_fn, sampler, exact_j)


def make_cylinder(radius, fluid_outside):
    axis = np.zeros(2)
    sign = 1.0 if fluid_outside else -1.0

    def implicit_fn(y):
        return sign * (np.linalg.norm(y[:2] - axis) - radius)

    def grad_fn(y):
        v = y[:2] - axis
        r = np.linalg.norm(v)
        return sign * np.array([v[0] / r, v[1] / r, 0.0])

    def sampler(rng, d):
        th = rng.uniform(0.0, 2.0 * np.pi)
        r = radius + sign * d
        return np.array([r * np.cos(th), r * np.sin(th), rng.uniform(-5.0, 5.0)])

    def exact_j(x, d):
        # principal curvatures are (sign/radius, 0): the product has ONE factor
        kap = sign / radius
        return abs((1 - kap * d) / (1 + kap * d))

    side = "convex" if fluid_outside else "concave"
    return (f"cylinder R={radius} ({side})", implicit_fn, grad_fn, sampler, exact_j)


def make_torus(ring_radius, tube_radius):
    def implicit_fn(y):
        q = np.hypot(y[0], y[1]) - ring_radius
        return np.hypot(q, y[2]) - tube_radius

    def grad_fn(y):
        return np.array(
            [
                (implicit_fn(y + FD_GRAD * e) - implicit_fn(y - FD_GRAD * e)) / (2.0 * FD_GRAD)
                for e in np.eye(3)
            ]
        )

    def sampler(rng, d):
        th = rng.uniform(0.0, 2.0 * np.pi)
        ph = rng.uniform(0.0, 2.0 * np.pi)
        rr = tube_radius + d
        return np.array(
            [
                (ring_radius + rr * np.cos(ph)) * np.cos(th),
                (ring_radius + rr * np.cos(ph)) * np.sin(th),
                rr * np.sin(ph),
            ]
        )

    def exact_j(x, d):
        # k_tube = 1/r everywhere; k_ring = cos(phi)/(R + r cos phi), which is
        # NEGATIVE on the inner half, so this surface is a saddle there.
        ph = np.arctan2(x[2], np.hypot(x[0], x[1]) - ring_radius)
        k1 = 1.0 / tube_radius
        k2 = np.cos(ph) / (ring_radius + tube_radius * np.cos(ph))
        return abs(((1 - k1 * d) / (1 + k1 * d)) * ((1 - k2 * d) / (1 + k2 * d)))

    return (f"torus R={ring_radius}, r={tube_radius}", implicit_fn, grad_fn, sampler, exact_j)


SURFACES = [
    make_sphere(8.0, True),
    make_sphere(6.0, False),
    make_cylinder(5.0, True),
    make_cylinder(6.0, False),
    make_torus(6.0, 2.5),
]


# --------------------------------------------------------------------------- #
# Checks
# --------------------------------------------------------------------------- #
def check_agreement(n_per_surface=5, seed=11):
    print("### Do the three routes agree?  (continuum, analytic implicit functions)\n")
    worst = 0.0
    for name, implicit_fn, grad_fn, sampler, exact_j in SURFACES:
        print(f"--- {name}")
        print(
            f"    {'d':>6} {'closed form':>13} {'|det grad R|':>14} "
            f"{'from 2H,K':>12} {'2H':>10} {'K':>11} {'worst rel':>11}"
        )
        rng = np.random.default_rng(seed)
        for _ in range(n_per_surface):
            d = rng.uniform(0.1, 1.8)
            x = sampler(rng, d)
            d_meas = abs(signed_distance(implicit_fn, grad_fn, x))
            ref = exact_j(x, d_meas)
            det = jacobian_by_determinant(implicit_fn, grad_fn, x)
            two_h, k = curvature_invariants(implicit_fn, grad_fn, x)
            inv, _ = jacobian_from_invariants(two_h, k, d_meas)
            rel = max(abs(det - ref), abs(inv - ref)) / ref
            worst = max(worst, rel)
            print(
                f"    {d_meas:6.3f} {ref:13.6f} {det:14.6f} {inv:12.6f} "
                f"{two_h:10.5f} {k:11.5f} {rel:11.2e}"
            )
        print()
    print(f"worst relative disagreement across all surfaces: {worst:.2e}")
    print(
        "PASS -- determinant and invariant routes both reproduce the closed form"
        if worst < 1e-3
        else "FAIL"
    )
    print(
        "\nNote what this does NOT show: every number above is continuum. On a mesh,\n"
        "grad(nhat) is a second difference of discretely sampled geometry. Whether K\n"
        "survives that is unmeasured and is the main risk in shipping the exact form."
    )


def check_exponent_trap(seed=5):
    """
    Reading (r'/r)^2 off a sphere and hardcoding SpaceDim-1 is the natural
    mistake. It is right only for a sphere.
    """
    print("\n\n### The D-1 exponent trap: ratio^(D-1) against the product form\n")
    print(
        "    The mistake is reading J = (r'/r)^2 off a sphere and implementing\n"
        "    ((1 - k d)/(1 + k d))^(SpaceDim-1) with a SINGLE curvature k. That is\n"
        "    right whenever k_1 == k_2 and wrong otherwise.\n"
    )
    print(
        f"    {'surface':>26} {'k_1':>8} {'k_2':>8} {'product':>10} "
        f"{'ratio^(D-1)':>12} {'error':>9}"
    )
    for name, implicit_fn, grad_fn, sampler, exact_j in SURFACES:
        rng = np.random.default_rng(seed)
        x = sampler(rng, 0.8)
        d = abs(signed_distance(implicit_fn, grad_fn, x))
        two_h, k = curvature_invariants(implicit_fn, grad_fn, x)
        k1, k2 = principal_curvatures(two_h, k)
        correct = jacobian_by_determinant(implicit_fn, grad_fn, x)
        dominant = k1 if abs(k1) >= abs(k2) else k2
        naive = abs(((1 - dominant * d) / (1 + dominant * d)) ** (3 - 1))
        print(
            f"    {name:>26} {k1:8.4f} {k2:8.4f} {correct:10.6f} "
            f"{naive:12.6f} {(naive - correct) / correct:+9.1%}"
        )
    print(
        "\n    The exponent is correct for a sphere and wrong for everything else.\n"
        "    In 2-D a circle has ONE principal curvature and SpaceDim-1 = 1, so the\n"
        "    exponent coincides with the product form exactly -- a 2-D-only test\n"
        "    cannot catch this. CI builds 2-D."
    )


def check_concave_singularity():
    """
    For concave geometry 2H < 0, so the denominator 1 + 2H d + K d^2 can vanish.
    The reflect band allows d up to 3*s_max ~ 2.6 dx, so this is reachable on a
    tight cavity and needs a guard rather than a division.
    """
    print("\n\n### Concave denominator: where does 1 + 2H d + K d^2 vanish?\n")
    print(
        f"    {'cavity radius':>15} {'k_1':>9} {'k_2':>9} "
        f"{'d at which den = 0':>20} {'reachable in band?':>21}"
    )
    band_max = 1.5 * np.sqrt(3.0)  # (3/2)*dx*sum|n_i|, worst case, dx = 1
    for radius in (10.0, 6.0, 4.0, 3.0, 2.5, 2.0):
        _, implicit_fn, grad_fn, sampler, _ = make_sphere(radius, False)
        rng = np.random.default_rng(1)
        x = sampler(rng, 0.5)
        two_h, k = curvature_invariants(implicit_fn, grad_fn, x)
        k1, k2 = principal_curvatures(two_h, k)
        # den = (1 + k_1 d)(1 + k_2 d), so it vanishes at d = -1/k_i for k_i < 0
        zeros = [-1.0 / kk for kk in (k1, k2) if kk < -1e-12]
        d0 = min(zeros) if zeros else float("inf")
        flag = "YES -- needs a guard" if d0 <= band_max else "no"
        print(f"    {radius:15.1f} {k1:9.4f} {k2:9.4f} {d0:20.3f} {flag:>21}")
    print(
        f"\n    Band reaches d = {band_max:.3f} dx on a body-diagonal normal. A cavity of\n"
        "    radius <~ 3 dx therefore puts the singularity inside the reflect band.\n"
        "    The image position is degenerate there too (it collapses to the centre),\n"
        "    so this is a genuine breakdown of the reflection, not just of J."
    )


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "agree"):
        check_agreement()
    if what in ("all", "exponent"):
        check_exponent_trap()
    if what in ("all", "singular"):
        check_concave_singularity()
