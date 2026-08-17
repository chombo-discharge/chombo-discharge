# DELETE BEFORE MERGING TO MAIN

These are throwaway Python prototypes used to design mirrored cut-cell deposition
for issue #29. They are not part of the build, are not run by CI, and must be
removed before this branch merges. They are committed only so the numbers behind
the design decisions can be reproduced and challenged while the work is in review.

## What is here

| Script | Geometry | What it establishes |
|---|---|---|
| `mirror_test.py` | planar EB, arbitrary normal | the mirror is unbiased at every kappa and every normal; NGP and plain CIC are not; containment and reach |
| `mirror_sphere.py` | sphere, convex and concave | the cell-constant plane is *not* the dominant error; the method needs a curvature Jacobian |
| `mirror_cylinder.py` | cylinder, convex and concave | the Jacobian is a product over principal curvatures, **not** an exponent `D-1` |

All three are Cartesian: uniform cubic cells, tensor-product CIC, 16^3 valid
cells with 2 ghost cells. Only the shape of the embedded boundary varies.
Volume fractions are analytic where a closed form exists (plane) and
sub-sampled otherwise (sphere, cylinder).

## Running

```
python3 mirror_test.py random 1600000
python3 mirror_sphere.py
python3 mirror_cylinder.py
```

Each prints per-kappa-bin mean deviation of the deposited density from the exact
answer of 1, with a standard error, plus a mass check
`sum(kappa*phi*dV) / sum(m)`.

## Headline numbers

* Planar EB, all normals: mirror is unbiased to 0.2%, mass ratio 1.00001.
  Plain CIC runs -0.98 to -0.11; NGP runs -0.998 to +0.09.
* Sphere: without a curvature Jacobian the mirror overshoots by +0.04 at
  R = 40 dx and +0.82 at R = 4 dx. With it, every radius lands inside 0.15%.
* Cylinder: the product form `prod_i (1 - k_i d)/(1 + k_i d)` lands inside
  0.003. Using `(r'/r)^(D-1)` instead over-corrects by -0.24 to +0.21.

## Two harness bugs found while doing this

Both produced confident, plausible, wrong numbers, and both were caught only by
a control case rather than by reading the code. Anyone re-running or extending
these should assume more of the same.

1. **Catastrophic cancellation in the plane-cube volume fraction.** The
   inclusion-exclusion sum adds terms of order `d^3` to produce an answer of
   order 1, so cells far from the plane returned `kappa ~ 1e-13` instead of
   exactly 0. Deep-solid cells then counted as fluid and the reported reach blew
   up to 21 cells. Fixed by deciding all-solid / all-fluid cells from the corner
   extrema first. A real `EBISBox` returns exact covered flags, so this is a
   prototype-only failure mode.

2. **Wrong fluid volume in the sphere harness.** `Vbox - Vsphere` is correct only
   when the sphere lies wholly inside the box; for a large radius it is hugely
   negative, so the particle weight was garbage. Caught by a control at
   `R = 200 dx`, where a sphere is effectively a plane and the answer therefore
   had to reproduce the planar result. Fixed by deriving the fluid volume from
   the sampling acceptance rate instead.

## Scoring caveat that is easy to get wrong

Score only cells at least **3** cells inside the sampled region. A cell centre on
the solid side is reconstructed almost entirely from images, i.e. from particles
near `R(c)`, which sits up to `sqrt(3)*dx` away. A 1-cell margin leaves those
cells under-supplied and fabricates a 1.7% deficit at small kappa that survives a
16x increase in particle count and looks exactly like a real bias. The
requirement is `1 + sqrt(3) = 2.73 -> 3`, and the measured residual vanishes
precisely there.

## Scope these do not cover

The normals here are analytic. A real `ebisbox.normal(vof)` is reconstructed from
discrete area fractions and degrades as kappa -> 0, which is exactly where the
mirror matters. These harnesses isolate the geometric approximation; the
discrete-normal reconstruction error can only be measured in the code.
