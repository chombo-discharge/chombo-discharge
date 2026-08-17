# DELETE BEFORE MERGING TO MAIN

These are throwaway Python prototypes used to design mirrored cut-cell deposition
for issue #29. They are not part of the build, are not run by CI, and must be
removed before this branch merges. They are committed only so the numbers behind
the design decisions can be reproduced and challenged while the work is in review.

## What is here

The first three were written for revision 1 of the plan. The `_ext` / `curvature`
files were written for revision 2, after an adversarial review refuted three of
revision 1's quantitative conclusions; each exists because the original harness
could not answer the question that decided a design choice.

| Script | Geometry | What it establishes |
|---|---|---|
| `mirror_test.py` | planar EB, arbitrary normal | the mirror is unbiased at every kappa and every normal; NGP and plain CIC are not; containment and reach |
| `mirror_sphere.py` | sphere, convex and concave | the reflection error on a curved surface, for three reflection sources, with no Jacobian |
| `mirror_cylinder.py` | cylinder, convex and concave | the Jacobian is a product over principal curvatures, **not** an exponent `D-1` |
| `mirror_curvature.py` | sphere, cylinder, **torus** | the exact Jacobian needs no Hessian: two curvature invariants suffice, and they agree with `\|det grad R\|` and with the closed form |
| `mirror_sphere_ext.py` | sphere, both signs, R = 4..40 | reflection source **crossed with** Jacobian; the mass check is an identity; kappa sub-sampling is not the explanation |
| `mirror_planar_ext.py` | planar EB | the reflect-criterion table; the band width, measured rather than asserted; margin sensitivity; the broken `lattice` mode |

All are Cartesian: uniform cubic cells, tensor-product CIC, 16^3 valid cells with
2 ghost cells. Only the shape of the embedded boundary varies. Volume fractions
are analytic where a closed form exists (plane) and sub-sampled otherwise.

## Running

```
python3 mirror_test.py random 1600000
python3 mirror_sphere.py
python3 mirror_cylinder.py

python3 mirror_curvature.py              # agreement, exponent trap, concave singularity
python3 mirror_sphere_ext.py radii       # the main cross
python3 mirror_sphere_ext.py nsub        # is the mass ratio a kappa artifact? (no)
python3 mirror_sphere_ext.py identity    # the mass check is an identity
python3 mirror_planar_ext.py criteria    # reflect-criterion table
python3 mirror_planar_ext.py band        # measured band width
python3 mirror_planar_ext.py margin      # scoring-margin sensitivity
python3 mirror_planar_ext.py lattice     # why `lattice` mode must not be used
```

Do **not** use `python3 mirror_test.py lattice`. See harness bug 3.

## Headline numbers

* Planar EB, all normals: the mirror is unbiased where CIC and NGP are not.
  Aggregated over cases the mean is within 0.2%; **per case and per bin the
  spread is -1.64% to +1.32%** at 1.6M particles. The aggregate does not bound
  the individual cases, and revision 1 of the plan quoted it as though it did.
* The reflect band is `(3/2)*dx*sum|n_i| = 3*s_max`, and that bound is tight to
  within 1%. The narrower `dx*(1 + sum|n_i|/2)` is short by up to 0.71 dx on
  tilted normals and agrees only for axis-aligned ones.
* Sphere, **exact** reflection, no Jacobian: +4% at R = 40 dx, +16% at 16,
  +40% at 8, +117% at 4. With the exact Jacobian every radius and both curvature
  signs land inside ~1.5%. No improvement to the reflection source removes this,
  and no Jacobian removes the reflection-source error.
* The cell-constant plane is **not** a minor term. With the Jacobian applied it
  costs 7.3% at R = 8 dx, 9.0% at 6 and 17.7% at 4, against ~1% for exact
  reflection. Revision 1 concluded the opposite from the `kappa <= 0.05` column
  alone, measured while the uncorrected curvature error ran to +1.17.
* `sum(kappa phi dV)/sum(m)` is an **identity**, equal to
  `mean_p[kappatilde(x_p) + J_p*kappatilde(R(x_p))]` to machine precision. On a
  plane it is 1 by construction; on a curved surface it is the kappa-weighted
  volume average of the density error. It is not independent evidence, and its
  curved-geometry values are not a kappa quadrature artifact.
* Cylinder: the product form lands inside 0.003. `(r'/r)^(D-1)` over-corrects a
  cylinder by -0.24 to +0.21 and a torus by -44%.

## Four harness bugs found while doing this

All produced confident, plausible, wrong numbers, and none was caught by reading
the code -- each was caught by a control or by the output contradicting its own
prose. Anyone re-running or extending these should assume more of the same.

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

3. **`mirror_test.py`'s `lattice` mode is unsound**, and its docstring recommends
   it as the shot-noise-free way to measure bias. At the default `sub=4` it
   reports the mirror as **+20.9% biased** in the kappa = 0.63 cut cells of the
   axis-aligned case, with a standard error of exactly `0.00000`. It is
   quantization of the fluid boundary by the sampling lattice: the deviation runs
   -25.5% at sub=2, +20.9% at 4, -0.9% at 8 and 16. Run
   `mirror_planar_ext.py lattice` to see it. The zero standard error makes this
   the most convincing of the four and the least true. Not fixed -- use `random`.

4. **Two bugs in `mirror_curvature.py` while writing it.** The exponent-trap
   check compared the image and the particle against their shared foot point, so
   its "ratio" was identically 1 and every surface scored a 0% error; and the
   concave-singularity check rejected the double root of `1 + 2H d + K d^2` as
   complex, because a `1e-9` imaginary tolerance is too tight for a discriminant
   that is analytically zero. Both reported clean output that contradicted the
   prose next to it. Both fixed by going through the principal curvatures
   extracted from the invariants rather than through geometry or `np.roots`.

## Scoring caveat that is easy to get wrong

Score only cells at least **3** cells inside the sampled region. A cell centre on
the solid side is reconstructed almost entirely from images, i.e. from particles
near `R(c)`, which sits up to `sqrt(3)*dx` away. A 1-cell margin leaves those
cells under-supplied and fabricates a 1.7% deficit at small kappa that survives a
16x increase in particle count and looks exactly like a real bias. The
requirement is `1 + sqrt(3) = 2.73 -> 3`.

Revision 1 of this file added "and the measured residual vanishes precisely
there". It does not. `mirror_planar_ext.py margin` shows the margin-1 artifact is
real for the two `s max` cases (-0.0177 -> -0.0009 and -0.0125 -> +0.0020), which
justifies margin >= 2, but there is no plateau at 3: on `near lo, tilted` the
kappa <= 0.05 deviation runs -0.007, -0.011, -0.019, -0.025 for margins 1 to 4
while the sample shrinks from 235 to 24 cut cells. Margin 3 is a defensible
conservative choice, not a measured floor.

## Known defects in these harnesses, not yet fixed

* `mirror_sphere.py` has **no Jacobian** and runs R = 4, 2.5, 6 only. Use
  `mirror_sphere_ext.py` for anything involving J or other radii.
* `mirror_cylinder.py` runs four cases, of which **R=12 produces zero particles**
  (the radius exceeds the box half-diagonal) and **R=8 produces zero scored
  cells** (its cut cells fall outside the margin-3 mask). Both fail silently and
  print an empty table.
* `mirror_test.py` prints `FAIL` at 1.6M particles because the mirror reaches 3
  cells against 2 nominal ghosts. That verdict is about source-patch deposition,
  which the plan superseded with remap; it is not a failure of the scheme.
* `mirror_test.py`'s `near lo edge` case scores **zero cut cells** at margin >= 3,
  because its only cut cells sit at index 2. One of the two cases added to probe
  near-domain-edge behaviour therefore tests nothing.

## Scope these do not cover

The normals here are analytic. A real `ebisbox.normal(vof)` is reconstructed from
discrete area fractions and degrades as kappa -> 0, which is exactly where the
mirror matters. These harnesses isolate the geometric approximation; the
discrete-normal reconstruction error can only be measured in the code.

**`mirror_curvature.py` is entirely continuum.** Its implicit functions are
analytic and its finite differences use steps of 1e-4 to 1e-6. `K` is a second
derivative of the geometry, and nothing here measures what `tr M` and
`0.5[(tr M)^2 - tr(M^2)]` do when `grad(nhat)` is differenced on the actual mesh,
at the actual dx, from a discretely sampled implicit function whose own error is
`O(dx^2)` -- at R = 6 dx, where the correction is worth 56%. If `K` is noise at
mesh scale, the exact Jacobian is no better than the linearized one. That can
only be measured in the code and it is the largest remaining risk in the scheme.

Everything is single-patch and single-level, so nothing here can see the
level-routing behaviour of `ParticleContainer::remap`, which promotes a particle
to the finest level containing its position.
