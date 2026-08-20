"""The tangent frame is part of the stored format -- and a frame-free form removes it.

Revision 6, written for review 6. PLAN.md 2.1 stores the shape operator as 3 Reals, "2x2
symmetric, tangent frame"; 3.3 fits it in "any orthonormal basis perpendicular to n_i" and
2.2 reflects in "t1, t2 spanning the plane perpendicular to n_c". Those must be the SAME
basis. No other harness in this directory can see a mismatch, because every one of them
uses an ISOTROPIC shape operator: mirror_source.quad_reflect takes a scalar kap_princ, and
S = c*I is frame-invariant. The fitted anisotropic 2x2 S from mirror_discrete_curvature is
never fed into a reflection anywhere.

This file closes that gap on the cheapest anisotropic surface -- a cylinder, principal
curvatures (1/R, 0) -- and measures two things:

  1. what a frame mismatch costs, and that tr S and det S (hence J, hence the 6.1 mass
     check) do not move when it happens;
  2. that a frame-FREE formulation, storing the shape operator as the 3x3 symmetric
     world-frame tensor Sw = P S P^T with P = [t1 t2], reproduces the frame formulation
     to machine precision:

         d    = eta + 0.5 * w . (Sw w)              w = p - x_c,  eta = n_c . w
         nhat = normalize(n_c + Sw w)
         2H   = tr Sw          K = 0.5[(tr Sw)^2 - tr(Sw^2)]

     Sw n_c = 0 by construction, so the eta direction contributes nothing and no tangent
     projection is needed at any point.

Run: python3 mirror_frame.py
"""
import numpy as np

RNG = np.random.default_rng(3)


def tangent_frame(n):
    """The rule both existing harnesses use (mirror_discrete_curvature.py:125-130,
    mirror_source.py:61-65): a = e_x unless |n_x| > 0.9, then e_y."""
    a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) <= 0.9 else np.array([0.0, 1.0, 0.0])
    t1 = np.cross(n, a)
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(n, t1)

    return t1, t2


def cylinder_shape_operator(n, t1, t2, R):
    """S = +dn_hat of a cylinder of radius R about z, in the (t1,t2) frame. The normal is
    invariant along z, so dn projects onto the circumferential direction only: c = 1/R
    around, 0 along the axis."""
    th = np.cross(np.array([0.0, 0.0, 1.0]), n)
    M = np.outer(th, th) / R

    S = np.zeros((2, 2))
    for i, ti in enumerate((t1, t2)):
        for j, tj in enumerate((t1, t2)):
            S[i, j] = ti @ (M @ tj)

    return S


def reflect_framed(p, xc, nc, S, t1, t2):
    """PLAN.md 2.2 as written: a 2x2 S read in an explicit tangent frame."""
    w = p - xc
    xi = np.array([t1 @ w, t2 @ w])
    d = nc @ w + 0.5 * xi @ S @ xi
    Sxi = S @ xi
    nh = nc + Sxi[0] * t1 + Sxi[1] * t2
    nh /= np.linalg.norm(nh)

    return d, p - 2.0 * d * nh


def reflect_frameless(p, xc, nc, Sw):
    """The frame-free form: Sw is the 3x3 symmetric world-frame shape operator."""
    w = p - xc
    d = nc @ w + 0.5 * w @ (Sw @ w)
    nh = nc + Sw @ w
    nh /= np.linalg.norm(nh)

    return d, p - 2.0 * d * nh


def run(R=6.0, ntrial=4000):
    print(f"### Cylinder R = {R:.1f} dx about z, fluid outside. c1 = 1/R, c2 = 0.")
    print("### Particles up to 2.6 dx out (PLAN.md 1.2's CIC reach) and +-2 dx tangentially.\n")

    rows = []
    for _ in range(ntrial):
        th = RNG.uniform(0.0, 2.0 * np.pi)
        xc = np.array([R * np.cos(th), R * np.sin(th), RNG.uniform(-3.0, 3.0)])
        nc = np.array([np.cos(th), np.sin(th), 0.0])

        t1, t2 = tangent_frame(nc)
        p = xc + RNG.uniform(0.05, 2.6) * nc + RNG.uniform(-2.0, 2.0) * t1 + RNG.uniform(-2.0, 2.0) * t2

        S = cylinder_shape_operator(nc, t1, t2, R)
        P = np.column_stack([t1, t2])
        Sw = P @ S @ P.T

        d_ref, img_ref = reflect_framed(p, xc, nc, S, t1, t2)
        d_free, img_free = reflect_frameless(p, xc, nc, Sw)

        # A DIFFERENT but equally legal orthonormal tangent basis: rotate within the plane.
        phi = RNG.uniform(0.0, 2.0 * np.pi)
        u1 = np.cos(phi) * t1 + np.sin(phi) * t2
        u2 = -np.sin(phi) * t1 + np.cos(phi) * t2
        d_mis, img_mis = reflect_framed(p, xc, nc, S, u1, u2)
        S_u = cylinder_shape_operator(nc, u1, u2, R)

        rows.append((abs(d_free - d_ref),
                     np.linalg.norm(img_free - img_ref),
                     abs(d_mis - d_ref),
                     np.linalg.norm(img_mis - img_ref),
                     abs(np.trace(S_u) - np.trace(S)),
                     abs(np.linalg.det(S_u) - np.linalg.det(S)),
                     np.linalg.norm(Sw @ nc),
                     abs(np.trace(Sw) - np.trace(S)),
                     abs(0.5 * (np.trace(Sw) ** 2 - np.trace(Sw @ Sw)) - np.linalg.det(S))))

    a = np.array(rows)
    labels = ["|dd|   frame-free vs framed",
              "|dimg| frame-free vs framed",
              "|dd|   MISMATCHED frame",
              "|dimg| MISMATCHED frame",
              "|d(tr S)|  under rotation",
              "|d(det S)| under rotation",
              "|Sw n_c|",
              "tr(Sw) - 2H",
              "0.5[(tr Sw)^2-tr(Sw^2)] - K"]

    print(f"    {'quantity':<32}{'median':>12}{'p95':>12}{'max':>12}   units")
    for i, lab in enumerate(labels):
        unit = "dx" if i < 4 else "1/dx or 1/dx^2"
        print(f"    {lab:<32}{np.median(a[:, i]):>12.3e}{np.percentile(a[:, i], 95):>12.3e}"
              f"{a[:, i].max():>12.3e}   {unit}")

    print("\n### Read: rows 1-2 say the frame-free form is EXACT, not an approximation.")
    print("### Rows 3-4 say a frame mismatch moves the image by a median of %.2f dx and up"
          % np.median(a[:, 3]))
    print("### to %.2f dx, while rows 5-6 say the Jacobian invariants -- and therefore J and"
          % a[:, 3].max())
    print("### the 6.1 mass check -- do not move at all. That is the silent failure.")


if __name__ == "__main__":
    run()
