/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PolyUtils.cpp
  @brief  Implementation of CD_PolyUtils.cpp
  @author Robert Marskar
*/

// Chombo includes
#include <PolyGeom.H>

// Our includes
#include <CD_PolyUtils.H>
#include <CD_NamespaceHeader.H>

RealVect
PolyUtils::brentRootFinder(const RefCountedPtr<BaseIF>& a_impFunc, const RealVect& a_point1, const RealVect& a_point2)
{
  const Real tol = PolyGeom::getTolerance();

  const unsigned int maxIter = 100;
#if defined(CH_USE_DOUBLE)
  const Real eps = 3.0e-15;
#elif defined(CH_USE_FLOAT)
  const Real eps = 3.0e-7;
#else
#error Unknown Chombo precision
#endif

  RealVect aPt, bPt, cPt, unitTangent;
  Real     fa, fb, fc;
  Real     a, b, c, d, e;
  Real     p, q, r, s;
  Real     tol1, xm;

  // Translate the problem onto [RealVect::Zero, a_point2 - a_point1] so that we can solve for only the distance ( > 0.)
  // All function evaluations should use these coordinates
  aPt         = RealVect::Zero;
  bPt         = a_point2 - a_point1;
  unitTangent = bPt / bPt.vectorLength();

  // Starting points and values
  fa = -a_impFunc->value(aPt + a_point1);
  fb = -a_impFunc->value(bPt + a_point1);

  // Check for stupid input
  if (fa * fb > 0) {
    pout() << "fa " << fa << " fb " << fb << endl;
    MayDay::Abort("brentRootFinder. Root must be bracketed, but instead the supplied end points have the same sign.");
  }

  // Init to be safe
  a = 0.;
  b = 0.;
  c = 0.;
  d = 0.;
  e = 0.;

  cPt = bPt;
  c   = b;
  fc  = fb;

  for (int i = 0; i < maxIter; i++) {

    // If this triggers, the root is not on [b,c] and we can adjust the bounding interval
    if (fb * fc > 0) {
      cPt = aPt;
      c   = aPt.vectorLength();
      fc  = fa;
      d   = b - a;
      e   = d;
    }

    if (Abs(fc) < Abs(fb)) {
      aPt = bPt;
      a   = aPt.vectorLength();
      bPt = cPt;
      b   = bPt.vectorLength();
      cPt = aPt;
      c   = cPt.vectorLength();
      fa  = fb;
      fb  = fc;
      fc  = fa;
    }

    // Convergence check
    tol1 = 2.0 * eps * abs(bPt.vectorLength()) + 0.5 * tol;
    xm   = 0.5 * (c - bPt.vectorLength());
    if (Abs(xm) <= tol1 || fb == 0.0) {
      break;
    }

    if (Abs(e) >= tol1 && Abs(fa) > Abs(fb)) {
      // Attempt inverse quadratic interpolation
      s = fb / fa;
      if (aPt.vectorLength() == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      }
      else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1.));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }

      if (p > 0)
        q = -q;

      p = Abs(p);

      if (2.0 * p < std::min((Real)3.0 * xm * q - std::abs(tol1 * q), std::abs(e * q))) {
        // Accept interpolation
        e = d;
        d = p / q;
      }
      else {
        // Interpolation failed, use bisection
        d = xm;
        e = d;
      }
    }
    else {
      // Bounds decreasing too slowly, use bisection
      d = xm;
      e = d;
    }

    // Move last best guess to a
    aPt = bPt;
    a   = aPt.vectorLength();
    fa  = fb;

    // Evaluate new trial root
    if (Abs(d) > tol1) {
      bPt = bPt + d * unitTangent;
      b   = bPt.vectorLength();
    }
    else {
      if (xm < 0) {
        bPt = bPt - tol1 * unitTangent;
      }
      else {
        bPt = bPt + tol1 * unitTangent;
      }
      b = bPt.vectorLength();
    }

    fb = -a_impFunc->value(bPt + a_point1);
    //    std::cout << "Iteration: " << i << "\t Function value: " << fb << std::endl;
  }

  // if(i >= maxIter){
  //   cerr << "IntersectionUtils::BrentRootFinder: exceeding maximum iterations: " << maxIter << endl;
  // }

  return bPt + a_point1;
}

#include <CD_NamespaceFooter.H>
