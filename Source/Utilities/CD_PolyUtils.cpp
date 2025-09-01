/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PolyUtils.cpp
  @brief  Implementation of CD_PolyUtils.cpp
  @author Robert Marskar
*/

// Std includes
#include <cmath>
#include <limits>

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

Real
PolyUtils::brentSolve(const Real a_point1, const Real a_point2, const std::function<Real(const Real x)> a_func) noexcept
{
  CH_TIME("PolyUtils::bretSolve");

  constexpr Real tol     = 1E-10;
  constexpr Real eps     = std::numeric_limits<Real>::epsilon();
  constexpr int  maxIter = 500;

  Real root = std::numeric_limits<Real>::max();

  Real fa = a_func(a_point1);
  Real fb = a_func(a_point2);

  if (fa * fb > 0.0) {
    MayDay::Warning("PolyUtils::brentSolve -- abort because alpha is > 0 everywhere");

    return root;
  }

  Real a  = a_point1;
  Real b  = a_point2;
  Real c  = a;
  Real fc = fa;
  Real d  = b - a;
  Real e  = d;

  Real x  = b;
  Real fx = fb;

  for (int iter = 0; iter < maxIter; iter++) {
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c  = a;
      fc = fa;
      d  = b - a;
      e  = d;
    }

    if (std::abs(fc) < std::abs(fb)) {
      a = b;
      b = c;
      c = a;

      fa = fb;
      fb = fc;
      fc = fa;
    }

    const Real tol1 = 2 * eps * std::abs(b) + 0.5 * tol;
    const Real m    = 0.5 * (c - b);

    if (std::abs(m) <= tol1 || fb == 0.0) {
      return b;
    }

    if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
      const Real s = fb / fa;

      Real p;
      Real q;

      if (a == c) {
        p = 2 * m * s;
        q = 1.0 - s;
      }
      else {
        const Real q1 = fa / fc;
        const Real r  = fb / fc;

        p = s * (2 * m * q1 * (q1 - r) - (b - a) * (r - 1));
        q = (q1 - 1) * (r - 1) * (s - 1);
      }

      if (p > 0) {
        q = -q;
      }

      p = std::abs(p);

      if (2 * p < std::min(3 * m * q - std::abs(tol1 * q), std::abs(e * q))) {
        e = d;
        d = p / q;
      }
      else {
        d = m;
        e = d;
      }
    }
    else {
      d = m;
      e = d;
    }

    a  = b;
    fa = fb;

    if (std::abs(d) > tol1) {
      b += d;
    }
    else {
      b += (m > 0) ? tol1 : -tol1;
    }

    fb = a_func(b);
  }

  return root;
}

#include <CD_NamespaceFooter.H>
