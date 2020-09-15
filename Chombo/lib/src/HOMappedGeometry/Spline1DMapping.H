#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SPLINE1DMAPPING_H_
#define _SPLINE1DMAPPING_H_

#include "REAL.H"
#include "Vector.H"

#include "NamespaceHeader.H"


// Class to calculate a vertical mapping with splines
// Point data read from ParmParse initially, parameterized on [0,1]
class Spline1DMapping
{
  private:
    Real m_tol;
    int m_verbosity;
    int m_numPoints;
    static const int s_order;
    Real* m_ptrPoints;
    Vector<Real> m_points;
    Vector<Real> m_deriv1;
    Vector<Real> m_deriv2;
    Vector<Real> m_deriv3;

  public:
    Spline1DMapping();
    
    ~Spline1DMapping();

    Spline1DMapping(Vector<Real>& a_points);

    // Read from ParmParse with tag "splineMapping.xxx"
    void readParmParseInputs();

    void print();

    const Real* getPoints();

    int getNumPoints();

    /**
     *  Using the spline to interpolate at \xi \in [0,1]
     */
    Real getValue(const Real a_xi);

    /**
     *  Using the spline to interpolate a derivative at \xi \in [0,1]
     */
    Real getDerivative(const Real a_xi);

    void buildFromPoints(Vector<Real>& a_newPoints);

    /**
     *  Reinterpolates this mapping to a new size, using the existing 
     *  spline to interpolate new values.
     */
    void reinterpolate(int a_newNumPoints);

  protected:
 
    // Code adapted from: triadiagonal
    // http://www.pcs.cnu.edu/~bbradie/cpp/interp.C
    void tridiagonal(Vector<Real>& a_lower, Vector<Real>& a_diag, Vector<Real>& a_upper,
        Vector<Real>& a_rhs);

    // Code adapted from: cubic_nak
    // http://www.pcs.cnu.edu/~bbradie/cpp/interp.C
    void buildSpline();

};

#include "NamespaceFooter.H"

#endif
