#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Spline1DMapping.H"

#include <iostream>
#include "parstream.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

// Static for the maximum polynomial term
const int Spline1DMapping::s_order = 3;

Spline1DMapping::Spline1DMapping() 
{
// Set reasonable defaults
  m_tol = 1e-12;
  m_ptrPoints = NULL;
}


Spline1DMapping::~Spline1DMapping() 
{
  delete m_ptrPoints;
}


Spline1DMapping::Spline1DMapping(Vector<Real>& a_points)
{
  buildFromPoints(a_points);
}


// Read from ParmParse with tag "splineMapping.xxx"
void Spline1DMapping::readParmParseInputs()
{
  // Here's our parser.
  ParmParse pp("Spline1DMapping");

  // This determines the amount of diagnositic output generated
  pp.query("verbosity", m_verbosity);
  CH_assert(m_verbosity >= 0);

  // Number of grid points to read in
  int numPoints;
  pp.get("num_points", numPoints);

  // Interface locations between levels
  Vector<Real> points;
  int found = pp.queryarr("points", points, 0, numPoints);
  // If the points aren't provided, then use evenly spaced
  if (found == 0)
  {
    MayDay::Warning("Spline1DMapping using evenly-spaced vertical mapping!!");
    Real dxi = 1.0 / (Real) (numPoints - 1);
    for (int i=0; i < numPoints; i++)
      points[i] = dxi * (Real) i;
  }

  // Build the spline coefs
  buildFromPoints(points);
}


void Spline1DMapping::print()
{
  if (m_verbosity == 0) return;

  pout() << "Spline1DMapping: " << endl;
  pout() << "  verbosity = " << m_verbosity << endl;
  
  pout() << "  num_points = " << m_numPoints << endl;
  pout() << "  m_points = ";
  for (int i = 0; i < m_points.size(); ++i)
  {
    pout() << m_points[i] << " ";
  }
  pout() << endl;

  if (m_verbosity < 3) return;

  pout() << "  Cubic fit by interval = " << endl;
  for (int i = 0; i < m_points.size()-1; ++i)
  {
    pout() << "  " << i << ": " << 
      m_points[i] << " " <<
      m_deriv1[i] << " " <<
      m_deriv2[i] << " " <<
      m_deriv3[i] << " " << endl;
  }
}


const Real* Spline1DMapping::getPoints()
{
  CH_assert(m_ptrPoints);
  return m_ptrPoints;
}


int Spline1DMapping::getNumPoints()
{
  return m_numPoints;
}


void Spline1DMapping::buildFromPoints(Vector<Real>& a_newPoints)
{
  CH_assert(a_newPoints.size() > 4);

  // Make sure they're in order and separated
  for (int i = 0; i < (a_newPoints.size()-1); ++i)
    CH_assert((a_newPoints[i+1] - a_newPoints[i]) > m_tol);

  // Save the new data
  m_numPoints = a_newPoints.size();
  m_points = a_newPoints;

  // Allocate the derivative arrays - these live on the intervals
  m_deriv1 = Vector<Real>(m_numPoints,0);
  m_deriv2 = Vector<Real>(m_numPoints,0);
  m_deriv3 = Vector<Real>(m_numPoints,0);

  buildSpline();

  // Copy over into real array of points
  delete m_ptrPoints;
  m_ptrPoints = new Real[m_numPoints];
  for (int i = 0; i < m_numPoints; ++i)
    m_ptrPoints[i] = m_points[i];
}


/**
 *  Reinterpolates this mapping to a new size, using the existing 
 *  spline to interpolate new values.
 */
void Spline1DMapping::reinterpolate(int a_newNumPoints)
{
  CH_assert(m_numPoints > 0);
  CH_assert((a_newNumPoints > 0));
  // Calculate the new number of points
  Vector<Real> newPoints(a_newNumPoints,0);

  // First and last points must coincide
  newPoints[0] = m_points[0];
  newPoints[a_newNumPoints-1] = m_points[m_numPoints-1];

  // Loop through the new points, evaluate from the old spline
  Real h = 1.0 / (Real) (m_numPoints-1);
  Real newh = 1.0 / (Real) (a_newNumPoints-1);
  int j = 0;
  for (int i = 1; i < a_newNumPoints-1; ++i)
  {
    // Calculate the new mapped coordinate location on [0,1]
    Real x = newh * (Real) i;
    // Find the corresponding interval in the old spline
    for (; x > h*(Real) j; ++j);
    --j;
    // Calculate the new location
    Real dx = x - h*(Real) j;
    newPoints[i] = m_points[j] + 
      dx*(m_deriv1[j] + dx*(m_deriv2[j] + dx*m_deriv3[j]));

    //pout() << "newPoints[ " << i << "] = " << newPoints[i] <<
    //  " from old interval: " << j << endl;
  }

  // Reallocate memory and rebuild the spline
  buildFromPoints(newPoints);
}

 
// Code adapted from: triadiagonal
// http://www.pcs.cnu.edu/~bbradie/cpp/interp.C
void Spline1DMapping::tridiagonal(
  Vector<Real>& a_lower, Vector<Real>& a_diag, Vector<Real>& a_upper,
  Vector<Real>& a_rhs)
{
  int n = a_diag.size();
  for (int i = 0; i < n-1; ++i)
  {
     a_upper[i] /= a_diag[i];
     a_diag[i+1] -= a_lower[i]*a_upper[i];
  }

  a_rhs[0] /= a_diag[0];
  for (int i = 1; i < n; i++ )
     a_rhs[i] = ( a_rhs[i] - a_lower[i-1] * a_rhs[i-1] ) / a_diag[i];
     
  for (int i = n-2; i >= 0; i-- )
     a_rhs[i] -= a_rhs[i+1] * a_upper[i];
}
                    

// Code adapted from: cubic_nak
// http://www.pcs.cnu.edu/~bbradie/cpp/interp.C
void Spline1DMapping::buildSpline()
{
  int n = m_numPoints;
  Real h = 1.0 / (Real) (n-1);
  Vector<Real> dd(m_numPoints);
  Vector<Real> dl(m_numPoints, h);
  Vector<Real> du(m_numPoints, h);
 
  Vector<Real>& b = m_deriv1;
  Vector<Real>& c = m_deriv2;
  Vector<Real>& d = m_deriv3;
  Vector<Real>& f = m_points;

  for (int i = 0; i < n-2; i++) 
  {
      dd[i] = 4.0 * h;
      c[i]  = ( 3.0 / h ) * ( f[i+2] - 2.0*f[i+1] + f[i]);
  }
  dd[0] += 2.0*h;
  dd[n-3] += 2.0*h;
  du[0] -= h;
  dl[n-4] -= h;
 
  tridiagonal(dl, dd, du, c );
 
  for (int i = n-3; i >= 0; i--)
      c[i+1] = c[i];
  c[0] = 2.0 * c[1] - c[2];
  c[n-1] = 2.0 * c[n-2] - c[n-3];
  for (int i = 0; i < n-1; i++)
  {
      d[i] = (c[i+1] - c[i]) / (3.0 * h);
      b[i] = (f[i+1] - f[i]) / h - h*(c[i+1] + 2.0*c[i])/3.0;
  }
}

//
// Evaluate the spline to calculate value at a point
Real Spline1DMapping::getValue(const Real a_xi)
{
  CH_assert((a_xi >= 0.0 - m_tol) && (a_xi <= 1.0+m_tol));
  int n = m_numPoints-1;
  Real dxi = 1.0 / (Real) n;

  // Calculate the index of the xi value
  int i = floor((a_xi + m_tol)/ dxi);
  i = Max(i, 0); // Make sure you don't get a negative index
  CH_assert((i >= 0) && (i <= n));

  // Calculate the value location from the beginning of the interval
  Real ixi = a_xi - dxi*(Real) i;
  Real val = m_points[i] + 
        ixi*(m_deriv1[i] + ixi*(m_deriv2[i] + ixi*m_deriv3[i]));

  return val;
}


// Evaluate the spline to calculate the mapping derivative at a point
Real Spline1DMapping::getDerivative(const Real a_xi)
{
  CH_assert((a_xi >= 0.0 - m_tol) && (a_xi <= 1.0+m_tol));
  int n = m_numPoints-1;
  Real dxi = 1.0 / (Real) n;

  // Calculate the index of the xi value
  int i = floor((a_xi - m_tol)/ dxi);
  i = Max(i, 0); // Make sure you don't get a negative index
  CH_assert((i >= 0) && (i < n));

  // Calculate the value location from the beginning of the interval
  Real ixi = a_xi - dxi*(Real) i;
  Real deriv = m_deriv1[i] + ixi*(2.0*m_deriv2[i] + 3.0*m_deriv3[i]*ixi);

  return deriv;
}

#include "NamespaceFooter.H"
