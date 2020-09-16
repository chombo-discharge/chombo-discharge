#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GaussianQuadrature.H"

#include "NamespaceHeader.H"

GaussianQuadrature::GaussianQuadrature(int a_numberOfPoints)
{
  m_coefficients.resize(a_numberOfPoints);
  // initialize locations to zero
  for (int i=0; i<m_coefficients.size(); i++)
    {
      m_coefficients[i].location = RealVect::Zero;
    }

  // this needs to be dimension-dependent
  if (SpaceDim < 3)
    {
      // in 1D or 2D, just make this a single point.
      m_coefficients.resize(1);
      m_coefficients[0].location = RealVect::Zero;
      m_coefficients[0].weight = 1.0;
    }
  else if (SpaceDim == 3)
    {
      if (a_numberOfPoints == 1)
        {
          // simplest possible case (midpoint rule)
          m_coefficients[0].location[0] = 0.0;
          m_coefficients[0].weight = 2.0;
        }
      else if (a_numberOfPoints == 2)
        {
          m_coefficients[0].location[0] =  -sqrt(1.0/3.0);
          m_coefficients[0].weight = 1.0;

          m_coefficients[1].location[0] = sqrt(1.0/3.0);
          m_coefficients[1].weight = 1.0;
        }
      else if (a_numberOfPoints == 3)
        {
          m_coefficients[0].location[0] =  -sqrt(0.6);
          m_coefficients[0].weight = 5.0/9.0;

          m_coefficients[1].location[0] = 0.0;
          m_coefficients[1].weight = 8.0/9.0;

          m_coefficients[2].location[0] = + sqrt(0.6);
          m_coefficients[2].weight = 5.0/9.0;
        }
      else if (a_numberOfPoints == 4)
        {
          // we might want to re-examine these
          m_coefficients[0].location[0] = -0.861136312;
          m_coefficients[0].weight = 0.347854845;

          m_coefficients[1].location[0] = -0.339981044;
          m_coefficients[1].weight = 0.652145155;

          m_coefficients[2].location[0] = 0.339981044;
          m_coefficients[2].weight = 0.652145155;

          m_coefficients[3].location[0] = 0.861136312;
          m_coefficients[3].weight = 0.347854845;

        }
      else if (a_numberOfPoints == 5)
        {

          m_coefficients[0].location[0] = 0.0;
          m_coefficients[0].weight = 0.56888889;

          m_coefficients[1].location[0] = -0.53846931;
          m_coefficients[1].weight = 0.47862867;

          m_coefficients[2].location[0] = 0.53846931;
          m_coefficients[2].weight = 0.47862867;

          m_coefficients[3].location[0] = -0.90617985;
          m_coefficients[3].weight = 0.23692689;

          m_coefficients[4].location[0] = 0.90617985;
          m_coefficients[4].weight = 0.23692689;

        }
      else
        {
          MayDay::Error("undefined number of points for gaussian quadrature");
        }
    }
  else
    {
      // undefined SpaceDim
      MayDay::Error("GaussianQuadrature undefined for CH_SPACEDIM Dimensions");
    }
}

/// returns quadrature coefficients
const Vector<QuadratureElement >
GaussianQuadrature::coefficients(int dir1,
                                 int dir2) const
{
  Vector<QuadratureElement> returnVal(m_coefficients.size());

  // copy m_quadrature into returnVal, adjusting components
  // to account for dir1 and dir2
  // what we want is for the quadratur locations to be zero in
  // dir1 and dir2, since we only want to vary locations
  // in directions which are neither dir1 or dir2.
  // returnDir[i] specifies which coefficient components
  // goes into the return location component i.
  IntVect returnDir;
  returnDir[dir2] = SpaceDim-1;
  if (SpaceDim > 1) returnDir[dir1] = SpaceDim-2;
  int tempIndex = 0;
  for (int comp=0; comp<SpaceDim; comp++)
    {
      if (comp != dir1 && comp != dir2)
        {
          returnDir[comp] = tempIndex;
          tempIndex++;
        }
    }

  for (int i=0; i<m_coefficients.size(); i++)
    {
      for (int dir=0; dir<SpaceDim; dir++)
        {
          returnVal[i].location[dir]  = m_coefficients[i].location[returnDir[dir]];
        }
      returnVal[i].weight = m_coefficients[i].weight;
    }

  return returnVal;

}


/* this is essentially the length (area, volume, etc) of an edge
   given interval length a_length and dir1, dir2 the same as in the
   coefficients function
*/
Real
GaussianQuadrature::weightMult(const RealVect& a_length, int dir1, int dir2) const
{
  Real mult = 1.0;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (dir != dir1 && dir != dir2)
        {
          mult *= a_length[dir]/2.0;
        }
    }
  return mult;
}


/// make a copy of myself
EdgeQuadrature*
GaussianQuadrature::new_quadrature() const
{
  GaussianQuadrature* new_quad = new GaussianQuadrature(m_coefficients.size());
  return static_cast<EdgeQuadrature*>(new_quad);
}

#include "NamespaceFooter.H"
