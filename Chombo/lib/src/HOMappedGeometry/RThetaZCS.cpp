#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NewCoordSys.H"
#include "RThetaZCS.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

/// constructor
RThetaZCS::RThetaZCS(const RealVect& a_dx, const RealVect a_stretch,
                     const RealVect& a_origin)
  : m_stretch(a_stretch), m_origin(a_origin)
{
  m_Pi = 4.0*atan(1.0);
  m_dx = a_dx;
}
/**
   Destructor.
*/
RThetaZCS::~RThetaZCS()
{
}

/// given coordinate in mapped space, return its location in real space
RealVect
RThetaZCS::realCoord(const RealVect& a_Xi) const
{
  RealVect localXi = a_Xi;
  localXi *= m_stretch;
  localXi += m_origin;

  RealVect xyzLoc;
  D_TERM6(xyzLoc[0] = localXi[0]*cos(localXi[1]);,
          xyzLoc[1] = localXi[0]*sin(localXi[1]);,
          xyzLoc[2] = localXi[2];,
          xyzLoc[3] = localXi[3];,
          xyzLoc[4] = localXi[4];,
          xyzLoc[5] = localXi[5];)
  return xyzLoc;
}

/// given coordinate in real space, return its location in the mapped space
RealVect RThetaZCS::mappedCoord(const RealVect& a_x) const
{
  RealVect mappedXi;

  D_TERM6(mappedXi[0] = sqrt(a_x[0]*a_x[0] + a_x[1]*a_x[1]);,
          mappedXi[1] = acos(a_x[0]/mappedXi[0]);
          if (a_x[1] < 0) mappedXi[1] = 2.0*m_Pi-mappedXi[1];,
          mappedXi[2] = a_x[2];,
          mappedXi[3] = a_x[3];,
          mappedXi[4] = a_x[4];,
          mappedXi[5] = a_x[5];)

    return (mappedXi-m_origin)/m_stretch;

}

// Calculate the derivative of each coordinate vector
Real RThetaZCS::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  RealVect localXi = a_Xi;
  if (localXi[0] < 0.)
    { // negative radius: make radius positive, and flip the angle
      localXi[0] = -localXi[0];
      localXi[1] += Pi;
    }
  Real value = 0.0;
  // this is for test purposes
  //RealVect rThetaZ = mappedCoord(a_X);
  RealVect X = realCoord(localXi);
  Real testValue=value;
  Real difference=0;
  Real eps = 1.0e-7;
  if (a_dirX == 0)
    {
      if (a_dirXi == 0)
        {
          // dx/dr = cos(theta) = x/r
          testValue = X[0]/sqrt(X[0]*X[0] + X[1]*X[1]);

          value = m_stretch[0]*cos(localXi[1]);
          if (localXi[0] < 0.)
            {
              value *= -1.;
            }
        }
      else if (a_dirXi == 1)
        {
          // dx/dtheta = -rsin(theta) = -y
          testValue = -X[1];
          value = -m_stretch[1]*localXi[0]*sin(localXi[1]);
        }
      // else value stays 0

    }
  else if (a_dirX == 1)
    {
      if (a_dirXi == 0)
        {
          // dy/dr = sin(theta) = y/r
          testValue = X[1]/sqrt(X[0]*X[0] + X[1]*X[1]);
          value = m_stretch[0]*sin(localXi[1]);
          if (localXi[0] < 0.)
            {
              value *= -1.;
            }
        }
      else if (a_dirXi == 1)
        {
          // dy/dtheta = rcos(theta) = x
          testValue = X[0];
          value = m_stretch[1]*localXi[0]*cos(localXi[1]);
        }
      // else value stays 0
    }
  else if (a_dirX == 2)
    {
      if (a_dirXi == 2)
        {
          value = 1.0;
          value *= m_stretch[a_dirXi];
          testValue = value;
        }

    }
  else
    {
      MayDay::Error("Bad dirX in RZThetaCoordSys::dXdXi");
    }

  difference = value - testValue;
  if (Abs(difference) > eps)
    {
      MayDay::Warning("possible bad value for dX/dXi");
    }

  return value;

}


// -- begin factory implementations ---------------------------

RThetaZCSFactory::RThetaZCSFactory(const RealVect& a_stretch, 
                                   const RealVect& a_origin)
{
  m_stretch = a_stretch;
  m_origin = a_origin;
}


NewCoordSys*
RThetaZCSFactory::getCoordSys(const ProblemDomain& a_domain,
                              const RealVect& a_dx) const
{
  RThetaZCS* newCSPtr = new RThetaZCS(a_dx, m_stretch, m_origin);
  return static_cast< NewCoordSys* >(newCSPtr);
}


#include "NamespaceFooter.H"
