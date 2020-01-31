#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "RThetaZCoordSys.H"

#include "UsingNamespace.H"

RThetaZCoordSys::RThetaZCoordSys()
{
  m_Pi = 4.0*atan(1.0);
  m_isDefined = false;
}

RThetaZCoordSys::~RThetaZCoordSys()
{
  if (m_dNdXi.size() > 0)
    {
      for (int i=0; i<m_dNdXi.size(); i++)
        {
          if (m_dNdXi[i] != NULL)
            {
              delete m_dNdXi[i];
              m_dNdXi[i] = NULL;
            }
        }
    }
}

void
RThetaZCoordSys::define(const DisjointBoxLayout& a_grids,
                        const ProblemDomain& a_domain,
                        const RealVect& a_cellSpacing,
                        const IntVect& a_ghostVect)
{
  FourthOrderCoordSys::define(a_grids, a_domain,
                              a_cellSpacing, a_ghostVect);


#if 0
  // not sure what I need this for
  if (m_dNdXi.size() > 0)
    {
      for (int i=0; i<m_dNdXi.size(); i++)
        {
          if (m_dNdXi[i] != NULL)
            {
              delete m_dNdXi[i];
              m_dNdXi[i] = NULL;
            }
        }
    }

  m_dNdXi.resize(SpaceDim, NULL);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      m_dNdXi[dir] = new LevelData<FluxBox>(a_grids,
                                            SpaceDim, a_ghostVect);
    }
#endif

  m_isDefined = true;
}



void
RThetaZCoordSys::regrid(const DisjointBoxLayout& a_newGrids)
{
  // can't call regrid unless we've already "gridded"
  CH_assert(m_isDefined);

  // for now, just re-do define. May want to make this more
  // efficient in the future
  define(a_newGrids, m_domain, m_dx, m_ghostVect);

}

RealVect
RThetaZCoordSys::realCoord(const RealVect& a_X) const
{

  RealVect xyzLoc;
  D_TERM(xyzLoc[0] = a_X[0]*cos(a_X[1]);,
         xyzLoc[1] = a_X[0]*sin(a_X[1]);,
         xyzLoc[2] = a_X[2];)

    return xyzLoc;
}

RealVect
RThetaZCoordSys::mappedCoord(const RealVect& a_x) const
{
  RealVect mappedXi;

  D_TERM(mappedXi[0] = sqrt(a_x[0]*a_x[0] + a_x[1]*a_x[1]);,
         mappedXi[1] = acos(a_x[0]/mappedXi[0]);
         if (a_x[1] < 0) mappedXi[1] = 2.0*m_Pi-mappedXi[1];,
         mappedXi[2] = a_x[2];)

  return m_stretch*mappedXi;

}


Real
RThetaZCoordSys::pointwiseJ(const RealVect& a_X) const
{
   return FourthOrderCoordSys::pointwiseJ( a_X );
}


void
RThetaZCoordSys::mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                                      LevelData<FluxBox>& a_F)
{
  FourthOrderCoordSys::mappedGridDivergence(a_divF, a_F);
}



Real
RThetaZCoordSys::dXdXi(const RealVect& a_X, int a_dirX, int a_dirXi) const
{
  Real value = 0.0;
  // this is for test purposes
  RealVect rThetaZ = mappedCoord(a_X);
  Real testValue=value;
  Real difference=0;
  Real eps = 1.0e-7;
  if (a_dirX == 0)
    {
      if (a_dirXi == 0)
        {
          // dx/dr = cos(theta) = x/r
          value = a_X[0]/sqrt(a_X[0]*a_X[0] + a_X[1]*a_X[1]);

          testValue = cos(rThetaZ[1]);

        }
      else if (a_dirXi == 1)
        {
          // dx/dtheta = -rsin(theta) = -y
          value = -a_X[1];
          testValue = -rThetaZ[0]*sin(rThetaZ[1]);
        }
      // else value stays 0

    }
  else if (a_dirX == 1)
    {
      if (a_dirXi == 0)
        {
          // dy/dr = sin(theta) = y/r
          value = a_X[1]/sqrt(a_X[0]*a_X[0] + a_X[1]*a_X[1]);
          testValue = sin(rThetaZ[1]);
        }
      else if (a_dirXi == 1)
        {
          // dy/dtheta = rcos(theta) = x
          value = a_X[0];
          testValue = rThetaZ[0]*cos(rThetaZ[1]);
        }
      // else value stays 0
    }
  else if (a_dirX == 2)
    {
      if (a_dirXi == 2)
        {
          value = 1.0;
          testValue = value;
        }

    }
  else
    {
      MayDay::Error("Bad dirX in RZThetaCoordSys::dXdXi");
    }

  difference = value - testValue;
  if (difference > eps)
    {
      MayDay::Warning("possible bad value for dX/dXi");
    }

  return value;
}



// -- begin factory implementations ---------------------------

// default volInterval is all directions
RThetaZCoordSysFactory::RThetaZCoordSysFactory(const ProblemDomain& a_baseDomain,
                                               const Vector<int>& a_vectRefRatios,
                                               const RealVect& a_baseDx,
                                               const RealVect& a_stretch)
  : m_volInterval(0,SpaceDim-1)
{
  int numLevels = a_vectRefRatios.size() +1;

  m_stretch = a_stretch;

  // single-level case (if a single invalid refRatio is passed in)
  if (numLevels == 2 && a_vectRefRatios[0] <= 0)
    {
      numLevels = 1;
    }

  m_vectDomains.resize(numLevels);
  m_dxVect.resize(numLevels);
  m_vectRefRatios = a_vectRefRatios;

  m_vectDomains[0] = a_baseDomain;
  m_dxVect[0] = a_baseDx;

  for (int lev=1; lev<numLevels; lev++)
    {
      m_vectDomains[lev] = m_vectDomains[lev-1];
      m_vectDomains[lev].refine(m_vectRefRatios[lev-1]);
      m_dxVect[lev] = m_dxVect[lev-1];
      m_dxVect[lev] /= m_vectRefRatios[lev-1];
    }

}



CoordSys<FArrayBox, FluxBox>*
RThetaZCoordSysFactory::getCoordSys(const DisjointBoxLayout& a_grids,
                                    const ProblemDomain& a_levelDomain,
                                    const IntVect& a_ghostVect) const
{
  // first identify level
  int level = -1;
  for (int lev=0; lev<m_vectDomains.size(); lev++)
    {
      if (m_vectDomains[lev] == a_levelDomain)
        {
          level = lev;
        }
    }

  if (level < 0)
    {
      pout() << "attempted to match bad domain: " << a_levelDomain
             << endl;
      MayDay::Error("Invalid level in RThetaZCoordSysFactory::getCoordSys");
    }

  RThetaZCoordSys* newCSPtr = new RThetaZCoordSys;

  newCSPtr->m_volInterval = m_volInterval;
  newCSPtr->stretch(m_stretch);

  newCSPtr->define(a_grids, m_vectDomains[level],
                   m_dxVect[level], a_ghostVect);


  return static_cast< CoordSys<FArrayBox,FluxBox>* >(newCSPtr);

}

