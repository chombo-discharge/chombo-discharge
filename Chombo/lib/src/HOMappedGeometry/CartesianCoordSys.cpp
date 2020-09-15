#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CartesianCoordSys.H"

#include "NamespaceHeader.H"

CartesianCoordSys::~CartesianCoordSys()
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
CartesianCoordSys::define(const DisjointBoxLayout& a_grids,
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
CartesianCoordSys::regrid(const DisjointBoxLayout& a_newGrids)
{
  // can't call regrid unless we've already "gridded"
  CH_assert(m_isDefined);

  // for now, just re-do define. May want to make this more
  // efficient in the future
  define(a_newGrids, m_domain, m_dx, m_ghostVect);

}

RealVect
CartesianCoordSys::realCoord(const RealVect& a_X) const
{
  // this is the simplest mapping ever
   return m_stretch*a_X;
}

RealVect
CartesianCoordSys::mappedCoord(const RealVect& a_x) const
{
  // simplest possible mapping
  return a_x/m_stretch;
}


Real
CartesianCoordSys::pointwiseJ(const RealVect& a_X) const
{
   return FourthOrderCoordSys::pointwiseJ( a_X );
}


void
CartesianCoordSys::mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                                        const LevelData<FluxBox>& a_F) const
{
  FourthOrderCoordSys::mappedGridDivergence(a_divF, a_F);
}



Real
CartesianCoordSys::dXdXi(const RealVect& a_X, int a_dirX, int a_dirXi) const
{
  // simplest possible mapping
  Real value = 0.0;
  if (a_dirX == a_dirXi) value = m_stretch[a_dirXi];
  return value;
}



// -- begin factory implementations ---------------------------

CartesianCoordSysFactory::CartesianCoordSysFactory(const ProblemDomain& a_baseDomain,
                                                   const Vector<int>& a_vectRefRatios,
                                                   const RealVect& a_baseDx,
                                                   const RealVect& a_stretch)
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
CartesianCoordSysFactory::getCoordSys(const DisjointBoxLayout& a_grids,
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
      MayDay::Error("Invalid level in CartesianCoordSysFactory::getCoordSys");
    }

  CartesianCoordSys* newCSPtr = new CartesianCoordSys;

  newCSPtr->stretch(m_stretch);

  newCSPtr->define(a_grids, m_vectDomains[level],
                   m_dxVect[level], a_ghostVect);


  return static_cast< CoordSys<FArrayBox,FluxBox>* >(newCSPtr);

}

#include "NamespaceFooter.H"
