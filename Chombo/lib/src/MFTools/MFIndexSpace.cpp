#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MFIndexSpace.H"
#include "NamespaceHeader.H"

MFIndexSpace::MFIndexSpace()
{
}

MFIndexSpace::~MFIndexSpace()
{
  for (int i=0; i<m_ebis.size(); i++)
  {
    delete m_ebis[i];
  }
}

void
MFIndexSpace::define(const Box                     & a_domain,
                     const RealVect                & a_origin,
                     const Real                    & a_dx,
                     const Vector<GeometryService*>& a_geoservers,
                     int                             a_nCellMax,
                     int                             a_maxCoarsenings,
                     bool                            a_fixOnlyFirstPhaseRegNextToMultiValued)
{
  for (int i=0; i<m_ebis.size(); i++)
  {
    delete m_ebis[i];
  }
  int m_nCellMax = a_nCellMax;
  if (a_nCellMax <= 0)
    {
      if (SpaceDim == 2)
        {
          m_nCellMax = 64;
        }
      else
        {
          m_nCellMax = 16;
        }
    }
  //hack because breaking up boxes seems to be broken
  // EBIndexSpace::s_MFSingleBox = true;
  // m_nCellMax = 1024;

  int phases = a_geoservers.size();

  m_ebis.resize(phases);
  Vector<EBISLevel*> elevels(phases, NULL);
  Vector<EBISLevel*> elevelsFine(phases, NULL);
  Vector<EBISLevel*> elevOne;
  Vector<EBISLevel*> elevTwo;

  if (phases != 2)
    {
      pout() << "MFIndexSpace is hardwired for two fluids in many places." << endl;
      pout() << "This error message is to force this assumption to be true." << endl;
      MayDay::Error("What, you don't like two phase flow?");
    }
  int nlevels = 0;
  bool deeper = true;

  while (deeper)
  {
    for (int i=0; i<phases; i++)
      {
        //if a_fixOnlyFirstPhaseRegNextToMultiValued == true, this logic only fixes regularNextToMultiValued on phase[0]
        bool fixRegularNextToMultiValued = true;

        if (a_fixOnlyFirstPhaseRegNextToMultiValued)
          {
            if (i != 0)
              {
                fixRegularNextToMultiValued = false;
              }
          }

        if (nlevels == 0)
        {
          m_ebis[i] = new EBIndexSpace();
          m_ebis[i]->setCellMax(m_nCellMax);
          CH_assert(a_geoservers[i] != NULL);

          elevels[i] = m_ebis[i]->buildFirstLevel(a_domain,a_origin,a_dx,
                                                  *(a_geoservers[i]),m_nCellMax, a_maxCoarsenings,fixRegularNextToMultiValued);
          //no way this works for more than two phase
          elevels[i]->setBoundaryPhase( i % phases);
        }
        else
          {
            elevels[i] = m_ebis[i]->buildNextLevel(*a_geoservers[i],fixRegularNextToMultiValued);
          }

        if (elevels[i] == NULL)
        {
          deeper = false;
          break;
        }
        elevels[i]->m_phase = i;

        //remember we are only dealing with TWO phases so
        //i can just decide to stitch when i==1 (both phases have been constructed)
        if (i == 1)
          {
            elevels[1]->reconcileIrreg(*(elevels[0]));
            elevels[1]->levelStitch(*(elevels[0]), elevelsFine[1], elevelsFine[0]);
          }
      }
    if (deeper)
    {
      elevelsFine = elevels;
      nlevels++;
    }
  }

  for (int i=0; i<phases; i++)
    {
      m_ebis[i]->resetLevels(nlevels);
    }
}

void MFIndexSpace::fillEBISLayout(Vector<EBISLayout>& a_ebis,
                                  const DisjointBoxLayout& a_grids,
                                  const Box& a_domain,
                                  const int & nghost) const
{
  CH_assert(a_ebis.size() == m_ebis.size());
  for (int i=0; i<m_ebis.size(); i++)
  {
    m_ebis[i]->fillEBISLayout(a_ebis[i], a_grids,a_domain,nghost);
  }
}

void MFIndexSpace::fillEBISLayout(EBISLayout& a_ebis, int phase,
                                  const DisjointBoxLayout& a_grids,
                                  const Box& a_domain,
                                  const int & nghost) const
{

    m_ebis[phase]->fillEBISLayout(a_ebis, a_grids,a_domain,nghost);

}

IntVectSet
MFIndexSpace::interfaceRegion(int depth) const
{
  IntVectSet rtn;
  for (int i=0; i<m_ebis.size(); i++)
    {
      rtn |= m_ebis[i]->irregCells(depth);
    }
  return rtn;
}
#include "NamespaceFooter.H"
