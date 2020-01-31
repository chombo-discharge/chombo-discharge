#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GaussianAnalyticMappedTaggingStrategy.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
GaussianAnalyticMappedTaggingStrategy::
GaussianAnalyticMappedTaggingStrategy():
  AMRLevelMappedTaggingStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
GaussianAnalyticMappedTaggingStrategy::
~GaussianAnalyticMappedTaggingStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
IntVectSet
GaussianAnalyticMappedTaggingStrategy::
cellsToRefine(int a_level,
              const ProblemDomain& a_problem_domain,
              const DisjointBoxLayout& a_grids,
              LevelGridMetrics* a_gridMetrics,
              const MOLPhysics& a_physics,
              const MappedLevelData& a_data,
              Real a_dx,
              Real a_time) const
{
  IntVectSet localTags;

  // This is special code for analytically refining the advection of a gaussian
  // profile (see paper ...)
  RealVect p(D_DECL(0.5, 0.5, 0.5));
  RealVect av(D_DECL(1.0, 0.5, 0.0));
  p += a_time*av;
  IntVect offset(D_DECL((int)p[0], (int)p[1], (int)p[2]));
  D_TERM(p[0] -= (Real)offset[0];,
      p[1] -= (Real)offset[1];,
      p[2] -= (Real)offset[2];)
    Real r = 0.;
  switch (a_level)
  {
    case 0:
      r = 0.35;
      //      r = 0.425;
      break;
    case 1:
      r = 0.225;
      break;
  }
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
  {
    const Box& bx = a_grids[dit];
    const NewCoordSys *const coordSys = a_gridMetrics->getCoordSys(bx);
    ShiftIterator shiftIt = a_problem_domain.shiftIterator();

    for (BoxIterator bit(bx); bit.ok(); ++bit)
    {
      RealVect Xi = bit();
      Xi += 0.5;
      Xi *= a_dx;
      const RealVect X = coordSys->realCoord(Xi);
      RealVect len = X - p;
      if (len.vectorLength() <= r)
      {
        localTags |= bit();
      }
      else
      {
        for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
        {
          len = X;
          len += shiftIt();
          len -= p;
          if (len.vectorLength() <= r)
          {
            localTags |= bit();
            break;
          }
        }
      }
    }
  }
  return localTags;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AMRLevelMappedTaggingStrategy*
GaussianAnalyticMappedTaggingStrategy::
clone() const
{
  return new GaussianAnalyticMappedTaggingStrategy();
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

