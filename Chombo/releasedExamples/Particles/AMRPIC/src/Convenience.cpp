#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Convenience.H"
#include "conveniencefort_F.H"

#include "NamespaceHeader.H"

/// Interpolate the data to an intermediate time
void interpolateInTime(LevelData<FArrayBox>&       a_phi,
                       const LevelData<FArrayBox>& a_phiOld,
                       const LevelData<FArrayBox>& a_phiNew,
                       const Real&                 a_time,
                       const Real&                 a_tNew,
                       const Real&                 a_tOld,
                       const Real&                 a_dt,
                       const Interval&             a_srcComp,
                       const Interval&             a_dstComp,
                       const IntVect&              a_ghost)
{
  CH_assert(a_srcComp.size()==a_dstComp.size());
  CH_assert(a_dt >= 0.0);

  const int srcComp = a_srcComp.begin();
  const int dstComp = a_dstComp.begin();
  const int numComp = a_dstComp.size();

  const Real eps = 0.01* a_dt;

  const DisjointBoxLayout& grids = a_phi.getBoxes();
  if (a_tNew-a_time<=eps)   // case alpha=1
    {
      for (DataIterator di=a_phi.dataIterator(); di.ok(); ++di)
        {
          const Box box = grow(grids.get(di),a_ghost);
          a_phi[di].copy(a_phiNew[di],box,srcComp,box,dstComp,numComp);
        }
    }
  else if (a_time-a_tOld<=eps)   // case alpha=0
    {
      for (DataIterator di=a_phi.dataIterator(); di.ok(); ++di)
        {
          const Box box = grow(grids.get(di),a_ghost);
          a_phi[di].copy(a_phiOld[di],box,srcComp,box,dstComp,numComp);
        }
    }
  else
    {
      CH_assert(a_tNew-a_tOld>=eps && a_time>a_tOld);
      const Real alpha = (a_time-a_tOld)/(a_tNew-a_tOld);

      for (DataIterator di=a_phi.dataIterator(); di.ok(); ++di)
        {
          const FArrayBox& phiNew = a_phiNew[di];
          const FArrayBox& phiOld = a_phiOld[di];

          FArrayBox& phi = a_phi[di];
          const Box box = grow(grids.get(di),a_ghost);

          CH_assert(phiNew.box().contains(box));
          CH_assert(phiOld.box().contains(box));

          FORT_INTERPOLATEINTIME( CHF_FRA(phi),
                                  CHF_CONST_FRA(phiNew),
                                  CHF_CONST_FRA(phiOld),
                                  CHF_CONST_REAL(alpha),
                                  CHF_CONST_INT(srcComp),
                                  CHF_CONST_INT(dstComp),
                                  CHF_CONST_INT(numComp),
                                  CHF_BOX(box));
        }
    }
}

#include "NamespaceFooter.H"
