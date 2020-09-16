#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
using std::pow;
using std::sqrt;
#ifdef CH_MPI
#include "mpi.h"
#endif

#include "computeNormWeighted.H"
#include "computeNorm.H"
#include "LayoutIterator.H"
#include "BoxLayoutData.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

Real normWeighted(const BoxLayoutData<FArrayBox>& a_layout,
                  const BoxLayoutData<FArrayBox>& a_J,
                  const Interval& a_intvl,
                  const int& a_p)
{
  const BoxLayout validRegion = a_layout.boxLayout();
  DataIterator it = validRegion.dataIterator();
  int ncomp = a_intvl.size();
  int layoutComp = a_intvl.begin();
  CH_assert(ncomp == 1); // simpler this way
  Real normOut;
  if (a_p == 0)  // max norm
    {
      normOut = norm(a_layout, a_intvl, 0);
    }
  else
    {
      Real pwr = Real(a_p);
      normOut = 0.;
      for (;it.ok(); ++it)
        {
          // Without weighting,
          // normOut += a_layout[it()].norm(validRegion[it()], 1,
          // interval.begin(), interval.size());
          const Box& bx = validRegion[it];
          const FArrayBox& layoutFab = a_layout[it];
          const FArrayBox& JFab = a_J[it];
          if (a_p == 1) // abs sum norm
            {
              FArrayBox weightedFab(bx, 1);
              weightedFab.copy(layoutFab, bx, layoutComp, bx, 0, 1);
              weightedFab.mult(JFab, bx, 0, 0, 1);
              normOut += weightedFab.norm(bx, 1, 0, 1);
            }
          else if (a_p == 2)
            {
              for (BoxIterator bit(bx); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  Real layoutVal = layoutFab(iv, layoutComp);
                  Real weightVal = JFab(iv, 0);
                  normOut += layoutVal * layoutVal * weightVal;
                }
            }
          else // not 0 or 1 or 2
            {
              for (BoxIterator bit(bx); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  Real layoutVal = layoutFab(iv, layoutComp);
                  Real weightVal = JFab(iv, 0);
                  normOut += pow(layoutVal, pwr) * weightVal;
                }
            }
        }
#ifdef CH_MPI
      Real recv;
      int result = MPI_Allreduce(&normOut, &recv, 1, MPI_CH_REAL,
                                 MPI_SUM, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
        {
          //bark!!!
          MayDay::Error("sorry, but I had a communcation error on norm");
        }
      normOut = recv;
#endif
      if (a_p == 2)
        {
          normOut = sqrt(normOut);
        }
      else if (a_p != 1) // if (a_p == 1) continue
        {
          Real invpwr = 1.0 / pwr;
          normOut = pow(normOut, invpwr);
        }
    }
  return normOut;
}


Real computeNormWeighted(const Vector<LevelData<FArrayBox>* >& a_phi,
                         const Vector<LevelData<FArrayBox>* >& a_J,
                         const Vector<int>&                    a_nRefFine,
                         const Real                            a_dxCrse,
                         const Interval                        a_comps,
                         const int                             a_p,
                         const int                             a_lBase)
{
  Real normOut;
  if (a_p == 0)
    { // max norm does not require a_J
      normOut = computeNorm(a_phi, a_nRefFine, a_dxCrse, a_comps, a_p, a_lBase);
    }
  else
    {
      Real dxLevel = a_dxCrse;
      int numLevels = a_phi.size();
      Real normLevel;

      // it is often the case that while a_phi has many possible
      // levels, only a subset of them are defined -- check that
      // just to be sure
      if (a_phi[numLevels-1] == NULL)
        {
          int lev = numLevels-1;
          while (a_phi[lev] == NULL)
            {
              lev--;
            }
          numLevels = lev+1;
        }

      normOut = 0.0;

      // loop over levels
      for (int lev = a_lBase; lev < numLevels; lev++)
        {
          int refRatio = -1;
          //in case there are extra levels which are not defined
          if (a_phi[lev] != NULL)
            {
              CH_assert(a_phi[lev]->isDefined());
              LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
              LevelData<FArrayBox>& thisJ = *(a_J[lev]);
              const DisjointBoxLayout* finerGridsPtr = NULL;
              if (lev < numLevels-1)
                {
                  finerGridsPtr = &(a_phi[lev+1]->getBoxes());
                  refRatio = a_nRefFine[lev];
                }
              normLevel = computeNormWeighted(thisPhi, thisJ, finerGridsPtr, refRatio,
                                              dxLevel, a_comps, a_p);
              Real pwr = a_p;
              normLevel = pow(normLevel, pwr);
              normOut += normLevel;
            }
          // update  dxLevel
          dxLevel = dxLevel/refRatio;
        }
      // shouldn't need to do broadcast/gather thing

      // now take the inverse power of the norm
      Real pwr = Real(a_p);
      Real invp = 1.0 / pwr;
      normOut = pow(normOut, invp);
    }
  return normOut;
}

Real computeNormWeighted(const LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>& a_J,
                         const DisjointBoxLayout*    a_finerGridsPtr,
                         const int                   a_nRefFine,
                         const Real                  a_dx,
                         const Interval              a_comps,
                         const int                   a_p)
{
  Real normLevel = 0.0;
  if (a_p == 0)
    { // max norm does not require a_J
      normLevel = computeNorm(a_phi, a_finerGridsPtr, a_nRefFine, a_dx, a_comps, a_p);
    }
  else
    {
      const DisjointBoxLayout& levelGrids = a_phi.getBoxes();
      LevelData<FArrayBox> temp(levelGrids, a_comps.size());
      int ncomp = a_comps.size();
      Interval tempComps(0, ncomp-1);
      DataIterator dit = a_phi.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox& tempFab = temp[dit];
          Box bx = tempFab.box();
          tempFab.copy(bx, tempComps, bx, a_phi[dit], a_comps);
          if (a_finerGridsPtr != NULL)
            {
              LayoutIterator litFine = a_finerGridsPtr->layoutIterator();
              // now loop over fine boxes and set covered regions to 0
              for (litFine.reset(); litFine.ok(); ++litFine)
                {
                  Box coveredBox(a_finerGridsPtr->get(litFine()));
                  coveredBox.coarsen(a_nRefFine);
                  coveredBox &= bx;
                  if (!coveredBox.isEmpty())
                    {
                      tempFab.setVal(0.0, coveredBox, 0, ncomp);
                    }
                }  // end loop over fine-grid boxes
            } // end if there is a finer level
        } // end loop over this level's grids
      // compute unscaled norm
      normLevel = normWeighted(temp, a_J, tempComps, a_p);
      // now scale it
      Real exponent = SpaceDim;
      exponent /= a_p;
      Real scale = pow(a_dx, exponent);
      normLevel *= scale;
    }
  return normLevel;
}

// anisotropic version
Real computeNormWeighted(const Vector<LevelData<FArrayBox>* >& a_phi,
                         const Vector<LevelData<FArrayBox>* >& a_J,
                         const Vector<IntVect>&                a_nRefFine,
                         const RealVect&                       a_dxCrse,
                         const Interval                        a_comps,
                         const int                             a_p,
                         const int                             a_lBase)
{
  Real normOut;
  if (a_p == 0)
    { // max norm does not require a_J
      normOut = computeNorm(a_phi, a_nRefFine, a_dxCrse, a_comps, a_p, a_lBase);
    }
  else
    {
      Real pwr = Real(a_p);
      int numLevels = a_phi.size();
      Real normLevel;
      RealVect dxLevel = a_dxCrse;
      // it is often the case that while a_phi has many possible
      // levels, only a subset of them are defined -- check that
      // just to be sure
      if (a_phi[numLevels-1] == NULL)
        {
          int lev = numLevels-1;
          while (a_phi[lev] == NULL)
            {
              lev--;
            }
          numLevels = lev+1;
        }
      // loop over levels
      for (int lev = a_lBase; lev < numLevels; lev++)
        {
          IntVect refRatio = -IntVect::Unit;
          //in case there are extra levels which are not defined
          if (a_phi[lev] != NULL)
            {
              CH_assert(a_phi[lev]->isDefined());
              const LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
              const LevelData<FArrayBox>& thisJ = *(a_J[lev]);
              const DisjointBoxLayout* finerGridsPtr = NULL;
              if (lev < numLevels-1)
                {
                  finerGridsPtr = &(a_phi[lev+1]->getBoxes());
                  refRatio = a_nRefFine[lev];
                }
              normLevel = computeNormWeighted(thisPhi, thisJ, finerGridsPtr, refRatio,
                                              dxLevel, a_comps, a_p);
              normLevel = pow(normLevel, pwr);
              normOut += normLevel;
            }
          // update  dxLevel
          dxLevel = dxLevel/RealVect(refRatio);
        }
      // shouldn't need to do broadcast/gather thing
      // now take the inverse power of the norm
      if (a_p != 1)
        {
          Real invp = 1.0 / pwr;
          normOut = pow(normOut, invp);
        }
    }
  return normOut;
}

// anisotropic version
Real computeNormWeighted(const LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>& a_J,
                         const DisjointBoxLayout*    a_finerGridsPtr,
                         const IntVect&              a_nRefFine,
                         const RealVect&             a_dx,
                         const Interval              a_comps,
                         const int                   a_p)
{
  Real normLevel;
  if (a_p == 0)
    { // max norm does not require a_J
      normLevel = computeNorm(a_phi, a_finerGridsPtr, a_nRefFine, a_dx, a_comps, a_p);
    }
  else
    {
      Real pwr = Real(a_p);
      const DisjointBoxLayout& levelGrids = a_phi.getBoxes();
      LevelData<FArrayBox> temp(levelGrids, a_comps.size());
      int ncomp = a_comps.size();
      Interval tempComps(0,ncomp-1);
      DataIterator dit = a_phi.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox& tempFab = temp[dit];
          Box bx = tempFab.box();
          tempFab.copy(bx, tempComps, bx, a_phi[dit], a_comps);
          if (a_finerGridsPtr != NULL)
            {
              LayoutIterator litFine = a_finerGridsPtr->layoutIterator();
              // now loop over fine boxes and set covered regions to 0
              for (litFine.reset(); litFine.ok(); ++litFine)
                {
                  Box coveredBox(a_finerGridsPtr->get(litFine()));
                  coveredBox.coarsen(a_nRefFine);
                  coveredBox &= bx;
                  if (!coveredBox.isEmpty())
                    {
                      tempFab.setVal(0.0, coveredBox, 0, ncomp);
                    }
                }  // end loop over fine-grid boxes
            } // end if there is a finer level
        } // end loop over this level's grids
      
      // compute unscaled norm
      // FIXME
      normLevel = normWeighted(temp, a_J, tempComps, a_p);
      // now scale it
      Real dV = a_dx.product();
      Real exponent = 1. / pwr;
      Real scale = pow(dV, exponent);
      normLevel *= scale;
    }

  return normLevel;
}

#include "NamespaceFooter.H"
