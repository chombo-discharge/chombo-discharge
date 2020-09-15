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

#include "LayoutIterator.H"
#include "BoxLayoutData.H"

#include "computeMappedNorm.H"
#include "NamespaceHeader.H"

Real computeNorm(const Vector<LevelData<FArrayBox>* >& a_phi,
                 const Vector<int>&                    a_nRefFine,
                 const Vector<CoordSys<FArrayBox,FluxBox>* >& a_coordSys,
                 const Interval                        a_comps,
                 const int                             a_p,
                 const int                             a_lBase)
{
  int numLevels = a_phi.size();
  Real norm;
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

  norm = 0.0;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
  {
    int refRatio = -1;
    //in case there are extra levels which are not defined
    if (a_phi[lev] != NULL)
    {
      CH_assert(a_phi[lev]->isDefined());

      LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
      const DisjointBoxLayout* finerGridsPtr = NULL;

      if (lev < numLevels-1)
      {
        finerGridsPtr = &(a_phi[lev+1]->getBoxes());
        refRatio = a_nRefFine[lev];
      }

      normLevel = computeNorm(thisPhi, finerGridsPtr, refRatio,
                              a_coordSys[lev], a_comps, a_p);

      if (a_p != 0)
      {
        Real p = a_p;
        normLevel = pow(normLevel,p);

        norm += normLevel;
      }
      else if (normLevel > norm)
      {
        norm = normLevel;
      }
    }

  }

  // shouldn't need to do broadcast/gather thing

  // now take the inverse power of the norm
  if (a_p != 0)
  {
    Real invp = 1.0/a_p;
    norm = pow(norm,invp);
  }

  return norm;
}

Real computeNorm(const LevelData<FArrayBox>& a_phi,
                 const DisjointBoxLayout*    a_finerGridsPtr,
                 const int                   a_nRefFine,
                 const CoordSys<FArrayBox,FluxBox>* a_coordSys,
                 const Interval              a_comps,
                 const int                   a_p)
{
  Real normLevel;
  normLevel = 0.0;

  const DisjointBoxLayout& levelGrids = a_phi.getBoxes();

  LevelData<FArrayBox> temp(levelGrids, a_comps.size());
  const LevelData<FArrayBox>& vol = a_coordSys->getCellVolumes();

  Interval tempComps(0,a_comps.size()-1);

  DataIterator dit = a_phi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    temp[dit()].copy(temp[dit()].box(), tempComps,
                     temp[dit()].box(), a_phi[dit()], a_comps);

    // multiply by cell volumes if not maxnorm
    if (a_p != 0)
      {
        for (int comp=0; comp<=tempComps.end(); comp++)
          {
            temp[dit].mult(vol[dit],temp[dit].box(),
                           0,comp,1);
          }
      }

    if (a_finerGridsPtr != NULL)
    {
      LayoutIterator litFine = a_finerGridsPtr->layoutIterator();

      // now loop over fine boxes and set covered regions to 0
      for (litFine.reset(); litFine.ok(); ++litFine)
      {
        Box coveredBox(a_finerGridsPtr->get(litFine()));
        coveredBox.coarsen(a_nRefFine);
        coveredBox &= temp[dit()].box();

        if (!coveredBox.isEmpty())
        {
          temp[dit()].setVal(0.0, coveredBox, 0, tempComps.size());
        }
      }  // end loop over fine-grid boxes
    } // end if there is a finer level
  } // end loop over this level's grids

  // compute unscaled norm
  normLevel = norm(temp, tempComps, a_p);

  return normLevel;
}


#include "NamespaceFooter.H"
