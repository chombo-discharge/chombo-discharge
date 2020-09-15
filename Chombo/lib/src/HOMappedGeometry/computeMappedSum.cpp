#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LayoutIterator.H"

#include "computeMappedSum.H"
#include "NamespaceHeader.H"

Real computeSum(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<int>&                    a_nRefFine,
                const Vector<CoordSys<FArrayBox,FluxBox>* >& a_coordSys,
                const Interval&                       a_comps,
                const int&                            a_lBase)
{
  Real volume;

  return computeSum(volume,a_phi,a_nRefFine,a_coordSys,a_comps,a_lBase);
}

Real computeSum(Real&                                 a_volume,
                const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<int>&                    a_nRefFine,
                const Vector<CoordSys<FArrayBox,FluxBox>* >&   a_coordSys,
                const Interval&                       a_comps,
                const int&                            a_lBase)
{
  int numLevels = a_phi.size();
  Real sum, sumLevel;
  Real vol, volLevel;

  sum = 0.0;
  vol = 0.0;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
  {
    int refRatio = -1;
    //in case there are extra levels which are not defined
    if (a_phi[lev] != NULL)
    {
      CH_assert(a_phi[lev]->isDefined());

      LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
      const DisjointBoxLayout* finerGridsPtr;

      if (lev < numLevels-1)
      {
        finerGridsPtr = &(a_phi[lev+1]->getBoxes());
        refRatio = a_nRefFine[lev];
      }
      else
      {
        finerGridsPtr = NULL;
      }

      sumLevel = computeSum(volLevel,thisPhi,finerGridsPtr,refRatio,
                            a_coordSys[lev],a_comps,false);
      sum += sumLevel;
      vol += volLevel;


    }
  }

  // do broadcast/gather thing here
#ifdef CH_MPI
  Real recv;
  int result = MPI_Allreduce(&sum, &recv, 1, MPI_CH_REAL,
                             MPI_SUM, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Sorry, but I had a communication error in computeSum");
  }

  sum = recv;

  result = MPI_Allreduce(&vol, &recv, 1, MPI_CH_REAL,
                         MPI_SUM, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Sorry, but I had a communication error in computeSum");
  }

  vol = recv;
#endif

  a_volume = vol;

  return sum;
}

Real computeSum(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGridsPtr,
                const int&                  a_nRefFine,
                const CoordSys<FArrayBox,FluxBox>*    a_coordSysPtr,
                const Interval&             a_comps,
                const bool                  a_global)
{
  Real volume;

  return computeSum(volume,a_phi,a_finerGridsPtr,a_nRefFine,a_coordSysPtr,
                    a_comps, a_global);

}

Real computeSum(Real&                       a_volume,
                const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGridsPtr,
                const int&                  a_nRefFine,
                const CoordSys<FArrayBox,FluxBox>*    a_coordSys,
                const Interval&             a_comps,
                const bool                  a_global)
{
  Real sum, sumLevel;
  Real vol, volLevel;

  sum = 0.0;
  vol = 0.0;

  const DisjointBoxLayout& levelGrids = a_phi.getBoxes();

  LevelData<FArrayBox> temp(levelGrids,a_comps.size());
  const LevelData<FArrayBox>& cellVol = a_coordSys->getCellVolumes();
  LevelData<FArrayBox> volTemp(levelGrids, 1);

  Interval tempComps(0,a_comps.size()-1);
  Interval volComps(0,0);

  DataIterator dit = a_phi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    Box copyBox = temp[dit()].box();
    temp[dit()].copy(copyBox,tempComps,copyBox,a_phi[dit()],a_comps);
    volTemp[dit].copy(copyBox,volComps, copyBox, cellVol[dit], volComps);

    // now multiply each component by cell volumes
    for (int comp=tempComps.begin(); comp <= tempComps.end(); comp++)
      {
        temp[dit].mult(cellVol[dit],copyBox,0, comp, 1);
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
          temp[dit()].setVal(0.0,coveredBox,0,tempComps.size());
          volTemp[dit()].setVal(0.0,coveredBox,0,1);
        }
      }  // end loop over fine-grid boxes
    } // end if there is a finer level

    sumLevel = temp[dit()].sum(0,tempComps.size());
    sum += sumLevel;

    volLevel = volTemp[dit()].sum(0,1);
    vol += volLevel;
  } // end loop over this level's grids


  // do broadcast/gather thing here
#ifdef CH_MPI
  if (a_global)
  {
    Real recv;
    int result = MPI_Allreduce(&sum, &recv, 1, MPI_CH_REAL,
                               MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in computeSum");
    }

    sum = recv;

    result = MPI_Allreduce(&vol, &recv, 1, MPI_CH_REAL,
                           MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in computeSum");
    }

    vol = recv;
  }
#endif

  a_volume = vol;

  return sum;
}

#include "NamespaceFooter.H"
