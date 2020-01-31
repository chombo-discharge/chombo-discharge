#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include "Ancillae.H"
#include "IntVectSet.H"
#include "SelfGravF_F.H"

// pass a copy cause MPI_Allreduce does not guarantee it's constantness
Real globalMax(Real a_localVal)
{
#ifdef CH_MPI
  Real maxVal;
  int result = MPI_Allreduce(&a_localVal, &maxVal, 1, MPI_CH_REAL,
                             MPI_MAX, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
  {
    MayDay::Error("MPI communcation error in maxVal");
  }
  return maxVal;
#else
  return a_localVal;
#endif
}

// pass a copy cause MPI_Allreduce does not guarantee it's constantness
Real globalMin(Real a_localVal)
{
#ifdef CH_MPI
  Real minVal;
  int result = MPI_Allreduce(&a_localVal, &minVal, 1, MPI_CH_REAL,
                             MPI_MIN, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
  {
    MayDay::Error("MPI communcation error in minVal");
  }
  return minVal;
#else
  return a_localVal;
#endif
}


//
Real levelMin(const LevelData<FArrayBox>& a_levFab,
              const int                   a_comp)
{
  CH_assert(a_comp>=0);

  Real minVal = 1.e12;
  const DisjointBoxLayout& grids = a_levFab.getBoxes();

  for (DataIterator dit=a_levFab.dataIterator(); dit.ok(); ++dit)
  {
    Real tmpMin = a_levFab[dit()].min(grids.get(dit()),a_comp);
    if (tmpMin<minVal) minVal=tmpMin;
  }

  return globalMin(minVal);
}

//
Real levelMax(const LevelData<FArrayBox>& a_levFab,
              const int                   a_comp)
{
  CH_assert(a_comp>=0);

  Real maxVal = zero;
  const DisjointBoxLayout& grids = a_levFab.getBoxes();

  for (DataIterator dit=a_levFab.dataIterator(); dit.ok(); ++dit)
  {
    Real tmpMax = a_levFab[dit()].max(grids.get(dit()),a_comp);
    if (tmpMax>maxVal) maxVal=tmpMax;
  }

  return globalMax(maxVal);
}

Real globalAverage(const LevelData<FArrayBox>& a_levFab,
                   const int                   a_indx,
                   const bool                  a_perVol)
{
  const DisjointBoxLayout& grids = a_levFab.getBoxes();
  Real mean = zero;
  for (DataIterator di = grids.dataIterator(); di.ok(); ++di)
  {
    const FArrayBox& fab = a_levFab[di()];
    const Box&       box = grids.get(di());
    mean     += fab.sum(box,a_indx);
  }

  Real gmean;
#ifdef CH_MPI
  int result0 = MPI_Allreduce(&mean, &gmean, 1, MPI_CH_REAL,
                             MPI_SUM, Chombo_MPI::comm);
  if (result0 != MPI_SUCCESS)
  {
    MayDay::Error("MPI communcation error in globalAverage");
  }
#else
  gmean = mean;
#endif

  if (a_perVol)
  {
    long totCells= 0;
    for (LayoutIterator lit = grids.layoutIterator(); lit.ok(); ++lit)
    {
      totCells += grids[lit()].numPts();
    }
    if (totCells>0)
      gmean /= (Real)totCells;
    else
      gmean =zero;
  }

  return gmean;
}

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
  CH_assert((a_tNew-a_tOld)>=zero);

  const int srcComp = a_srcComp.begin();
  const int dstComp = a_dstComp.begin();
  const int numComp = a_dstComp.size();

  DataIterator dit = a_phi.dataIterator();
  const DisjointBoxLayout& grids = a_phi.getBoxes();

  const Real alpha = (a_time-a_tOld)/(a_tNew-a_tOld);

  const Real eps = 0.01* a_dt;
  if (Abs(one-alpha) < eps) // case alpha=1
  {
    for (dit.begin(); dit.ok(); ++dit)
    {
      Box box = grids.get(dit());
      box.grow(a_ghost);
      a_phi[dit()].copy(a_phiNew[dit()],box,srcComp,box,dstComp,numComp);
    }
  }
  else if (Abs(alpha) < eps) // case alpha=0
  {
    for (dit.begin(); dit.ok(); ++dit)
    {
      Box box = grids.get(dit());
      box.grow(a_ghost);
      a_phi[dit()].copy(a_phiOld[dit()],box,srcComp,box,dstComp,numComp);
    }
  }
  else
  {
    for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& phiNew = a_phiNew[dit()];
      const FArrayBox& phiOld = a_phiOld[dit()];

      FArrayBox& phi = a_phi[dit()];
      Box box = grids.get(dit());
      box.grow(a_ghost);

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


void computePlotVars(LevelData<FArrayBox>&       a_UPlot,
                     const LevelData<FArrayBox>& a_UNew)
{
  // Make sure everything is defined
  CH_assert(a_UPlot.isDefined());
  CH_assert(a_UNew.isDefined());

  // Beginning of loop through patches/grids.
  for (DataIterator dit = a_UNew.dataIterator(); dit.ok(); ++dit)
  {
    FArrayBox& UPlot = a_UPlot[dit()];
    const FArrayBox& UNew = a_UNew[dit()];

    Box box = UPlot.box();

    FORT_COMPUTEPLOTVARS( CHF_FRA(UPlot),
                          CHF_CONST_FRA(UNew),
                          CHF_BOX(box));
  }
}


void computeDeltaU(const LevelData<FArrayBox>& a_UNew,
                   const LevelData<FArrayBox>& a_UOld)
{
  // Make sure everything is defined
  CH_assert(a_UNew.isDefined());
  CH_assert(a_UOld.isDefined());

  int nComp = a_UNew.interval().size();
  Vector<Real> deltaU(nComp,zero);

  // Beginning of loop through patches/grids.
  for (DataIterator dit = a_UNew.dataIterator(); dit.ok(); ++dit)
  {
    // The current box
    Box box = a_UNew[dit()].box();

    const FArrayBox& UNew = a_UNew[dit()];
    const FArrayBox& UOld = a_UOld[dit()];

    FArrayBox dU(UNew.box(),UNew.nComp());
    dU.copy(UNew);
    dU -= UOld;

    for (int i=0; i<nComp; ++i)
    {
      deltaU[i] += dU.norm(box,2,i,1);
    }
  }

  pout() << "AMRLevelSelfGravity:: deltaU "  << '\n';
  for (int i=0; i<nComp; ++i)
  {
    pout() << "AMRLevelSelfGravity:: component "<< i << ": "<< deltaU[i] << '\n';
  }
}

void bufferBoxes(DisjointBoxLayout&       a_bufferGrids,
                 const DisjointBoxLayout& a_grids,
                 const DisjointBoxLayout& a_crseGrids,
                 const int                a_bufferSize,
                 const int                a_refRatio)
{
  Vector<Box> bufBoxVec;
  Vector<int> procIDVec;

  DisjointBoxLayout crseFineGrids;
  coarsen(crseFineGrids,a_grids,a_refRatio);

  pout() << " buffer Grids " << '\n';

  LayoutIterator fli(a_grids.layoutIterator());
  for (LayoutIterator cli = a_crseGrids.layoutIterator(); cli.ok(); ++cli)
  {
    // coarse box
    const Box& crseBox = a_crseGrids.get(cli());
    IntVectSet crseIVS(crseBox);

    // subtract iv's of all fine boxes
    for (fli.begin(); fli.ok(); ++fli)
    {
      // coarsened box
      const Box& crseFineBox = crseFineGrids.get(fli());
      // const int procID = a_crseGrids.procID(cli());

      // optimize search for intersections
      if (crseBox.smallEnd(0) > crseFineBox.bigEnd(0))
      {
        continue;
      }
      else if (crseBox.bigEnd(0) < crseFineBox.smallEnd(0))
      {
        break;
      }

      crseIVS -= crseFineBox;
    }

    pout() << " crseBox " << crseBox << '\n';

    crseIVS.convert();
    Vector<Box> crseBoxSet = crseIVS.boxes();

    for (int ib=0; ib<crseBoxSet.size(); ib++)
    {

      pout() << " ib " << ib << " crseBoxSet " << crseBoxSet[ib] << '\n';

      // grow a bufferSize to intersect the buffer region
      crseBoxSet[ib].grow(a_bufferSize);


      for (fli.begin(); fli.ok(); ++fli)
      {
        // current box
        const Box& crseFineBox = crseFineGrids.get(fli());

        //        pout() << " crseFineBox " << crseFineBox << '\n';


        // optimize search for intersections
        if (crseFineBox.bigEnd(0) < crseBoxSet[ib].smallEnd(0))
        {
          continue;
        }
        else if (crseFineBox.smallEnd(0) > crseBoxSet[ib].bigEnd(0))
        {
          break;
        }

        Box bufBox = crseBoxSet[ib] & crseFineBox;

        if (!bufBox.isEmpty())
        {

          pout() << " bufBox " << bufBox << '\n';

          const int procID = a_grids.procID(fli());
          bufBox.refine(a_refRatio);
          bufBoxVec.push_back(bufBox);
          procIDVec.push_back(procID);
        }
      }
    }
  }
  a_bufferGrids.define(bufBoxVec,procIDVec);
}

