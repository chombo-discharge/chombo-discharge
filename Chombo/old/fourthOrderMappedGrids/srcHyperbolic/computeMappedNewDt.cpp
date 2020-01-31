#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "computeMappedNewDt.H"
#include "DataIterator.H"
#include "EdgeToCell.H"

#include "NamespaceHeader.H"

// helper function to compute max save dt for advection
Real computeMappedNewDt(const LevelData<FluxBox>& a_faceVel,
                        const MultiBlockCoordSys* a_coordSysPtr,
                        Real a_cfl)
{
  const DisjointBoxLayout grids = a_faceVel.getBoxes();
  // CH_assert(grids == a_coordSysPtr->grids());

  Real maxVelLoc = 0.;
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bx = grids[dit];
      const FluxBox& faceVelFlub = a_faceVel[dit];

      /*
        Get face-centered N, on bx grown by 1.
      */
      Box grownBox1 = grow(bx, 1);
      int nFaceMetricTerms = SpaceDim*SpaceDim;
      FluxBox faceMetricTerms(grownBox1, nFaceMetricTerms);

      const NewCoordSys* coordSysBlockPtr =
        a_coordSysPtr->getCoordSys(bx);

      coordSysBlockPtr->getN(faceMetricTerms, grownBox1);

      // compute area-weighted normal velocities -- 2nd order should
      // be good enough
      bool doFourthOrder = false;
      // LevelData<FluxBox> faceNormalVel(grids, 1);
      FluxBox faceNormalVelFlub(bx, 1);
      coordSysBlockPtr->computeMetricTermProductAverage(faceNormalVelFlub,
                                                        faceVelFlub,
                                                        faceMetricTerms,
                                                        bx,
                                                        doFourthOrder);

      // now average back to cell-centers and divide by cell volumes.
      // instead of averaging face->cell we pick
      // the max absolute value on the two opposing faces.
      // const LevelData<FArrayBox>& cellVolumes = a_coordSysPtr->getCellVolumes();
      // Need cell volume on base box only.
      FArrayBox cellVolumesFab(bx, 1);
      coordSysBlockPtr->cellVol(cellVolumesFab, faceMetricTerms, bx);

      FArrayBox cellVelFab(bx, 1);
      for (int dir=0; dir<SpaceDim; dir++)
        {
          int faceComp = 0;
          int cellComp = 0;
          // take max face velocity on the two faces
          // use absolute value
          faceNormalVelFlub[dir].abs();
          EdgeToCellMax(faceNormalVelFlub, faceComp,
                        cellVelFab, cellComp,
                        dir);
          cellVelFab.divide(cellVolumesFab, 0, 0, 1);

          // now compute maxVelLoc on this box
          // note that this is essentially an estimate of max(vel/dx)
          Real thisMax = cellVelFab.norm(0,0,1);
          if (thisMax > maxVelLoc) maxVelLoc = thisMax;
        }

    }

  Real maxVel = maxVelLoc;
#ifdef CH_MPI
  MPI_Allreduce(&maxVelLoc, &maxVel, 1, MPI_CH_REAL,
                MPI_MAX, Chombo_MPI::comm);
#endif

  Real smallVal = 1.0e-15;
  Real newDt = 0;
  if (maxVel > smallVal)
    {
      newDt = a_cfl/maxVel;
    }

  return newDt;

}

#include "NamespaceFooter.H"
