#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "computeMappedDt.H"
#include "DataIterator.H"
#include "EdgeToCell.H"

#include "NamespaceHeader.H"

// helper function to compute max save dt for advection
Real computeMappedDt(const LevelData<FluxBox>& a_faceVel,
                     const FourthOrderCoordSys* a_coordSysPtr,
                     Real a_cfl,
                     const bool velocities_are_physical)
{
  const DisjointBoxLayout grids = a_faceVel.getBoxes();
  CH_assert(grids == a_coordSysPtr->grids());

  DataIterator dit = grids.dataIterator();

  // compute area-weighted normal velocities -- 2nd order should
  // be good enough
  bool doFourthOrder = false;
  LevelData<FluxBox> faceNormalVel(grids, 1);

  if (!velocities_are_physical && SpaceDim > 1)
     {
        a_coordSysPtr->computeMetricTermProductAverage(faceNormalVel,
                                                       a_faceVel,
                                                       doFourthOrder);
     }
  else
     {
        const RealVect& dx = a_coordSysPtr->dx();

        RealVect face_area;
        for (int dir=0; dir<SpaceDim; ++dir)
        {
           face_area[dir] = 1.;
           for (int tdir=0; tdir<SpaceDim; ++tdir)
           {
              if (tdir != dir) face_area[dir] *= dx[tdir];
           }
        }

        for (dit.begin(); dit.ok(); ++dit)
        {
           const FluxBox& thisFaceVel = a_faceVel[dit];
           FluxBox& thisNormalVel = faceNormalVel[dit];
           for (int dir=0; dir<SpaceDim; ++dir)
           {
              FArrayBox& thisNormalVel_dir = thisNormalVel[dir];
              thisNormalVel_dir.copy(thisFaceVel[dir],dir,0,1);
              thisNormalVel_dir *= face_area[dir];
           }
        }
     }

  // now average back to cell-centers and divide by cell volumes.
  // instead of averaging face->cell we pick
  // the max absolute value on the two opposing faces.
  const LevelData<FArrayBox>& cellVolumes = a_coordSysPtr->getCellVolumes();

  Real maxVelLoc = 0;
  //  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox cellVel(grids[dit], 1);
      for (int dir=0; dir<SpaceDim; dir++)
        {
          int faceComp = 0;
          int cellComp = 0;
          // take max face velocity on the two faces
          // use absolute value
          faceNormalVel[dit][dir].abs();
          EdgeToCellMax(faceNormalVel[dit], faceComp,
                        cellVel, cellComp,
                        dir);
          cellVel.divide(cellVolumes[dit],
                         0,0, 1);

          // now compute maxVelLoc on this box
          // note that this is essentially an estimate of max(vel/dx)
          Real thisMax = cellVel.norm(0,0,1);
          if (thisMax > maxVelLoc) maxVelLoc = thisMax;
        }

    }

  Real maxVel = maxVelLoc;
#ifdef CH_MPI
  MPI_Allreduce(&maxVelLoc, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
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
