#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "mappedDCU.H"
#include "mappedDCUF_F.H"
#include "mappedAdvectionFlux.H"

#include "NamespaceHeader.H"

// function to manage computation of DCU fluxes
void
computeDCUFlux(LevelData<FluxBox>& a_DCUflux,
               const LevelData<FArrayBox>& a_oldPhi,
               const LevelData<FluxBox>& a_faceVel,
               FourthOrderCoordSys* a_coordSysPtr,
               Real a_time,
               Real a_dt)
{

  const DisjointBoxLayout& grids = a_DCUflux.getBoxes();

  // get geometry info
  LevelData<FluxBox> NTvel(grids, 1, a_faceVel.ghostVect());
  const LevelData<FluxBox>& N = a_coordSysPtr->getFaceMetricTerms();
  const LevelData<FArrayBox>& J = a_coordSysPtr->getJ();
  const LevelData<FluxBox>& NJinverse = a_coordSysPtr->getNJinverse();

  if (a_faceVel.nComp() == SpaceDim)
  {
     // can probably get away with doing this to 2nd order
     // as well
     a_coordSysPtr->computeMetricTermProductAverage(NTvel,
                                                    a_faceVel);
  }
  else
  {
     // if a_faceVel only has a single component, then assume
     // that it's already N^T*vel
     DataIterator dit = NTvel.dataIterator();
     for (dit.begin(); dit.ok(); ++dit)
     {
        // this is a bit inefficient, but it is easier than
        // doing some sort of aliasing...
        NTvel[dit].copy(a_faceVel[dit]);
     }
  } // end if a_faceVel has a single component

  LevelData<FArrayBox> phi_cg(grids, a_oldPhi.nComp(), a_oldPhi.ghostVect() );
  cellUJToCellU( phi_cg, a_oldPhi, a_coordSysPtr );
  phi_cg.exchange();

  DataIterator dit = phi_cg.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
     const Box& thisBox = grids[dit];
     const FArrayBox& thisPhi = phi_cg[dit];
     FluxBox& thisNTvel = NTvel[dit];
     FluxBox& thisFlux = a_DCUflux[dit];

     for (int dir(0); dir<SpaceDim; dir++)
     {
        FORT_DCUFLUX(CHF_FRA(thisFlux[dir]),
                     CHF_BOX(thisFlux[dir].box()),
                     CHF_CONST_FRA(thisPhi),
                     CHF_CONST_FRA(thisNTvel[dir]),
                     CHF_INT(dir),
                     CHF_REAL(a_dt));
     }
  } // end loop over boxes to compute DCU flux.

}

#include "NamespaceFooter.H"
