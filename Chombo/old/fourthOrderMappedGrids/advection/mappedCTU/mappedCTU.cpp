#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "mappedCTU.H"
#include "MappedAdvectionPhysics.H"
#include "EdgeToCell.H"
#include "PositivityFixF_F.H"

#include "NamespaceHeader.H"

// function to manage computation of CTU fluxes
void
computeCTUFlux(LevelData<FluxBox>& a_CTUflux,
               const LevelData<FArrayBox>& a_oldPhi,
               const LevelData<FluxBox>& a_faceVel,
               FourthOrderCoordSys* a_coordSysPtr,
               PatchGodunov& a_CTUpatchGod,
               Real a_time,
               Real a_dt)
{

  const DisjointBoxLayout& levelGrids = a_CTUflux.getBoxes();

  MappedAdvectionPhysics* advectPhys =
    dynamic_cast<MappedAdvectionPhysics*>(a_CTUpatchGod.getGodunovPhysicsPtr());


  if (advectPhys == NULL)
    {
      MayDay::Error("AMRLevelAdvect::computeLowOrderFlux -- unable to upcast GodunovPhysics to MappedAdvectionPhysics for CTU operator");
    }


  a_CTUpatchGod.setCurrentTime(a_time);

  // get geometry info
  LevelData<FluxBox> NTvel(levelGrids, 1, a_faceVel.ghostVect());
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


  DataIterator dit = a_oldPhi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& thisBox = levelGrids[dit];
      const FArrayBox& thisOldPhi = a_oldPhi[dit];
      FluxBox& thisNTvel = NTvel[dit];

      FArrayBox averagedCellVel(thisNTvel.box(), SpaceDim);
      EdgeToCell(thisNTvel, averagedCellVel);

      // local storage...
      Box srcBox(thisBox);
      srcBox.grow(1);
      FArrayBox srcFab(srcBox, a_oldPhi.nComp());
      srcFab.setVal(0.0);
      FluxBox& thisFlux = a_CTUflux[dit];
      FluxBox thisPhiHalf(thisFlux.box(), thisFlux.nComp());

      a_CTUpatchGod.setCurrentBox(thisBox);

      advectPhys->setCellNTvelPtr(&averagedCellVel);
      advectPhys->setNTvelPtr(&thisNTvel);

      advectPhys->setMappingParameters(&(J[dit]),
                                       &(N[dit]),
                                       &(NJinverse[dit]));

      a_CTUpatchGod.computeWHalf(thisPhiHalf, thisOldPhi, srcFab,
                                 a_dt, thisBox);

      // check for positivity and fix if necessary
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FORT_POSITIVITYFIX(CHF_FRA(thisPhiHalf[dir]),
                             CHF_BOX(thisPhiHalf[dir].box()));
        }

      // now compute flux. We should be able to use the
      // MappedAdvectionPhysics getFlux function for this...
      for (int dir=0; dir<SpaceDim; dir++)
        {
          advectPhys->getFlux(thisFlux[dir],
                              thisPhiHalf[dir],
                              dir,
                              thisFlux[dir].box());
        }
    } // end loop over boxes to compute CTU flux.

}

#include "NamespaceFooter.H"
