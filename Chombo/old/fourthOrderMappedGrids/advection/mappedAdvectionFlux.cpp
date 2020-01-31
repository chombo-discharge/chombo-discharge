#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif



#include "mappedAdvectionFlux.H"
#include "mappedLimiter.H"
#include "FourthOrderUtil.H"
#include "BoxIterator.H"

#include "mappedAdvectionF_F.H"
#include "FourthOrderUtil.H"

#include "NamespaceHeader.H"


void
computeMappedFourthOrderFlux(LevelData<FluxBox>& a_compCoordFlux,
                             LevelData<FluxBox>& a_faceVel,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>* a_cellAdvVel,
                             const FourthOrderCoordSys* a_FOCS,
                             const RealVect& a_dx,
                             bool a_limitFaceValues
                             )
{
  const DisjointBoxLayout grids = a_phi.getBoxes();

  // Compute computational-grid cell averages
  LevelData<FArrayBox> phi_cg(grids, a_phi.nComp(), a_phi.ghostVect() );
  cellUJToCellU( phi_cg, a_phi, a_FOCS );
  phi_cg.exchange();
//  LevelData<FArrayBox> vel_cg( grids, a_cellAdvVel->nComp(), a_cellAdvVel->ghostVect() );
//  cellUJToCellU( vel_cg, *a_cellAdvVel );

  // construct appropriately accurate face-averages of phi and advVel
  LevelData<FluxBox> facePhi( grids, phi_cg.nComp(), phi_cg.ghostVect() );
  fourthOrderCellToFace(facePhi, phi_cg );

  //computeFaceAverages( facePhi, phi_cg );
  facePhi.exchange();

//  computeFaceAverages( faceVel, vel_cg );
  fourthOrderCellToFace( a_faceVel, *a_cellAdvVel );
  //computeFaceAverages( a_faceVel, *a_cellAdvVel );
  a_faceVel.exchange();

  // if we're limiting the face-centered values, do it here
  if (a_limitFaceValues)
  {
     applyMappedLimiter(facePhi, phi_cg, a_faceVel, a_FOCS, a_dx);
  }

  // compute computational-space fluxes -
  //   need all three fluxes on each face, so compCoordFlux has
  //   dimension (SpaceDim * nComp)
  CH_assert(a_compCoordFlux.nComp() == SpaceDim*phi_cg.nComp());

  computeCompFaceFluxes( a_compCoordFlux, facePhi, a_faceVel, true );
  a_compCoordFlux.exchange();
}


void
cellUJToCellU(
              LevelData<FArrayBox>& a_cell_u,
              const LevelData<FArrayBox>& a_cell_uJ,
              const FourthOrderCoordSys* a_FOCS)
  {
   int ncomp = a_cell_u.nComp();
   CH_assert(a_cell_uJ.nComp() == ncomp);

   // Obtain J
   const LevelData<FArrayBox>& J = a_FOCS->getJ();

   DataIterator dit = a_cell_uJ.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FArrayBox& thisU = a_cell_u[dit];
      const FArrayBox& thisUJ = a_cell_uJ[dit];
      const FArrayBox& thisJ = J[dit];
      FArrayBox thisJinv(thisJ.box(), 1);
      thisJinv.copy(thisJ);
      thisJinv.invert(1.0);

      Box intersectBox(thisUJ.box());
      intersectBox &= thisJinv.box();
      intersectBox &= thisU.box();

      thisU.copy(thisUJ, intersectBox);
      BoxIterator ibit(intersectBox);
      for (ibit.begin(); ibit.ok(); ++ibit)
      {
         const IntVect& iv = ibit();
         const Real val = thisJinv(iv);
         for (int n=0; n<ncomp; n++)
            thisU(iv,n) *= val;
      }

      // compute box over which we can do this
      Box gradIntersectBox(thisU.box());
      gradIntersectBox.grow(-1);
      gradIntersectBox &= thisJ.box();
      FArrayBox gradProduct(gradIntersectBox, ncomp);
      gradProduct.setVal( 0.0 );

      // since dx's cancel, use 1 here
      Real fakeDx = 1.0;
      Real factor = -1.0/12.0;
      for (int dir=0; dir<SpaceDim; dir++)
        {
          // call fortran wrapper function to handle multidim case
          incrementGradProduct(gradProduct,
                               thisU,
                               thisJ,
                               gradIntersectBox,
                               fakeDx,
                               factor,
                               dir);

        }

      BoxIterator jbit(gradIntersectBox);
      for (jbit.begin(); jbit.ok(); ++jbit)
      {
         const IntVect& iv = jbit();
         const Real val = thisJinv(iv);
         for (int n=0; n<ncomp; n++)
            gradProduct(iv,n) *= val;
      }

      thisU.plus(gradProduct, gradIntersectBox, 0, 0, ncomp);
   } // end loop over boxes
}




void
computeCompFaceFluxes( LevelData<FluxBox>& a_uTimesV,
                       const LevelData<FluxBox>& a_u,
                       const LevelData<FluxBox>& a_v,
                       bool a_useFourthOrder)
{
   // Compute the SpaceDim-by-nComp face-averaged fluxes in computational
   // space, where a_v is the SpaceDim-dimensional velocity vector and
   // a_u is the nComp-dim state vector
   int ncomp = a_u.nComp();
   CH_assert(a_v.nComp() == SpaceDim);
   CH_assert(a_uTimesV.nComp() == SpaceDim * ncomp);

   // loop over boxes
   DataIterator dit = a_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FluxBox& thisUV = a_uTimesV[dit];
      const FluxBox& thisU = a_u[dit];
      const FluxBox& thisV = a_v[dit];

      // loop over faces (index "d" in notes)
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
      {
         FArrayBox& thisUVdir = thisUV[faceDir];
         const FArrayBox& thisUdir = thisU[faceDir];
         const FArrayBox& thisVdir = thisV[faceDir];

         // compute <u_p><v_d> tensor
         Box intersectBox(thisUdir.box());
         intersectBox &= thisVdir.box();
         intersectBox &= thisUVdir.box();

         thisUVdir.setVal(0.0);
         FORT_INCREMENTFACEPROD(CHF_FRA(thisUVdir),
                                CHF_CONST_FRA(thisUdir),
                                CHF_CONST_FRA(thisVdir),
                                CHF_BOX(intersectBox));

         if (a_useFourthOrder)
         {
            // now increment with product of tangential gradients
            // (sum over index "d'" in notes)
            for (int tanDir=0; tanDir<SpaceDim; tanDir++)
            {
               if (tanDir!=faceDir)
               {
                  // since dx's cancel, use 1 here
                  Real fakeDx = 1.0;
                  Real factor = 1.0/12.0;

                  // compute box over which we can do this
                  Box gradIntersectBox(thisUdir.box());
                  gradIntersectBox &= thisVdir.box();
                  gradIntersectBox.grow(tanDir,-1);
                  gradIntersectBox &= thisUVdir.box();

                  FORT_INCREMENTFACEPRODGRAD(CHF_FRA(thisUVdir),
                                             CHF_CONST_FRA(thisUdir),
                                             CHF_CONST_FRA(thisVdir),
                                             CHF_BOX(gradIntersectBox),
                                             CHF_REAL(fakeDx),
                                             CHF_REAL(factor),
                                             CHF_INT(tanDir));

               } // end if tangential direction
            } // end loop over tanDir
         } // end if fourth order
      } // end loop over faceDir
   } // end loop over boxes
}


#include "NamespaceFooter.H"
