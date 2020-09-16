#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FourthOrderUtil.H"
#include "LoHiSide.H"
#include "FourthOrderUtilF_F.H"

#include "NamespaceHeader.H"


/// -----------------------------------------------------------
void fourthOrderAverage(LevelData<FluxBox>& a_phi,
                        int a_sgn)
{
  CH_TIME("fourthOrderAverage flux");

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisPhi = a_phi[dit];
      fourthOrderAverageFace(thisPhi, a_sgn);
    } // end loop over boxes
}


/// -----------------------------------------------------------
void fourthOrderAverageFace(FluxBox& a_phi,
                            int a_sgn)
{
  CH_TIME("fourtOrderAverageFace");

  // default a_sgn = +1 for convolving face-centered to face-averaged;
  // a_sign = -1 for deconvolving face-averaged to face-centered.
  Real factor = a_sgn * (1.0/24.0);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& thisPhiDir = a_phi[dir];
      // derivBox shrunk by 1 in tangential directions
      Box derivBox(thisPhiDir.box());
      derivBox.grow(-1);
      derivBox.grow(dir, 1);
      FArrayBox tempLap(derivBox, thisPhiDir.nComp());
      tempLap.setVal(0.0);

      for (int tanDir=0; tanDir <SpaceDim; tanDir++)
        {
          if (tanDir != dir)
            {
              FORT_INCREMENTLAPLACIAN(CHF_FRA(tempLap),
                                      CHF_CONST_FRA(thisPhiDir),
                                      CHF_BOX(derivBox),
                                      CHF_CONST_INT(tanDir),
                                      CHF_CONST_REAL(factor));

            } // end if tangential direction
        } // end loop over tangential directions
      thisPhiDir.plus(tempLap, derivBox, 0, 0, tempLap.nComp());
    } // end loop over FluxBox face directions
}


/// -----------------------------------------------------------
void fourthOrderAverageCenterFace(FArrayBox& a_phiAvg,
                                  const FArrayBox& a_phiCen,
                                  int a_dir,
                                  int a_sgn)
{
  CH_TIME("fourtOrderAverageCenterFace");

  int ncomp = a_phiCen.nComp();
  // derivBox shrunk by 1 in tangential directions
  IntVect shrinkVect = BASISV(a_dir) - IntVect::Unit;
  Box derivBox(a_phiCen.box());
  derivBox.grow(shrinkVect);

  FArrayBox tempLap(derivBox, a_phiAvg.nComp());
  tempLap.setVal(0.);
  Real factor = a_sgn * (1.0/24.0);
  for (int tanDir=0; tanDir <SpaceDim; tanDir++)
    {
      if (tanDir != a_dir)
        {
          // tempLap += factor * D^2_{tanDir}(a_phiCen)
          FORT_INCREMENTLAPLACIAN(CHF_FRA(tempLap),
                                  CHF_CONST_FRA(a_phiCen),
                                  CHF_BOX(derivBox),
                                  CHF_CONST_INT(tanDir),
                                  CHF_CONST_REAL(factor));
        } // end if tangential direction
    } // end loop over tangential directions
  Box bxIntersect = derivBox & a_phiAvg.box();
  a_phiAvg.copy(a_phiCen);
  a_phiAvg.plus(tempLap, bxIntersect, 0, 0, ncomp);
}


/// -----------------------------------------------------------
void fourthOrderAverage(LevelData<FArrayBox>& a_phi,
                        int a_sgn)
{
  CH_TIME("fourtOrderAverage");

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisPhi = a_phi[dit];

      fourthOrderAverageCell(thisPhi, a_sgn);
    } // end loop over boxes
}


/// -----------------------------------------------------------
void fourthOrderAverage(LevelData<FArrayBox>& a_phi,
                        const ProblemDomain&  a_domain,
                        int a_sgn)
{
  CH_TIME("fourtOrderAverage domain");

  const DisjointBoxLayout& layout = a_phi.disjointBoxLayout();
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisPhi = a_phi[dit];
      const Box& bx = layout[dit];

      fourthOrderAverageCell(thisPhi, a_domain, bx, a_sgn);
    } // end loop over boxes
}


/// -----------------------------------------------------------
void fourthOrderAverageCell(FArrayBox& a_phi,
                            int a_sgn)
{
  CH_TIME("fourthOrderAverageCell");

  // derivBox shrunk by 1
  Box derivBox(a_phi.box());
  derivBox.grow(-1);

  FArrayBox tempLap(derivBox, a_phi.nComp());
  tempLap.setVal(0.0);

  Real factor = a_sgn * (1.0/24.0);
  for (int dir=0; dir <SpaceDim; dir++)
    {
      FORT_INCREMENTLAPLACIAN(CHF_FRA(tempLap),
                              CHF_CONST_FRA(a_phi),
                              CHF_BOX(derivBox),
                              CHF_CONST_INT(dir),
                              CHF_CONST_REAL(factor));
    } // end loop over directions
  a_phi.plus(tempLap, derivBox, 0, 0, tempLap.nComp());
}


/// -----------------------------------------------------------
void fourthOrderAverageCell(FArrayBox& a_phi,
                            const ProblemDomain&  a_domain,
                            const Box&            a_bx,
                            int a_sgn)
{
  CH_TIME("fourthOrderAverageCell domain");

  CH_assert(a_phi.box().contains(a_bx));
  Box derivBox(a_bx);
  derivBox.grow(1);
  derivBox &= a_domain;
  CH_assert(a_phi.box().contains(derivBox));

  FArrayBox tempLap(a_bx, a_phi.nComp());
  tempLap.setVal(0.0);

  Real factor = a_sgn * (1.0/24.0);
  for (int dir=0; dir <SpaceDim; dir++)
    {
      Box derivDirBox(a_bx);
      if (!a_domain.isPeriodic(dir) &&
          (a_domain.domainBox().smallEnd(dir) == a_bx.smallEnd(dir)))
        {
          // one-sided derivative, low end
          Box loBox(a_bx);
          loBox.setRange(dir, a_bx.smallEnd(dir));
          FORT_INCREMENTLOSIDELAPLACIAN(CHF_FRA(tempLap),
                                        CHF_FRA(a_phi),
                                        CHF_BOX(loBox),
                                        CHF_INT(dir),
                                        CHF_REAL(factor));
          derivDirBox.growLo(dir, -1);
        }
      if (!a_domain.isPeriodic(dir) &&
          (a_domain.domainBox().bigEnd(dir) == a_bx.bigEnd(dir)))
        {
          // one-sided derivative, high end
          Box hiBox(a_bx);
          hiBox.setRange(dir, a_bx.bigEnd(dir));
          FORT_INCREMENTHISIDELAPLACIAN(CHF_FRA(tempLap),
                                        CHF_FRA(a_phi),
                                        CHF_BOX(hiBox),
                                        CHF_INT(dir),
                                        CHF_REAL(factor));
          derivDirBox.growHi(dir, -1);
        }
      FORT_INCREMENTLAPLACIAN(CHF_FRA(tempLap),
                              CHF_CONST_FRA(a_phi),
                              CHF_BOX(derivDirBox),
                              CHF_CONST_INT(dir),
                              CHF_CONST_REAL(factor));
    } // end loop over directions
  a_phi.plus(tempLap, a_bx, 0, 0, tempLap.nComp());
}


/// -----------------------------------------------------------
void fourthOrderCellToFace(LevelData<FluxBox>& a_faceData,
                           const LevelData<FArrayBox>& a_cellData)
{
  CH_TIME("fourthOrderCellToFace level");

  const DisjointBoxLayout& grids = a_faceData.getBoxes();

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      fourthOrderCellToFace(a_faceData[dit], a_cellData[dit]);
    } // end loop over grid boxes
}


/// -----------------------------------------------------------
void fourthOrderCellToFace(FluxBox& a_faceData,
                           const FArrayBox& a_cellData)
{
  CH_TIME("fourthOrderCellToFace");

  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& thisFaceDataDir = a_faceData[dir];
      // compute box over which we can do the averaging
      Box testBox(a_cellData.box());
      // need two cells in each direction for averaging
      testBox.grow(dir, -2);
      testBox.surroundingNodes(dir);
      testBox &= thisFaceDataDir.box();
      FORT_CELLTOFACE4THORDER(CHF_FRA(thisFaceDataDir),
                              CHF_CONST_FRA(a_cellData),
                              CHF_BOX(testBox),
                              CHF_INT(dir));
    }
}


/// -----------------------------------------------------------
void fourthOrderCellToFace(FArrayBox& a_faceData,
                           const FArrayBox& a_cellData,
                           const Box& a_faceBox,
                           int a_dir)
{
  CH_TIME("fourthOrderCellToFace box");

  CH_assert(a_faceData.box().contains(a_faceBox));
  Box testBox(a_faceBox);
  testBox.enclosedCells(a_dir);
  testBox.grow(a_dir, 2);
  CH_assert(a_cellData.box().contains(testBox));


  FORT_CELLTOFACE4THORDER(CHF_FRA(a_faceData),
                          CHF_CONST_FRA(a_cellData),
                          CHF_BOX(a_faceBox),
                          CHF_INT(a_dir));


}


/// -----------------------------------------------------------
void fourthOrderCellToFaceCenters(LevelData<FluxBox>& a_faceData,
                                  const LevelData<FArrayBox>& a_cellData)
{
  CH_TIME("fourthOrderCellToFaceCenters level");

  const DisjointBoxLayout& grids = a_faceData.getBoxes();

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      fourthOrderCellToFaceCenters(a_faceData[dit], a_cellData[dit]);
    } // end loop over grid boxes
}


/// -----------------------------------------------------------
void fourthOrderCellToFaceCenters(FluxBox& a_faceData,
                                  const FArrayBox& a_cellData)
{
  CH_TIME("fourthOrderCellToFaceCenters");

  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& thisFaceDataDir = a_faceData[dir];
      // compute box over which we can do the averaging
      Box testBox(a_cellData.box());
      // need two cells in each direction for averaging
      testBox.grow(dir, -2);
      testBox.surroundingNodes(dir);
      testBox &= thisFaceDataDir.box();
      FORT_CELLTOFACE4THORDERCENTERS(CHF_FRA(thisFaceDataDir),
                                     CHF_CONST_FRA(a_cellData),
                                     CHF_BOX(testBox),
                                     CHF_INT(dir));
    }
}


/// -----------------------------------------------------------
void fourthOrderCellToFaceCenters(FArrayBox& a_faceData,
                                  const FArrayBox& a_cellData,
                                  const Box& a_faceBox,
                                  int a_dir)
{
  CH_TIME("fourthOrderCellToFaceCenters box");

  CH_assert(a_faceData.box().contains(a_faceBox));
  Box testBox(a_faceBox);
  testBox.enclosedCells(a_dir);
  testBox.grow(a_dir, 2);
  CH_assert(a_cellData.box().contains(testBox));

  FORT_CELLTOFACE4THORDERCENTERS(CHF_FRA(a_faceData),
                                 CHF_CONST_FRA(a_cellData),
                                 CHF_BOX(a_faceBox),
                                 CHF_INT(a_dir));
}


/// -----------------------------------------------------------
// compute fourth-order approximation to u*v
void fourthOrderMultFace(FluxBox& a_uTimesV,
                         const FluxBox& a_u,
                         const FluxBox& a_v,
                         int a_sgn)
{
  CH_TIME("fourthOrderMultFace flux");

  int ncomp = a_u.nComp();
  CH_assert(a_v.nComp() == ncomp);
  CH_assert(a_uTimesV.nComp() == ncomp);

  // since dx's cancel, use 1 here
  Real fakeDx = 1.0;
  Real factor = a_sgn * (1.0/12.0); //

  for (int faceDir=0; faceDir<SpaceDim; faceDir++)
     {
        FArrayBox& thisUVdir = a_uTimesV[faceDir];
        const FArrayBox& thisUdir = a_u[faceDir];
        const FArrayBox& thisVdir = a_v[faceDir];

        // first multiply <u><v>
        Box intersectBox(thisUdir.box());
        intersectBox &= thisVdir.box();
        intersectBox &= thisUVdir.box();

        thisUVdir.copy(thisUdir, intersectBox);
        thisUVdir.mult(thisVdir, intersectBox, 0, 0, ncomp);

        // now increment with product of tangential gradients
        for (int tanDir =0; tanDir<SpaceDim; tanDir++)
           {
              if (tanDir != faceDir)
                 {
                    // compute box over which we can do this
                    Box gradIntersectBox(thisUdir.box());
                    gradIntersectBox &= thisVdir.box();
                    gradIntersectBox.grow(tanDir,-1);
                    gradIntersectBox &= thisUVdir.box();

                    FORT_INCREMENTGRADPROD(CHF_FRA(thisUVdir),
                                           CHF_CONST_FRA(thisUdir),
                                           CHF_CONST_FRA(thisVdir),
                                           CHF_BOX(gradIntersectBox),
                                           CHF_REAL(fakeDx),
                                           CHF_REAL(factor),
                                           CHF_INT(tanDir));

                 } // end if tangential direction
           } // end loop over tanDir
     } // end loop over faceDir
}


/// -----------------------------------------------------------
// compute fourth-order approximation to u*v
void fourthOrderMultFace(FArrayBox& a_uTimesV,
                         const FArrayBox& a_u,
                         const FArrayBox& a_uForDx,
                         const FArrayBox& a_v,
                         const FArrayBox& a_vForDx,
                         int a_faceDir,
                         int a_sgn)
{
  CH_TIME("fourthOrderMultFace");

  int ncomp = a_u.nComp();
  CH_assert(a_v.nComp() == ncomp);
  CH_assert(a_uTimesV.nComp() == ncomp);

  // since dx's cancel, use 1 here
  Real fakeDx = 1.0;
  Real factor = a_sgn * (1.0/12.0); //

  // first multiply <u><v>
  Box intersectBox(a_u.box());
  intersectBox &= a_v.box();
  intersectBox &= a_uTimesV.box();

  a_uTimesV.copy(a_u, intersectBox);
  a_uTimesV.mult(a_v, intersectBox, 0, 0, ncomp);

  // now increment with product of tangential gradients
  for (int tanDir =0; tanDir<SpaceDim; tanDir++)
    {
      if (tanDir != a_faceDir)
        {
          // compute box over which we can do this
          Box gradIntersectBox(a_u.box());
          gradIntersectBox &= a_v.box();
          gradIntersectBox.grow(tanDir,-1);
          gradIntersectBox &= a_uTimesV.box();

          FORT_INCREMENTGRADPROD(CHF_FRA(a_uTimesV),
                                 CHF_CONST_FRA(a_uForDx),
                                 CHF_CONST_FRA(a_vForDx),
                                 CHF_BOX(gradIntersectBox),
                                 CHF_REAL(fakeDx),
                                 CHF_REAL(factor),
                                 CHF_INT(tanDir));

        } // end if tangential direction
    } // end loop over tanDir
}

/// -----------------------------------------------------------
// compute fourth-order approximation to u*v
void fourthOrderMult(LevelData<FluxBox>& a_uTimesV,
                     const LevelData<FluxBox>& a_u,
                     const LevelData<FluxBox>& a_v,
                     int a_sgn)
{
  CH_TIME("fourthOrderMult");

  int ncomp = a_u.nComp();
  CH_assert(a_v.nComp() == ncomp);
  CH_assert(a_uTimesV.nComp() == ncomp);

  DataIterator dit = a_u.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
       fourthOrderMultFace(a_uTimesV[dit], a_u[dit], a_v[dit], a_sgn);
    }
}

/// -----------------------------------------------------------
// compute fourth-order approximation to u*v
void fourthOrderCellProd(LevelData<FArrayBox>& a_uTimesV,
                         const LevelData<FArrayBox>& a_u,
                         const LevelData<FArrayBox>& a_v)
{
  CH_TIME("fourthOrderCellProd level");

   int ncomp = a_u.nComp();
   CH_assert(a_v.nComp() == ncomp);
   CH_assert(a_uTimesV.nComp() == ncomp);

   DataIterator dit = a_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FArrayBox& thisUV = a_uTimesV[dit];
      const FArrayBox& thisU = a_u[dit];
      const FArrayBox& thisV = a_v[dit];

      // first multiply <u><v>
      Box intersectBox(thisU.box());
      intersectBox &= thisV.box();
      intersectBox &= thisUV.box();

      thisUV.copy(thisU, intersectBox);
      thisUV.mult(thisV, intersectBox, 0, 0, ncomp);

      // now increment with product of tangential gradients
      for (int dir=0; dir<SpaceDim; dir++)
      {
         // since dx's cancel, use 1 here
         Real fakeDx = 1.0;
         Real factor = 1.0/12.0;

         // compute box over which we can do this
         Box gradIntersectBox(thisU.box());
         gradIntersectBox &= thisV.box();
         gradIntersectBox.grow(dir,-1);
         gradIntersectBox &= thisUV.box();

         FORT_INCREMENTGRADPROD(CHF_FRA(thisUV),
                                CHF_CONST_FRA(thisU),
                                CHF_CONST_FRA(thisV),
                                CHF_BOX(gradIntersectBox),
                                CHF_REAL(fakeDx),
                                CHF_REAL(factor),
                                CHF_INT(dir));

      } // end loop over tanDir
   } // end loop over boxes
}

/*--------------------------------------------------------------------*/
//   Computes <FG> from <F> and <G> using 1-sided derivatives at
//   boundaries
/**  Gradients of \<F\> are determined using 1-sided derivatives at
 *   the domain boundary.  \<G\> is usually \<J\> and since it is
 *   known outside the boundary, centered differences are used
 *   everywhere.  You can override this with the a_oneSidedDG
 *   argument.
 *   \param[out] a_FG   \<FG\>
 *   \param[in]  a_F    \<F\>, may have multiple components
 *   \param[in]  a_G    \<G\>, only has 1 component
 *   \param[in]  a_box  Where to compute
 *   \param[in]  a_domain
 *                      The problem domain
 *   \param[in]  a_oneSidedDG
 *                      T - Compute d\<G\> using one-sided
 *                          differences near the boundary (same as
 *                          d\<F\>).
 *                      F - Use centered differences everywhere
 *                          (default)
 *   \note
 *   <ul>
 *     <li>
 *          \f[
 *             \langle FG\rangle = \langle F\rangle\langle G\rangle +
 *                \frac(h^2)(12)\left(
 *                   G_0(\langle F\rangle) \cdot
 *                   G_0(\langle G\rangle)\right)
 *          \f]
 *          as described in Fourth-order Chombo doc.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void fourthOrderCellProd(FArrayBox&           a_FG,
                         const FArrayBox&     a_F,
                         const FArrayBox&     a_G,
                         const Box&           a_box,
                         const ProblemDomain& a_domain,
                         const bool           a_oneSidedDG,
                         Interval             a_fixedDims)
{
  CH_TIME("fourthOrderCellProd domain");

  int nComp = a_FG.nComp();
  CH_assert(a_F.nComp() == nComp);
  CH_assert(a_G.nComp() == 1);
  // <FG> is only computed in the domain.  This assertion forces the user to
  // understand this.
  CH_assert(a_domain.contains(a_box));

  IntVect interpUnit = IntVect::Unit;
  for (int dir = a_fixedDims.begin(); dir <= a_fixedDims.end(); dir++)
    {
      interpUnit[dir] = 0;
    }
  // Grown box needed to compute gradients
  Box grownBox = grow(a_box, interpUnit);
  Box grownBoxInDom = grownBox & a_domain;

  CH_assert(a_FG.box().contains(a_box));
  CH_assert(a_F.box().contains(grownBoxInDom));
  if (a_oneSidedDG)
    {
      CH_assert(a_G.box().contains(grownBoxInDom));
    }
  else
    {
      CH_assert(a_G.box().contains(grownBox));
    }

//--Initially set <FG> = <F><G> on a_box.

  a_FG.copy(a_F, a_box);
  for (int iComp = 0; iComp != nComp; ++iComp)
    {
      a_FG.mult(a_G, a_box, 0, iComp);
    }

//--Add products of gradients to a_FG

  FArrayBox gradF(a_box, nComp);
  Real factor = 1.0/12.0;
  for (int dir = 0; dir != SpaceDim; ++dir)
    if ( !a_fixedDims.contains(dir) )
    {
      int hasLo = 0;
      int hasHi = 0;
      Box loBox(a_box);
      Box hiBox(a_box);
      Box centerBox(a_box);

      // Get the gradients of F
      if (!a_domain.isPeriodic(dir) &&
          (a_domain.domainBox().smallEnd(dir) == a_box.smallEnd(dir)))
        {
          // One-sided derivatives, low end
          hasLo = 1;
          loBox.setRange(dir, a_box.smallEnd(dir));
          FORT_UDIVLOSIDEGRAD(CHF_FRA(gradF),
                              CHF_CONST_FRA(a_F),
                              CHF_BOX(loBox),
                              CHF_CONST_INT(dir));
          centerBox.growLo(dir, -1);
        }
      if (!a_domain.isPeriodic(dir) &&
          (a_domain.domainBox().bigEnd(dir) == a_box.bigEnd(dir)))
        {
          // One-sided derivative, high end
          hasHi = 1;
          hiBox.setRange(dir, a_box.bigEnd(dir));
          FORT_UDIVHISIDEGRAD(CHF_FRA(gradF),
                              CHF_CONST_FRA(a_F),
                              CHF_BOX(hiBox),
                              CHF_CONST_INT(dir));
          centerBox.growHi(dir, -1);
        }
      // Fill remaining using centered differences
      FORT_UDIVCENTERGRAD(CHF_FRA(gradF),
                          CHF_CONST_FRA(a_F),
                          CHF_BOX(centerBox),
                          CHF_CONST_INT(dir));

      // If gradients of G are all to be centered, then ignore lo and hi boxes
      if (!a_oneSidedDG)
        {
          hasLo = 0;
          hasHi = 0;
          centerBox = a_box;
        }

      // Add products of gradients to a_FG.  Gradients of a_G are computed
      // within.
      FORT_INCREMENTGRADPRODDU(CHF_FRA(a_FG),
                               CHF_CONST_FRA(gradF),
                               CHF_CONST_FRA1(a_G, 0),
                               CHF_BOX(loBox),
                               CHF_CONST_INT(hasLo),
                               CHF_BOX(hiBox),
                               CHF_CONST_INT(hasHi),
                               CHF_BOX(centerBox),
                               CHF_BOX(a_box),
                               CHF_CONST_REAL(factor),
                               CHF_CONST_INT(dir));
    }
}

// this is just a basic wrapper around the INCREMENTGRADPROD fortran
void
incrementGradProduct(FArrayBox& a_gradProduct,
                     const FArrayBox& a_u,
                     const FArrayBox& a_v,
                     const Box& a_region,
                     Real a_dx,
                     Real a_factor,
                     int a_dir)
{
  CH_TIME("incrementGradProduct");

  FORT_INCREMENTGRADPROD(CHF_FRA(a_gradProduct),
                         CHF_CONST_FRA(a_u),
                         CHF_CONST_FRA(a_v),
                         CHF_BOX(a_region),
                         CHF_REAL(a_dx),
                         CHF_REAL(a_factor),
                         CHF_INT(a_dir));
}

// Fourth-order extrapolate cell-averaged interior data to two ghost cell layers at domain boundaries
void
fourthOrderCellExtrapAtDomainBdry(FArrayBox&           a_data,
                                  const ProblemDomain& a_domain,
                                  const Box&           a_interiorbox)
{
  CH_TIME("fourthOrderCellExtrapAtDomainBdry");

  const Box& bx= a_data.box();
  CH_assert(bx.contains(a_interiorbox));

  const Box& domBox = a_domain.domainBox();

  for (int dir=0; dir<SpaceDim; ++dir)
  {
     // Periodic BC's trump all other BCs -- do nothing in that case
     if (!a_domain.isPeriodic(dir))
     {
        for (int side=-1; side<2; side+=2)
        {
           Box dstbox;
           bool isbc = false;

           switch(side)
              {
              case -1:
                 // Handle the low side
                 isbc = (a_interiorbox.smallEnd(dir) == domBox.smallEnd(dir));
                 if (isbc)
                 {
                    dstbox = adjCellLo(a_interiorbox, dir, 2);
                 }
                 break;
              case 1:
                 // Handle high side
                 isbc = (a_interiorbox.bigEnd(dir) == domBox.bigEnd(dir));
                 if (isbc)
                 {
                    dstbox = adjCellHi(a_interiorbox, dir, 2);
                 }
                 break;
              }

           if (isbc)
           {
              CH_assert(a_interiorbox.size(dir)>=4);
              CH_assert(bx.contains(dstbox));
              FORT_CELLEXTRAP4THORDER(
                 CHF_CONST_INT(dir),
                 CHF_CONST_INT(side),
                 CHF_BOX(a_interiorbox),
                 CHF_BOX(dstbox),
                 CHF_FRA(a_data)
                 );
           }
        }
     }
  }
}

// Fill interior ghost cells and boundary ghost cell with cell-averaged data via extrapolation
void
fourthOrderCellGhostFill(LevelData<FArrayBox>& a_data,
                         const ProblemDomain&  a_domain)
{
  CH_TIME("fourthOrderCellGhostFill");

   // Fill interior ghost cells
   a_data.exchange();

   // Fill physical boundary ghost cells
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
      {
         fourthOrderCellExtrapAtDomainBdry(a_data[dit], a_domain, grids[dit]);
      } // end loop over boxes
}

// Fill width-two boxes on low and high sides of a_srcbox (in all directions)
// with cell-averaged values extrapolated to fourth-order from values in a_srcbox.
void
fourthOrderCellExtrap(FArrayBox& a_data,
                      const Box& a_srcbox)
{
  CH_TIME("fourthOrderCellExtrap");

   const Box& bx = a_data.box();
   CH_assert(bx.contains(a_srcbox));

   for (int dir=0; dir<SpaceDim; ++dir)
   {
      CH_assert(a_srcbox.size(dir)>=4);
      for (int side=-1; side<2; side+=2)
      {

         Box dstbox;
         if (side==-1)
         {
            dstbox = adjCellLo(a_srcbox, dir, 2);
         }
         else
         {
            dstbox = adjCellHi(a_srcbox, dir, 2);
         }

         CH_assert(bx.contains(dstbox));
         FORT_CELLEXTRAP4THORDER(
            CHF_CONST_INT(dir),
            CHF_CONST_INT(side),
            CHF_BOX(a_srcbox),
            CHF_BOX(dstbox),
            CHF_FRA(a_data)
            );
      }
   }
}

// Fill width-one boxes on low and high sides of a_srcbox (in all directions)
// with values extrapolated to second-order from values in a_srcbox.
void
secondOrderCellExtrap(FArrayBox& a_data,
                      const Box& a_srcbox)
{
  CH_TIME("secondOrderCellExtrap");

   for (int dir=0; dir<SpaceDim; ++dir)
     {
       secondOrderCellExtrap(a_data, a_srcbox, dir);
     }
}

// Fill width-one boxes on low and high sides of a_srcbox
// in specified direction,
// with values extrapolated to second-order from values in a_srcbox.
void
secondOrderCellExtrap(FArrayBox& a_data,
                      const Box& a_srcbox,
                      int a_dir)
{
  CH_TIME("secondOrderCellExtrap dir");

   const Box& bx = a_data.box();
   CH_assert(bx.contains(a_srcbox));

   CH_assert(a_srcbox.size(a_dir)>=3);
   for (SideIterator sit; sit.ok(); ++sit)
     {
       Side::LoHiSide side = sit();
       int sideSign = sign(side);
       Box dstbox = adjCellBox(a_srcbox, a_dir, side, 1);
       CH_assert(bx.contains(dstbox));
       FORT_SECOND_ORDER_EXTRAP(
                                CHF_CONST_INT(a_dir),
                                CHF_CONST_INT(sideSign),
                                CHF_BOX(a_srcbox),
                                CHF_BOX(dstbox),
                                CHF_FRA(a_data)
                                );
     }
}

// Second-order extrapolate at domain boundary (in all directions)
void
secondOrderCellExtrapAtDomainBdry(FArrayBox&           a_data,
                                  const Box&           a_interiorbox,
                                  const ProblemDomain& a_domain)
{
  CH_TIME("secondOrderCellExtrapAtDomainBdry");

   const Box& bx = a_data.box();
   CH_assert(bx.contains(a_interiorbox));
   const Box& domBox = a_domain.domainBox();

   for (int dir=0; dir<SpaceDim; ++dir)
   {
      if (!a_domain.isPeriodic(dir))
      {
         CH_assert(a_interiorbox.size(dir)>=3);

         for (int side=-1; side<2; side+=2)
         {
            Box dstbox;
            bool isbc = false;

            switch(side)
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(dir) == domBox.smallEnd(dir));
               if (isbc)
               {
                  dstbox = adjCellLo(a_interiorbox,dir,1);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(dir) == domBox.bigEnd(dir));
               if (isbc)
               {
                  dstbox = adjCellHi(a_interiorbox,dir,1);
               }
               break;
            }

            if (isbc)
            {
               CH_assert(bx.contains(dstbox));
               FORT_SECOND_ORDER_EXTRAP(
                 CHF_CONST_INT(dir),
                 CHF_CONST_INT(side),
                 CHF_BOX(a_interiorbox),
                 CHF_BOX(dstbox),
                 CHF_FRA(a_data)
                 );
            }
         }
      }
   }
}

// Second-order extrapolate in transverse directions
void
secondOrderTransExtrap(FArrayBox& a_data,
                       const int  a_dir,
                       const Box& a_srcbox)
{
  CH_TIME("secondOrderTransExtrap");

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_srcbox));

   for (int tdir=0; tdir<SpaceDim; ++tdir)
   {
      if (tdir != a_dir)
      {
         CH_assert(a_srcbox.size(tdir)>=3);

         for (int side=-1; side<2; side+=2)
         {
            Box dstbox;

            switch(side)
            {
            case -1:
               dstbox = adjCellLo(a_srcbox,tdir,1);
               break;
            case 1:
               dstbox = adjCellHi(a_srcbox,tdir,1);
               break;
            }

            CH_assert(bx.contains(dstbox));
            FORT_SECOND_ORDER_EXTRAP(
               CHF_CONST_INT(tdir),
               CHF_CONST_INT(side),
               CHF_BOX(a_srcbox),
               CHF_BOX(dstbox),
               CHF_FRA(a_data)
               );
         }
      }
   }
}

// Second-order extrapolate in transverse directions at domain boundary
void
secondOrderTransExtrapAtDomainBdry(FArrayBox&           a_data,
                                   const int            a_dir,
                                   const Box&           a_interiorbox,
                                   const ProblemDomain& a_domain)
{
  CH_TIME("secondOrderTransExtrapAtDomainBdry");

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));
   const Box& domBox = a_domain.domainBox();

   for (int tdir=0; tdir<SpaceDim; ++tdir)
   {
      if (tdir != a_dir && !a_domain.isPeriodic(tdir))
      {
         CH_assert(a_interiorbox.size(tdir)>=3);

         for (int side=-1; side<2; side+=2)
         {
            Box dstbox;
            bool isbc = false;

            switch(side)
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(tdir) == domBox.smallEnd(tdir));
               if (isbc)
               {
                  dstbox = adjCellLo(a_interiorbox,tdir,1);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(tdir) == domBox.bigEnd(tdir));
               if (isbc)
               {
                  dstbox = adjCellHi(a_interiorbox,tdir,1);
               }
               break;
            }

            if (isbc)
            {
               CH_assert(bx.contains(dstbox));
               FORT_SECOND_ORDER_EXTRAP(
                  CHF_CONST_INT(tdir),
                  CHF_CONST_INT(side),
                  CHF_BOX(a_interiorbox),
                  CHF_BOX(dstbox),
                  CHF_FRA(a_data)
                  );
            }
         }
      }
   }
}

// Second-order extrapolate in transverse directions at domain boundary
void
secondOrderTransExtrapFacesAtDomainBdry(FArrayBox&           a_data,
                                        const int            a_dir,
                                        const Box&           a_interiorbox,
                                        const ProblemDomain& a_domain)
{
  CH_TIME("secondOrderTransExtrapFacesAtDomainBdry");

   const Box& bx= a_data.box();
   const Box& domBox = a_domain.domainBox();
   for (int tdir=0; tdir<SpaceDim; ++tdir)
   {
      if (tdir != a_dir && !a_domain.isPeriodic(tdir))
      {
         CH_assert(a_interiorbox.size(tdir)>=3);

         for (int side=-1; side<2; side+=2)
         {
            // dstbox will contain the a_dir-faces of a_data to be changed
            Box dstbox;
            bool isbc = false;

            switch(side)
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(tdir) == domBox.smallEnd(tdir));
               if (isbc)
               {
                  dstbox = adjCellLo(a_interiorbox,tdir,1);
                  dstbox.surroundingNodes(a_dir);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(tdir) == domBox.bigEnd(tdir));
               if (isbc)
               {
                  dstbox = adjCellHi(a_interiorbox,tdir,1);
                  dstbox.surroundingNodes(a_dir);
               }
               break;
            }
            // dstbox contains the a_dir-faces of a_data to be changed
            if (isbc)
            {
               CH_assert(bx.contains(dstbox));
               FORT_SECOND_ORDER_EXTRAP(
                  CHF_CONST_INT(tdir),
                  CHF_CONST_INT(side),
                  CHF_BOX(a_interiorbox),
                  CHF_BOX(dstbox),
                  CHF_FRA(a_data)
                  );
            }
         }
      }
   }
}

// At domain boundaries, fill width-one boxes on low and high sides of a_srcbox
// (in all directions) with values extrapolated to zeroth-order from values in
// a_srcbox.
void
zerothOrderCellExtrapAtDomainBdry(FArrayBox&           a_data,
                                  const Box&           a_interiorbox,
                                  const ProblemDomain& a_domain)
{
  CH_TIME("zerothOrderCelLExtrapAtDomainBdry");

   const Box& bx = a_data.box();
   CH_assert(bx.contains(grow(a_interiorbox, 1)));
   const Box& domBox = a_domain.domainBox();

   for (int dir=0; dir<SpaceDim; ++dir)
   {
      if (!a_domain.isPeriodic(dir))
      {
        // Handle the low side
        if (a_interiorbox.smallEnd(dir) == domBox.smallEnd(dir))
          {
            Box dstbox = adjCellLo(a_interiorbox,dir,1);
            Box srcbox(dstbox);
            srcbox.shift(dir, 1);
            a_data.copy(srcbox,
                        a_data.interval(),
                        dstbox,
                        a_data,
                        a_data.interval());
          }
        // Handle high side
        if (a_interiorbox.bigEnd(dir) == domBox.bigEnd(dir))
          {
            Box dstbox = adjCellHi(a_interiorbox,dir,1);
            Box srcbox(dstbox);
            srcbox.shift(dir, -1);
            a_data.copy(srcbox,
                        a_data.interval(),
                        dstbox,
                        a_data,
                        a_data.interval());
          }
      }
   }
}

/*--------------------------------------------------------------------*/
//   Determine transverse components of N by averaging from nearby
//   faces
/**  \param[out] a_Nface
 *                      Metrics on the faces in 'a_faceDir'
 *   \param[in]  a_NtransverseFace
 *                      Contains normal components of N on the faces
 *                      in the transverse directions
 *   \param[in]  a_box  Transverse terms will be update on all faces
 *                      in 'a_faceDir'
 *   \param[in]  a_faceDir
 *                      Face in this direction are updated...
 *   \param[in]  a_transverseDir
 *                      ... with transverse terms in this direction
 *                                                                     \verbatim
 *   x======o======o======x     o - bounds a_box
 *   |      H      H      |     x - bounds grownCellBox
 *   |      H      H      |     H - faces to update with transverse
 *   |      H      H      |         terms
 *   x======o======o======x     = - faces from which to average
 *                                  transverse terms                \endverbatim
 *  \note
 *  <ul>
 *    <li> The normal components in 'a_NtransverseFace' must be known
 *         in one cell outside a_box for averaging (see above picture)
 *         E.g., a_NtransverseFace.box() >=
 *         surroundingNodes(a_box.grow(a_faceDir, 1), a_transverseDir)
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
averageTransverseN(FArrayBox&       a_Nface,
                   const FArrayBox& a_NtransverseFace,
                   const Box&       a_box,
                   const int        a_faceDir,
                   const int        a_transverseDir)
{
  CH_TIME("averageTransverseN");

  CH_assert(a_faceDir != a_transverseDir);
  const Box faceBox = surroundingNodes(a_box, a_faceDir);
  CH_assert(a_Nface.box().contains(faceBox));
  Box grownCellBox(a_box);
  grownCellBox.grow(a_faceDir, 1);
  CH_assert(a_NtransverseFace.box().contains(
              surroundingNodes(grownCellBox, a_transverseDir)));

  FORT_AVERAGETRANSVERSEN(
    CHF_FRA(a_Nface),
    CHF_CONST_FRA(a_NtransverseFace),
    CHF_BOX(faceBox),
    CHF_BOX(grownCellBox),
    CHF_CONST_INT(a_faceDir),
    CHF_CONST_INT(a_transverseDir));
}

/*--------------------------------------------------------------------*/

void cellFGToCellF(FArrayBox& a_F,
                   const FArrayBox& a_FG, // may have multiple components
                   const FArrayBox& a_G,  // 1 component
                   const Box& a_box) // where we want to compute F
{
  CH_TIME("cellFGToCellF");

  // return <F> = <FG>/<G> - h^2/12 * (grad (<FG>/<G>)) dot (grad <G>) / <G>
  // as described in ESL paper:  Colella, Dorr, Hittinger, and Martin,
  // "High-order, finite-volume methods in mapped coordinates"

  int ncomp = a_F.nComp();
  CH_assert(a_FG.nComp() == ncomp);
  CH_assert(a_G.nComp() == 1);

  // grown box needed to compute gradients
  Box grownBox(a_box);
  grownBox.grow(1);

  CH_assert(a_F.box().contains(a_box));
  CH_assert(a_FG.box().contains(grownBox));
  CH_assert(a_G.box().contains(grownBox));


  // Set Ginv = 1.0 / a_G
  FArrayBox Ginv(grownBox, 1);
  Ginv.copy(a_G);
  Ginv.invert(1.0);

  FArrayBox temp(grownBox, ncomp);

  // Initially set <F> = <FG>/<G> on intersectBox.
  temp.copy(a_FG, grownBox);
  for (int comp = 0; comp < ncomp; comp++)
    {
      // Multiply a_F[comp] by Ginv[0].
      temp.mult(Ginv, grownBox, 0, comp);
    }

  // Set gradProduct = factor * grad(a_F) dot grad(a_G)
  //                 = -1/12 * ( h*grad(<FG>/<G>) ) dot ( h*grad(<G>) )

  a_F.setVal(0., a_box, 0, ncomp);
  // since dx's cancel, use 1 here
  Real fakeDx = 1.0;
  Real factor = -1.0/12.0; // (dfm - 3/22/09)--WAS 12 instead of 24 -- No, change back!
  for (int dir=0; dir<SpaceDim; dir++)
    {
      // This function is in fourthOrderMappedGrids/util/FourthOrderUtil.cpp
      // sets gradProduct += factor * d(a_F)/dx[dir] * d(a_G)/dx[dir]
      incrementGradProduct(a_F,
                           temp,
                           a_G,
                           a_box,
                           fakeDx,
                           factor,
                           dir);
    }

  for (int comp = 0; comp < ncomp; comp++)
    {
      // Multiply gradProduct[comp] by Ginv[0].
      a_F.mult(Ginv, a_box, 0, comp);
    }
  // Now F = gradProduct = factor * grad(a_F) dot grad(a_G) / a_G

  a_F.plus(temp, a_box, 0, 0, ncomp);
}

/*--------------------------------------------------------------------*/
//   Computes <F> from <FG> and <G> using 1-sided derivatives at
//   boundaries
/**  Gradients of \<F\> are determined from
 *   \f$\langle\bar(F)\rangle = \langle FG\rangle/\langle G\rangle\f$
 *   and with 1-sided derivatives at the domain boundary.  \<G\> is
 *   usually \<J\> and since it is known outside the boundary,
 *   centered differences are usually used everywhere.  You can
 *   override this with the a_oneSidedDG argument.
 *   \param[out] a_F    \<F\>
 *   \param[in]  a_FG   \<FG\>, may have multiple components
 *   \param[in]  a_G    \<G\>, only has 1 component
 *   \param[in]  a_box  Where to compute
 *   \param[in]  a_domain
 *                      The problem domain
 *   \param[in]  a_oneSidedDG
 *                      T - Compute d\<G\> using one-sided
 *                          differences near the boundary (same as
 *                          d\<F\>).
 *                      F - Use centered differences everywhere
 *                          (default)
 *   \note
 *   <ul>
 *     <li>
 *          \f[
 *             \langle F\rangle = \frac(1)(\langle G\rangle)\left[
 *                \langle FG\rangle - \frac(h^2)(12)\left(
 *                   G_0(\langle\bar(F)\rangle) \cdot
 *                   G_0(\langle G\rangle)\right)\right]
 *          \f]
 *          as described in ESL paper:  Colella, Dorr, Hittinger,
 *          and Martin, "High-order, finite-volume methods in mapped
 *          coordinates"
 *   </ul>
 *//*-----------------------------------------------------------------*/

void cellFGToCellF(FArrayBox&           a_F,
                   const FArrayBox&     a_FG,
                   const FArrayBox&     a_G,
                   const Box&           a_box,
                   const ProblemDomain& a_domain,
                   const bool           a_oneSidedDG,
                   Interval             a_fixedDims)
{
  CH_TIME("cellFGToCellF");

  int nComp = a_F.nComp();
  CH_assert(a_FG.nComp() == nComp);
  CH_assert(a_G.nComp() == 1);
  // <F> is only computed in the domain.  This assertion forces the user to
  // understand this.
  CH_assert(a_domain.contains(a_box));

  // Grown box needed to compute gradients
  Box grownBox   = grow(a_box, 1);
  Box grownBoxInDom = grownBox & a_domain;

  CH_assert(a_F.box().contains(a_box));
  CH_assert(a_FG.box().contains(grownBoxInDom));
  if (a_oneSidedDG)
    {
      CH_assert(a_G.box().contains(grownBoxInDom));
    }
  else
    {
      CH_assert(a_G.box().contains(grownBox));
    }

  // Set Ginv = 1.0 / a_G
  FArrayBox Ginv(grownBoxInDom, 1);
  Ginv.copy(a_G);
  Ginv.invert(1.0);

//--Initially set <\bar{F}> = <FG>/<G> on grownBoxInDom.

  FArrayBox barF(grownBoxInDom, nComp);
  barF.copy(a_FG, grownBoxInDom);
  for (int iComp = 0; iComp != nComp; ++iComp)
    {
      // Multiply a_F[iComp] by Ginv[0].
      barF.mult(Ginv, grownBoxInDom, 0, iComp);
    }

  // Start a_F with <FG>
  a_F.copy(a_FG, a_box);

//--Add products of gradients to a_F

  FArrayBox gradBarF(a_box, nComp);
  Real factor = -1.0/12.0;
  for (int dir = 0; dir != SpaceDim; ++dir)
    if ( !a_fixedDims.contains(dir) )
    {
      int hasLo = 0;
      int hasHi = 0;
      Box loBox(a_box);
      Box hiBox(a_box);
      Box centerBox(a_box);

      // Get the gradients of barF
      if (!a_domain.isPeriodic(dir) &&
          (a_domain.domainBox().smallEnd(dir) == a_box.smallEnd(dir)))
        {
          // One-sided derivatives, low end
          hasLo = 1;
          loBox.setRange(dir, a_box.smallEnd(dir));
          FORT_UDIVLOSIDEGRAD(CHF_FRA(gradBarF),
                              CHF_CONST_FRA(barF),
                              CHF_BOX(loBox),
                              CHF_CONST_INT(dir));
          centerBox.growLo(dir, -1);
        }
      if (!a_domain.isPeriodic(dir) &&
          (a_domain.domainBox().bigEnd(dir) == a_box.bigEnd(dir)))
        {
          // One-sided derivative, high end
          hasHi = 1;
          hiBox.setRange(dir, a_box.bigEnd(dir));
          FORT_UDIVHISIDEGRAD(CHF_FRA(gradBarF),
                              CHF_CONST_FRA(barF),
                              CHF_BOX(hiBox),
                              CHF_CONST_INT(dir));
          centerBox.growHi(dir, -1);
        }
      // Fill remaining using centered differences
      FORT_UDIVCENTERGRAD(CHF_FRA(gradBarF),
                          CHF_CONST_FRA(barF),
                          CHF_BOX(centerBox),
                          CHF_CONST_INT(dir));

      // If gradients of G are all to be centered, then ignore lo and hi boxes
      if (!a_oneSidedDG)
        {
          hasLo = 0;
          hasHi = 0;
          centerBox = a_box;
        }

      // Add products of gradients to a_F.  Gradients of a_G are computed
      // within.
      FORT_INCREMENTGRADPRODDU(CHF_FRA(a_F),
                               CHF_CONST_FRA(gradBarF),
                               CHF_CONST_FRA1(a_G, 0),
                               CHF_BOX(loBox),
                               CHF_CONST_INT(hasLo),
                               CHF_BOX(hiBox),
                               CHF_CONST_INT(hasHi),
                               CHF_BOX(centerBox),
                               CHF_BOX(a_box),
                               CHF_CONST_REAL(factor),
                               CHF_CONST_INT(dir));
    }

//--Finish a_F by dividing by a_G

  for (int iComp = 0; iComp != nComp; ++iComp)
    {
      // Multiply a_F[iComp] by Ginv[0].
      a_F.mult(Ginv, a_box, 0, iComp);
    }
}

#include "NamespaceFooter.H"
