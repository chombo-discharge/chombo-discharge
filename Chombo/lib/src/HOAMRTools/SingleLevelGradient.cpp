#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// SingleLevelGradient.cpp
// Dan Martin, Fri, Jan 14, 2000


#include "SingleLevelGradient.H"
#include "GradientF_F.H"
#include "EdgeToCell.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "Mask.H"

#include "NamespaceHeader.H"

#define NUMGRADGROW 1
// ----------------------------------------------------------





// ----------------------------------------------------------
void
SingleLevelGradient::levelGradientMAC(LevelData<FluxBox>& a_edgeGrad,
                                      const LevelData<FArrayBox>& a_phi,
                                      const Real a_dx)

{

  CH_assert(a_edgeGrad.nComp() >= a_phi.nComp());

#ifdef NEW_GRADIENT
  // build IVS to handle grid detals
  LayoutData<IntVectSet> gridIVS;
  int nGrow = NUMGRADGROW;
  createGridIVS(gridIVS, a_edgeGrad.getBoxes(), nGrow);
#endif

  DataIterator dit = a_phi.dataIterator();

  if (a_edgeGrad.nComp() == a_phi.nComp())
    {
      // loop over boxes and compute gradient
      for (dit.reset(); dit.ok(); ++dit)
        {
          const FArrayBox& thisPhi = a_phi[dit()];
          FluxBox& thisEdgeGrad = a_edgeGrad[dit()];

          // loop over directions
          for (int dir = 0; dir<SpaceDim; dir++)
            {
              FArrayBox& edgeGradDirFab = thisEdgeGrad[dir];
              const Box& edgeBox = edgeGradDirFab.box();
              // for one-component gradient, only do normal gradients
              int edgeDir = dir;
#ifdef NEW_GRADIENT
              singleBoxMacGrad(edgeGradDirFab,
                               thisPhi,
                               0, 0, thisPhi.nComp(),
                               edgeBox, a_dx, dir, edgeDir,
                               gridIVS[dit()]);

#else
              // corner boxes identify which cells need special
              // stencil for transverse derivatives
              // "bottom", "top", "left", and "right" are relative
              // to edgeDir=0, dir=1
              Box BLcorner(edgeBox);
              Box BRcorner(edgeBox);
              Box ULcorner(edgeBox);
              Box URcorner(edgeBox);

              // now loop over directions and set box edge extents
              for (int tempDir=0; tempDir<SpaceDim; ++tempDir)
                {
                  if ((tempDir==edgeDir))
                    {
                      BLcorner.setBig(tempDir,edgeBox.smallEnd(tempDir));
                      ULcorner.setBig(tempDir,edgeBox.smallEnd(tempDir));
                      BRcorner.setSmall(tempDir,edgeBox.bigEnd(tempDir));
                      URcorner.setSmall(tempDir,edgeBox.bigEnd(tempDir));
                    }
                  else if (tempDir == dir)
                    {
                      BLcorner.setBig(tempDir,edgeBox.smallEnd(tempDir));
                      BRcorner.setBig(tempDir,edgeBox.smallEnd(tempDir));
                      ULcorner.setSmall(tempDir,edgeBox.bigEnd(tempDir));
                      URcorner.setSmall(tempDir,edgeBox.bigEnd(tempDir));
                    }
                }

              // loop over components
              for (int comp=0; comp<a_phi.nComp(); comp++)
                {
                  // now compute gradient in direction dir on this grid
                  FORT_MACGRAD(CHF_FRA1(edgeGradDirFab,comp),
                               CHF_CONST_FRA1(thisPhi,comp),
                               CHF_BOX(edgeBox),
                               CHF_CONST_REAL(a_dx),
                               CHF_INT(dir),
                               CHF_INT(edgeDir),
                               CHF_BOX(BLcorner),
                               CHF_BOX(BRcorner),
                               CHF_BOX(ULcorner),
                               CHF_BOX(URcorner));

                }
#endif
            } // end loop over dir
        } // end loop over grids

    } // end if only doing normal dir
  else if (a_edgeGrad.nComp() == SpaceDim*(a_phi.nComp()))
    {
    // multicomponent gradPhi means that we also need to do
    // transverse directions.
      DataIterator dit = a_phi.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          const FArrayBox& thisPhi = a_phi[dit()];
          FluxBox& thisEdgeGrad = a_edgeGrad[dit()];

          // loop over edges
          for (int edgeDir=0; edgeDir<SpaceDim; edgeDir++)
            {
              FArrayBox& thisEdgeGradDirFab = thisEdgeGrad[edgeDir];
              const Box& edgeBox = thisEdgeGradDirFab.box();

              // loop over component directions in edgeGrad
              // (this will be direction of gradient)
              for (int dir=0; dir<SpaceDim; dir++)
                {
#ifdef NEW_GRADIENT

#else
                  // corner boxes identify which cells need special
                  // stencil for transverse derivatives
                  // "bottom", "top", "left", and "right" are relative
                  // to edgeDir=0, dir=1
                  Box BLcorner(edgeBox);
                  Box BRcorner(edgeBox);
                  Box ULcorner(edgeBox);
                  Box URcorner(edgeBox);

                  // now loop over directions and set box edge extents
                  for (int tempDir=0; tempDir<SpaceDim; ++tempDir)
                    {
                      if ((tempDir==edgeDir))
                        {
                          BLcorner.setBig(tempDir,edgeBox.smallEnd(tempDir));
                          ULcorner.setBig(tempDir,edgeBox.smallEnd(tempDir));
                          BRcorner.setSmall(tempDir,edgeBox.bigEnd(tempDir));
                          URcorner.setSmall(tempDir,edgeBox.bigEnd(tempDir));
                        }
                      else if (tempDir == dir)
                        {
                          BLcorner.setBig(tempDir,edgeBox.smallEnd(tempDir));
                          BRcorner.setBig(tempDir,edgeBox.smallEnd(tempDir));
                          ULcorner.setSmall(tempDir,edgeBox.bigEnd(tempDir));
                          URcorner.setSmall(tempDir,edgeBox.bigEnd(tempDir));
                        }
                    }
#endif
                  // loop over components in phi
                  for (int phiComp=0; phiComp<a_phi.nComp(); ++phiComp)
                    {
                      int gradComp = phiComp*SpaceDim + dir;

#ifdef NEW_GRADIENT
                      singleBoxMacGrad(thisEdgeGradDirFab, thisPhi,
                                       gradComp, phiComp, 1,
                                       edgeBox, a_dx, dir, edgeDir,
                                       gridIVS[dit()]);
#else
                      FORT_MACGRAD(CHF_FRA1(thisEdgeGradDirFab,gradComp),
                                   CHF_CONST_FRA1(thisPhi, phiComp),
                                   CHF_BOX(edgeBox),
                                   CHF_CONST_REAL(a_dx),
                                   CHF_INT(dir),
                                   CHF_INT(edgeDir),
                                   CHF_BOX(BLcorner),
                                   CHF_BOX(BRcorner),
                                   CHF_BOX(ULcorner),
                                   CHF_BOX(URcorner));
#endif

                    } // end loop over components
                } // end loop over gradient directions
            } // end loop over edge orientation

        } // end loop over grids
    }
  else
    {
      // bad number of components in either phi or gradPhi
      MayDay::Error("levelGradientMAC: bad number of components!");
    }

}

// ----------------------------------------------------------
void
SingleLevelGradient::levelGradientCC(LevelData<FArrayBox>& a_grad,
                                     const LevelData<FArrayBox>& a_phi,
                                     const Real a_dx)
{
  // in this case, assume that all BC's have already been set
  // can just loop over grids and call gradCC subroutine
  // directly
  DataIterator dit = a_grad.dataIterator();
  const DisjointBoxLayout& grids = a_grad.getBoxes();

  for (dit.reset(); dit.ok(); ++dit)
    {
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FORT_GRADCC(CHF_FRA1(a_grad[dit()],dir),
                      CHF_CONST_FRA1(a_phi[dit()],0),
                      CHF_BOX(grids[dit()]),
                      CHF_CONST_REAL(a_dx),
                      CHF_INT(dir));
        }
    }
}

#ifdef NEW_GRADIENT
// utility function to do the actual computin' (to reduce code duplication)
void
SingleLevelGradient::singleBoxMacGrad(FArrayBox& a_gradFab,
                                      const FArrayBox& a_phiFab,
                                      int a_gradComp,
                                      int a_phiComp,
                                      int a_numComp,
                                      const Box& a_edgeBox,
                                      Real a_dx,
                                      int a_dir,
                                      int a_edgeDir,
                                      const IntVectSet& a_gridIVS)
{

  // corner boxes identify which cells need special
  // stencil for transverse derivatives
  // "bottom", "top", "left", and "right" are relative
  // to edgeDir=0, dir=1
  Box BLcorner(a_edgeBox);
  Box BRcorner(a_edgeBox);
  Box ULcorner(a_edgeBox);
  Box URcorner(a_edgeBox);

  if (a_edgeDir != a_dir)
    {
      // now loop over directions and set box edge extents
      for (int tempDir=0; tempDir<SpaceDim; ++tempDir)
        {
          if ((tempDir==a_edgeDir))
            {
              BLcorner.setBig(tempDir,a_edgeBox.smallEnd(tempDir));
              ULcorner.setBig(tempDir,a_edgeBox.smallEnd(tempDir));
              BRcorner.setSmall(tempDir,a_edgeBox.bigEnd(tempDir));
              URcorner.setSmall(tempDir,a_edgeBox.bigEnd(tempDir));
            }
          else if (tempDir == a_dir)
            {
              BLcorner.setBig(tempDir,a_edgeBox.smallEnd(tempDir));
              BRcorner.setBig(tempDir,a_edgeBox.smallEnd(tempDir));
              ULcorner.setSmall(tempDir,a_edgeBox.bigEnd(tempDir));
              URcorner.setSmall(tempDir,a_edgeBox.bigEnd(tempDir));
            }
        }
    }

  // loop over components
  int phiComp = a_phiComp;
  for (int comp=a_gradComp; comp< a_gradComp+a_numComp; comp++)
    {
      // now compute gradient in direction dir on this grid
      FORT_NEWMACGRAD(CHF_FRA1(a_gradFab,comp),
                      CHF_CONST_FRA1(a_phiFab,phiComp),
                      CHF_BOX(a_edgeBox),
                      CHF_CONST_REAL(a_dx),
                      CHF_INT(a_dir),
                      CHF_INT(a_edgeDir));

      // now do corner cells specially (replace gradient values computed
      // in fortran with one-sided gradients where necessary)
      if (a_dir != a_edgeDir)
        {
          Real factor = 0.25/a_dx;
          IntVect dirBasis = BASISV(a_dir);
          IntVect edgeBasis = BASISV(a_edgeDir);

          // "lower left"
          {
            BoxIterator bit(BLcorner);
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect thisIV = bit();
                IntVect cornerVect = thisIV -dirBasis - edgeBasis;
                if (!a_gridIVS.contains(cornerVect))
                  {
                    a_gradFab(thisIV, comp) = factor*(3.0*(a_phiFab(thisIV+dirBasis, phiComp)
                                                          -a_phiFab(thisIV-dirBasis, phiComp))
                                                      -( a_phiFab(thisIV+dirBasis+edgeBasis, phiComp)
                                                         -a_phiFab(thisIV-dirBasis+edgeBasis, phiComp)));
                  }
              } // end loop over cells in BLcorner
          }  // end BLcorner context



          // "lower right"
          {
            BoxIterator bit(BRcorner);
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect thisIV = bit();
                IntVect cornerVect = thisIV -dirBasis;
                if (!a_gridIVS.contains(cornerVect))
                  {
                    a_gradFab(thisIV, comp) = factor*(3.0*(a_phiFab(thisIV+dirBasis-edgeBasis, phiComp)
                                                          -a_phiFab(thisIV-dirBasis-edgeBasis, phiComp))
                                                      -( a_phiFab(thisIV+dirBasis-2*edgeBasis, phiComp)
                                                         -a_phiFab(thisIV-dirBasis-2*edgeBasis, phiComp)));
                  }
              } // end loop over cells in BRcorner
          }  // end BRcorner context


          // "upper left"
          {
            BoxIterator bit(ULcorner);
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect thisIV = bit();
                IntVect cornerVect = thisIV +dirBasis - edgeBasis;
                if (!a_gridIVS.contains(cornerVect))
                  {
                    a_gradFab(thisIV, comp) = factor*(3.0*(a_phiFab(thisIV+dirBasis, phiComp)
                                                          -a_phiFab(thisIV-dirBasis, phiComp))
                                                      -( a_phiFab(thisIV+dirBasis+edgeBasis, phiComp)
                                                       -a_phiFab(thisIV-dirBasis+edgeBasis, phiComp)));
                  }
              } // end loop over cells in ULcorner
          }  // end ULcorner context



          // "upper right"
          {
            BoxIterator bit(URcorner);
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect thisIV = bit();
                IntVect cornerVect = thisIV +dirBasis;
                if (!a_gridIVS.contains(cornerVect))
                  {
                    a_gradFab(thisIV, comp) = factor*(3.0*(a_phiFab(thisIV+dirBasis-edgeBasis, phiComp)
                                                          -a_phiFab(thisIV-dirBasis-edgeBasis, phiComp))
                                                      -( a_phiFab(thisIV+dirBasis-2*edgeBasis, phiComp)
                                                         -a_phiFab(thisIV-dirBasis-2*edgeBasis, phiComp)));
                  }
              } // end loop over cells in URcorner
          }  // end URcorner context
        }
      phiComp++;
    }
}


void
SingleLevelGradient::createGridIVS(LayoutData<IntVectSet>& a_gridsIVS,
                                   const DisjointBoxLayout& a_grids,
                                   const int a_nGrow)
{

  a_gridsIVS.define(a_grids);
  DataIterator dit = a_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      createGridIVS(a_gridsIVS[dit()],
                    a_grids,
                    a_grids[dit()],
                    a_nGrow);
    }
}


// simple little function to turn DBL->IVS in region around
// a given box.
void
SingleLevelGradient::createGridIVS(IntVectSet& a_gridIVS,
                        const DisjointBoxLayout& a_grids,
                        const Box& a_localBox,
                        const int a_nGrow)
{
  LayoutIterator lit = a_grids.layoutIterator();
  // if this isn't already empty, then make it empty
  if (!a_gridIVS.isEmpty())
    {
      a_gridIVS.makeEmpty();
    }

  // this is the neigborhood around thisBox in which we're interested.
  Box checkBox(a_localBox);
  checkBox.grow(a_nGrow);

  bool done = false;
    for (lit.begin(); lit.ok() && !done; ++lit)
    {
      const Box& thisBox = a_grids.get(lit());
      Box intersectBox(thisBox);
      intersectBox &= checkBox;
      if (!intersectBox.isEmpty())
        {
          a_gridIVS |= intersectBox;
        }

      // to make sure we do boundaries correctly, also include
      // adjacent cells
      // also check for whether we can stop due to ordering
      done = true;
      for (int dir=0; dir<SpaceDim; dir++)
        {
          Box lobox = adjCellLo(thisBox, dir, 1);
          Box hibox = adjCellHi(thisBox, dir, 1);

          lobox &= checkBox;
          hibox &= checkBox;
          if (!lobox.isEmpty())
            {
              a_gridIVS |= lobox;
            }
          if (!hibox.isEmpty())
            {
              a_gridIVS |= hibox;
            }
          if (lobox.smallEnd() <= checkBox.bigEnd())
            {
              done = false;
            }
        }
    }
}


#endif

#include "NamespaceFooter.H"
