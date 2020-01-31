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

#include "LoHiSide.H"
#include "BoxIterator.H"
#include "HOExtrapF_F.H"
#include "BasicIBC.H"

#include "UsingNamespace.H"

  /// compute transverse flux boundary conditions in ghost-cell faces
  /** In the mapped-grid case, where the transverse derivative dN/dx_j
      doesn't vanish at the boundary, we need values for the fluxes on
      ghost-cell faces transverse to the boundary, in order to be able
      to compute dFlux/dx_j to dot with the metric-term gradients.
      In other words, in a (0:N-1,0:N-1) domain, fluxes on x-faces will
      need to be set for all faces where j=-1 and N.
      I think we can generally use second-order extrapolation, which
      produces O(h) gradients.  (dropping one order of accuracy at the
      boundary)
  */
void
BasicIBC::fluxBC(LevelData<FluxBox>& a_flux,
                 const ProblemDomain& a_domain,
                 const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                 Real a_dx,
                 Real a_time)
{
  const Box& domainBox = a_domain.domainBox();

  // I think we can get away with 2nd-order extrapolation here
  DataIterator dit = a_flux.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisFlux = a_flux[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisFluxDir = thisFlux[dir];
          for (int transDir = 0; transDir<SpaceDim; transDir++)
            {
              // do extrapolation in transverse directions
              if ((transDir != dir) && (!a_domain.isPeriodic(dir)))
                {

                  // do lo end
                  Box loGhostBox = adjCellLo(domainBox, transDir, 1);
                  loGhostBox.surroundingNodes(dir);
                  loGhostBox &= thisFluxDir.box();
                  int hiLo = -1;
                  FORT_QUADEXTRAP(CHF_FRA(thisFluxDir),
                                  CHF_BOX(loGhostBox),
                                  CHF_INT(transDir),
                                  CHF_INT(hiLo));

                  // now do hi end
                  Box hiGhostBox = adjCellHi(domainBox, transDir, 1);
                  hiGhostBox.surroundingNodes(dir);
                  hiGhostBox &= thisFluxDir.box();
                  hiLo = 1;
                  FORT_QUADEXTRAP(CHF_FRA(thisFluxDir),
                                  CHF_BOX(hiGhostBox),
                                  CHF_INT(transDir),
                                  CHF_INT(hiLo));
                }
            } // end loop over transverse directions
        } // end loop over face directions
    } // end loop over boxes
}
