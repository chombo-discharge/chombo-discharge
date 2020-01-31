#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "zalesakLimiter.H"
#include "SingleLevelDivergence.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

using std::min;

/// apply Zalesak limiter to enforce minVal
/** a_flux comes in with high-order flux and is modified in place using
    the low-order flux
*/
void
applyZalesakLimiter(LevelData<FluxBox>& a_flux,
                    const LevelData<FluxBox>& a_lowOrderFlux,
                    const LevelData<FArrayBox>& a_oldPhi,
                    Real a_minVal,
                    Real a_dx,
                    Real a_dt,
                    bool a_lowOrderFluxesOnly)
{
  // first test - just use low-order fluxes
  if (a_lowOrderFluxesOnly)
    {
      a_lowOrderFlux.copyTo(a_flux);
    }
  else
    {
      const DisjointBoxLayout& grids = a_flux.getBoxes();
      int ncomp = a_flux.nComp();
      CH_assert(a_oldPhi.nComp() == ncomp);
      CH_assert(a_lowOrderFlux.nComp() == ncomp);

      // actually compute limited values


      LevelData<FluxBox> antiDiffusiveFlux(grids,
                                           ncomp,
                                           a_flux.ghostVect());



      // don't think I need ghost cells for this guy...
      LevelData<FArrayBox> lowOrderUpdate(grids,
                                          ncomp);

      // implementation is simpler if we include a ghost cell on this guy
      LevelData<FArrayBox> limiterVal(grids, ncomp, IntVect::Unit);

      // compute low-order update
      Real fakeDx = 1.0;
      SingleLevelDivergence::levelDivergenceMAC(lowOrderUpdate,
                                                a_lowOrderFlux,
                                                fakeDx);

      Real hDinv = 1.0/(pow(a_dx, SpaceDim));
      DataIterator dit = a_oldPhi.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box gridBox = grids[dit];

          FArrayBox& thisLowOrderPhi = lowOrderUpdate[dit];
          const FArrayBox& thisOldPhi = a_oldPhi[dit];

          thisLowOrderPhi *= -a_dt;
          thisLowOrderPhi *= hDinv;
          thisLowOrderPhi += thisOldPhi;

          // also compute antidiffusive flux while we're here
          const FluxBox& thisLowOrderFlux = a_lowOrderFlux[dit];
          FluxBox& thisHighOrderFlux = a_flux[dit];

          FluxBox& thisAntiDiffusiveFlux = antiDiffusiveFlux[dit];
          thisAntiDiffusiveFlux.copy(thisHighOrderFlux);
          thisAntiDiffusiveFlux -= thisLowOrderFlux;

          // now compute sum of fluxes out of each cell
          FArrayBox sumFluxesOut(grids[dit], ncomp);
          sumFluxesOut.setVal(0.0);

          // will eventually want to move this to fortran
          for (int dir=0; dir<SpaceDim; dir++)
            {
              IntVect offset = BASISV(dir);
              FArrayBox& antiDiffusiveFluxDir = thisAntiDiffusiveFlux[dir];
              BoxIterator bit(gridBox);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  for (int comp=0; comp<ncomp; comp++)
                    {
                      // lowSide
                      if (antiDiffusiveFluxDir(iv,comp) < 0.0)
                        {
                          sumFluxesOut(iv,comp) -= antiDiffusiveFluxDir(iv,comp);
                        }
                      // hiSide
                      if (antiDiffusiveFluxDir(iv+offset,comp) > 0.0)
                        {
                          sumFluxesOut(iv,comp) += antiDiffusiveFluxDir(iv+offset,comp);
                        }
                    } // end loop over components
                } // end loop over cells
            } // end loop over directions

          // sumFluxesOut now contains the sum of outward-directed
          // anti-diffusivefluxes for each cell

          // now compute limiter values

          // set default value to one so that any ghost cells outside
          // valid domain have no effect
          FArrayBox& thisLimiterVal = limiterVal[dit];
          thisLimiterVal.setVal(1.0);

          Real cellVol = pow(a_dx, SpaceDim);
          //Real cellVol = a_dx;
          //#if 0
          BoxIterator bit(gridBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              for (int comp=0; comp<ncomp; comp++)
                {
                  // limiter value a la Zalesak
                  Real diff = (thisLowOrderPhi(iv,comp) - a_minVal)*cellVol;
                  if (sumFluxesOut(iv,comp) > 0)
                    {
                      thisLimiterVal(iv,comp) =min(1.0,
                                                   diff/(a_dt*sumFluxesOut(iv,comp)));
                    }
                  else
                    {
                      thisLimiterVal(iv,comp) = 0.0;
                    }
                }
            }
          //#endif
        }

      // fill in ghost cells appropriately
      limiterVal.exchange();

      // now limit fluxes
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = grids[dit];
          FArrayBox& thisLimiterVal = limiterVal[dit];
          FluxBox& thisAntiDiffusiveFlux = antiDiffusiveFlux[dit];

          for (int dir=0; dir<SpaceDim; dir++)
            {
              IntVect offset = BASISV(dir);
              FArrayBox& limitedFlux = thisAntiDiffusiveFlux[dir];
              FArrayBox& flux = a_flux[dit][dir];
              const FArrayBox& lowOrderFlux = a_lowOrderFlux[dit][dir];

              Box faceBox(gridBox);
              faceBox.surroundingNodes(dir);
              BoxIterator bit(faceBox);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  for (int comp=0; comp<ncomp; comp++)
                    {
                      // for each face, compute limiter value
                      IntVect iv = bit();
                      Real C = 0;
                      // rely on sign of antidiffusive flux to determine
                      // which limit value to use.
                      Real Aloc = limitedFlux(iv,comp);
                      // If we want to enforce max as well as min, this
                      // will become a min of the two limiter values
                      // A comes from right, use high-side limiter val
                      if (Aloc < 0)
                        {
                          C = thisLimiterVal(iv, comp);
                        }
                      else
                        {
                          // use low-side limiter val
                          C=thisLimiterVal(iv-offset, comp);
                        }

                      limitedFlux(iv,comp) = C*limitedFlux(iv,comp);

                      // finally, compute final limited flux value
                      flux(iv,comp) = lowOrderFlux(iv,comp) + limitedFlux(iv,comp);
                    }
                }
            } // end loop over flux directions
        } // end loop over grid boxes
    } // end if we're doing limiting
}

#include "NamespaceFooter.H"
