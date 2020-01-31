#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _ZALESAKLIMITER_H_
#define _ZALESAKLIMITER_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "FourthOrderCoordSys.H"
#include "mappedLimiter.H"

#include "NamespaceHeader.H"

  /// apply Zalesak limiter to enforce minVal
  /** a_flux comes in with high-order flux and is modified in place using
      the low-order flux.
      -- if lowOrderFluxesOnly is true, then simply copy lowOrderFlux into
         flux (useful for debugging, at times)
  */
void
applyZalesakLimiter(LevelData<FluxBox>& a_flux,
                    const LevelData<FluxBox>& a_lowOrderFlux,
                    const LevelData<FArrayBox>& a_oldPhi,
                    Real a_minVal,
                    Real a_dx,
                    Real a_dt,
                    bool a_lowOrderFluxesOnly);

#include "NamespaceFooter.H"

#endif
