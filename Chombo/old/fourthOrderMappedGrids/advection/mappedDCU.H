#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _MAPPEDDCU_H_
#define _MAPPEDDCU_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "FourthOrderCoordSys.H"

#include "NamespaceHeader.H"

/// function to compute mapped-grid DCU advective fluxes
/** note that this returns N^T vel phi at half-times on faces,
    if the flux = vel*phi...
*/
void
computeDCUFlux(LevelData<FluxBox>& a_DCUflux,
               const LevelData<FArrayBox>& a_oldPhi,
               const LevelData<FluxBox>& a_faceVel,
               FourthOrderCoordSys* a_coordSysPtr,
               Real a_time,
               Real a_dt);

#include "NamespaceFooter.H"

#endif
