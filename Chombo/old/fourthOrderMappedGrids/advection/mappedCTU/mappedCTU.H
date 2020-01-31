#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _MAPPEDCTU_H_
#define _MAPPEDCTU_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "FourthOrderCoordSys.H"
#include "PatchGodunov.H"

#include "NamespaceHeader.H"

/// function to compute mapped-grid CTU advective fluxes
/** note that this returns N^T vel phi at half-times on faces,
    if the flux = vel*phi...
*/
void
computeCTUFlux(LevelData<FluxBox>& a_CTUflux,
               const LevelData<FArrayBox>& a_oldPhi,
               const LevelData<FluxBox>& a_faceVel,
               FourthOrderCoordSys* a_coordSysPtr,
               PatchGodunov& a_CTUpatchGod,
               Real a_time,
               Real a_dt);

#include "NamespaceFooter.H"

#endif
