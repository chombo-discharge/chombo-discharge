#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COMPUTEMAPPEDDT_H_
#define _COMPUTEMAPPEDDT_H_

#include "FourthOrderCoordSys.H"
#include "LevelData.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"
// helper function to compute max save dt for advection
Real computeMappedDt(const LevelData<FluxBox>& a_faceVel,
                     const FourthOrderCoordSys* a_coordSysPtr,
                     Real a_cfl,
                     const bool velocities_are_physical=false);


#include "NamespaceFooter.H"

#endif
