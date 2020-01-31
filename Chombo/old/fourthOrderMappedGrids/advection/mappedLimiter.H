#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _MAPPEDLIMITER_H_
#define _MAPPEDLIMITER_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "FourthOrderCoordSys.H"

#include "NamespaceHeader.H"


void
applyMappedLimiter(LevelData<FluxBox>& a_facePhi,
                   const LevelData<FArrayBox>& a_cellPhi,
                   const LevelData<FluxBox>& a_faceVel,
                   const FourthOrderCoordSys* a_FOCS,
                   const RealVect& a_dx,
                   const bool velocites_are_physical=false);

#include "NamespaceFooter.H"

#endif
