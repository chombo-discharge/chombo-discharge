#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COMPUTEMAPPEDNEWDT_H_
#define _COMPUTEMAPPEDNEWDT_H_

#include "MultiBlockCoordSys.H"
#include "LevelData.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"

// helper function to compute max save dt for advection
Real computeMappedNewDt(const LevelData<FluxBox>& a_faceVel,
                        const MultiBlockCoordSys* a_coordSysPtr,
                        Real a_cfl);


#include "NamespaceFooter.H"

#endif
