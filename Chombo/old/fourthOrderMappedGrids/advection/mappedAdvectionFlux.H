#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _MAPPEDADVECTIONFLUX_H_
#define _MAPPEDADVECTIONFLUX_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "FourthOrderCoordSys.H"
#include "mappedLimiter.H"

#include "NamespaceHeader.H"

/// compute high-order flux for advection equation on mapped grids.
/** Assumes all physical boundary conditions have been set.
    -- also returns face-averaged face advection velocity in faceVel
 */
void
computeMappedFourthOrderFlux(LevelData<FluxBox>& a_compCoordFlux,
                             LevelData<FluxBox>& a_faceVel,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>* a_advVelPtr,
                             const FourthOrderCoordSys* a_FOCS,
                             const RealVect& a_dx,
                             bool a_limitFaceValues
                             );


void
cellUJToCellU(
              LevelData<FArrayBox>& a_cell_u,
              const LevelData<FArrayBox>& a_cell_uJ,
              const FourthOrderCoordSys* a_FOCS);


void
computeCompFaceFluxes( LevelData<FluxBox>& a_uTimesV,
                       const LevelData<FluxBox>& a_u,
                       const LevelData<FluxBox>& a_v,
                       bool a_useFourthOrder);

#include "NamespaceFooter.H"

#endif
