#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// SingleLevelGradient.H

#ifndef _SINGLELEVELGRADIENT_H_
#define _SINGLELEVELGRADIENT_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"


/// Class to encapsulate Gradient functions (both CC and face-centered)
/** (DFM-10/17/08) -- note that this version has been gutted to only have
    the most basic single-level Gradient because AMRTools doesn't compile
    in 4D yet */
class SingleLevelGradient
{

public:

  /// computes cell-centered, level-operator gradient of cell-centered phi
  /** in this case, assume that all relevant
      BC's (coarse-fine and physical) have already been set, so
      phi can be a const variable */
  static void levelGradientCC(
                              ///
                              LevelData<FArrayBox>& a_grad,
                              ///
                              const LevelData<FArrayBox>& a_phi,
                              ///
                              const Real a_dx);


  /// computes cell-centered composite gradient of cell-centered phi
  /**
      this one assumes that all ghost-cell values have already been set;
      if phiFinePtr != NULL, then also uses one-sided differencing to compute
      gradient on coarse side of corarse-fine interface.  note that gradient
      is only really defined in valid regions of grids. */
  static void compGradientCC(
                             ///
                             LevelData<FArrayBox>& a_Grad,
                             ///
                             const LevelData<FArrayBox>& a_phi,
                             ///
                             const LevelData<FArrayBox>* a_phiFinePtr,
                             ///
                             const Real a_dx,
                             ///
                             const int a_nRefFine,
                             ///
                             const ProblemDomain& a_dProblem);

  /// computes edge-centered level-operator gradient of cell-centered phi
  /**
      assumes _ALL_ ghost cell values have been preset (so phi can be const)
  */
  static void levelGradientMAC(
                               ///
                               LevelData<FluxBox>& a_edgeGrad,
                               ///
                               const LevelData<FArrayBox>& a_phi,
                               ///
                               const Real a_dx);


#ifdef NEW_GRADIENT
  /// utility function for internal use
  static void singleBoxMacGrad(FArrayBox& a_gradFab,
                               const FArrayBox& a_phiFab,
                               int a_gradComp,
                               int a_phiComp,
                               int a_numComp,
                               const Box& a_edgeBox,
                               Real a_dx,
                               int a_dir,
                               int a_edgeDir,
                               const IntVectSet& a_gridIVS);



  /// create grid IVS for all grids
  /** creates an IVS representation of the DisjointBoxLayout
      in the neighborhood around each box in the DBL.
  */
  static void createGridIVS(LayoutData<IntVectSet>& a_gridsIVS,
                            const DisjointBoxLayout& a_grids,
                            const int a_nGrow);


  /// internal utility function
  /** creates IVS representation of the DisjointBoxLayout
      in the neighborhood of localBox
  */
  static void createGridIVS(IntVectSet& a_gridIVS,
                            const DisjointBoxLayout& a_grids,
                            const Box& a_localBox,
                            const int a_nGrow = 1);

#endif

};

#include "NamespaceFooter.H"

#endif


