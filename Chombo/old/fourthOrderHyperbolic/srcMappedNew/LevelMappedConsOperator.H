#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELMAPPEDCONSOPERATOR_H_
#define _LEVELMAPPEDCONSOPERATOR_H_

// Original design:                     Phil Colella
// Original implementation:
// Implementation date:

#include "LevelConsOperator.H"
#include "MOLPhysics.H"
#include "NewFourthOrderCoordSys.H"

#include "NamespaceHeader.H"

/// Level Cons Operator - a class for evaluating the Laplacian.

/**
 */
class LevelMappedConsOperator : public LevelConsOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  LevelMappedConsOperator();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~LevelMappedConsOperator();

  /// set coordinate system object
  void setCoordSys(NewFourthOrderCoordSys* a_coordSysPtr);

  void cellUJToCellU(LevelData<FArrayBox>& a_Uavg,
                     const LevelData<FArrayBox>& a_UJavg) const;

protected:

  // put this in PatchMappedConsOperator
  // as well as computeFaceAreas computeFaceGradU computeNJInvGradUOnBox
  // computeNTNJInvGradUOnBox computeNormalDiffusiveFluxVector
  // secondOrderMappedGridDivergence divideOutCellVolume secondOrderLaplacian

  virtual void computeCompFaceFluxes(LevelData<FluxBox>& uTimesV,
                                     const LevelData<FluxBox>& u,
                                     const LevelData<FluxBox>& v) const;

  virtual void getPhysicalCellVolumes(LevelData<FArrayBox>&) const;

  virtual void computeFaceAverages(LevelData<FluxBox>& face_data,
                                   const LevelData<FArrayBox>& cell_data ) const;

  NewFourthOrderCoordSys* m_coordSysPtr;
};

#include "NamespaceFooter.H"

#endif
