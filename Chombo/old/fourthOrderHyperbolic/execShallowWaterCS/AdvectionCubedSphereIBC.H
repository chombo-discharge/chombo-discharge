#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTIONCUBEDSPHEREIBC_H_
#define _ADVECTIONCUBEDSPHEREIBC_H_

#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "ScalarFunction.H"
#include "VectorFunction.H"
#include "RefCountedPtr.H"
#include "PhysAdvectMappedIBC.H"

#include "NamespaceHeader.H"

class AdvectionCubedSphereIBC: public PhysAdvectMappedIBC
{
  public:


  //! Constructor with initial solution, advection velocity.
  /**
   */
  AdvectionCubedSphereIBC(RefCountedPtr<ScalarFunction> a_solution,
                          RefCountedPtr<VectorFunction> a_advectionVelocity);

  //! Destructor
  ~AdvectionCubedSphereIBC();

  // Overridden methods.

  PhysMappedIBC *new_physIBC(); // Factory method.

  void initialize(LevelData<FArrayBox>& a_U);

  void initializeWithJ(LevelData<FArrayBox>& a_U);

  void primBC(FArrayBox&            a_WGdnv,
              const FArrayBox&      a_Wextrap,
              const FArrayBox&      a_W,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time);

  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time);

  //! Return advection velocity.
  virtual void advVel(LevelData<FArrayBox>& a_advVel,
                      Real a_time = 0.);

  protected:

  void initializeUnified(LevelData<FArrayBox>& a_U,
                         bool a_includeJ);

  // Functions.
  RefCountedPtr<ScalarFunction> m_solution;
  RefCountedPtr<VectorFunction> m_velocity;

  bool m_gotAdvVel;
  LevelData<FArrayBox> m_advVel;

  private:

  // Forbidden.
  AdvectionCubedSphereIBC();
  AdvectionCubedSphereIBC(const AdvectionCubedSphereIBC&);
  AdvectionCubedSphereIBC& operator=(const AdvectionCubedSphereIBC&);
};

#include "NamespaceFooter.H"

#endif
