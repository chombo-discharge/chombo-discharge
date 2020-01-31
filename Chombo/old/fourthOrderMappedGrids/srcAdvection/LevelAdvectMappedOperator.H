#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELADVECTMAPPEDOPERATOR_H_
#define _LEVELADVECTMAPPEDOPERATOR_H_

#include "LevelMappedConsOperator.H"
#include "PatchAdvectMappedOperator.H"
#include "PhysAdvectMappedIBC.H"

#include "NamespaceHeader.H"

/// Level Cons Operator - a class for evaluating the Laplacian.

/**
 */
class LevelAdvectMappedOperator : public LevelMappedConsOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  LevelAdvectMappedOperator();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~LevelAdvectMappedOperator();

  virtual void define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                      const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                      const ProblemDomain&      a_domain,
                      const int&                a_refineCoarse,
                      const Real&               a_dx,
                      const MOLPhysics* const  a_molPhysics,
                      const int&                a_numStates,
                      const bool&               a_hasCoarser,
                      const bool&               a_hasFiner);

  virtual void defineAdvVel();

  virtual void evalRHS(
                       LevelData<FArrayBox>&       a_LofU,
                       LevelData<FArrayBox>&       a_U,
                       LevelFluxRegister&          a_finerFluxRegister,
                       LevelFluxRegister&          a_coarserFluxRegister,
                       const LevelData<FArrayBox>& a_UcoarseOld,
                       const Real&                 a_timeCoarseOld,
                       const LevelData<FArrayBox>& a_UcoarseNew,
                       const Real&                 a_timeCoarseNew,
                       Real                        a_time,
                       Real                        a_weight);

  virtual void setAdvVel(Real a_time);

  LevelData<FluxBox>& advVelFace()
  {
    return m_advVelFace;
  };

protected:

  virtual void setPatchIndex(const DataIndex&  a_ind) const;

  // Patch integrator
  // PatchAdvectMappedOperator m_patchConsOperator;

  PhysAdvectMappedIBC* m_physIBCPtr;

  // CELL-centered
  LevelData<FArrayBox> m_advVel;

  // FACE-centered
  LevelData<FluxBox> m_advVelFace;
};

#include "NamespaceFooter.H"

#endif
