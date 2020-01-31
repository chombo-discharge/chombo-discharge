#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELSHALLOWWATERMAPPEDOPERATOR_H_
#define _LEVELSHALLOWWATERMAPPEDOPERATOR_H_

#include "LevelMappedConsOperator.H"
#include "PatchShallowWaterMappedOperator.H"
#include "PhysShallowWaterMappedIBC.H"

#include "NamespaceHeader.H"

/// Level Cons Operator - a class for evaluating the Laplacian.

/**
 */
class LevelShallowWaterMappedOperator : public LevelMappedConsOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  LevelShallowWaterMappedOperator();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~LevelShallowWaterMappedOperator();

  virtual void define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                      const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                      const ProblemDomain&      a_domain,
                      const int&                a_refineCoarse,
                      const Real&               a_dx,
                      const MOLPhysics* const  a_molPhysics,
                      const int&                a_numStates,
                      const bool&               a_hasCoarser,
                      const bool&               a_hasFiner);

  virtual void defineMetricStuff();

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

  virtual void evalRHSpatches(
                              LevelData<FArrayBox>&       a_LofU,
                              const LevelData<FArrayBox>& a_U,
                              LevelFluxRegister&          a_finerFluxRegister,
                              LevelFluxRegister&          a_coarserFluxRegister,
                              Real                        a_weight);

  virtual void exchangeGhosts(LevelData<FArrayBox>&   a_U);

protected:

  virtual void setPatchIndex(const DataIndex&  a_ind) const;

  virtual void setCommonFlux(LevelData<FluxBox>&   a_flux,
                             const BlockRegister&  a_blockRegister) const;

  // Patch integrator
  // PatchShallowWaterMappedOperator m_patchConsOperator;

  PhysShallowWaterMappedIBC* m_physIBCPtr;

  LevelData<FluxBox> m_contravariantMetricFace;

  LevelData<FluxBox> m_orthoMatrix;
};

#include "NamespaceFooter.H"

#endif
