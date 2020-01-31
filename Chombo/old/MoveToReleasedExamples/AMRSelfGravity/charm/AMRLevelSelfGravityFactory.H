#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELSELFGRAVITYFACTORY_H_
#define _AMRLEVELSELFGRAVITYFACTORY_H_

#include "AMRLevelFactory.H"
#include "AMRLevelSelfGravity.H"

/// AMRLevelFactory-derived Factory object for AMRLevelSelfGravity
class AMRLevelSelfGravityFactory : public AMRLevelFactory
{
public:
  /// Null constructor
  AMRLevelSelfGravityFactory();

  /// Virtual constructor
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  virtual ~AMRLevelSelfGravityFactory();

  ///
  virtual void define(const Real&                 a_cfl,
                     const Real&                 a_domainLength,
                     const int&                  a_verbosity,
                     const int&                  a_tagBufferSize,
                     const int&                  a_maxInitRefLevel,
                     const Real&                 a_initialDtMultiplier,
                     const GodunovPhysics* const a_godunovPhysics,
                     const int&                  a_normalPredOrder,
                     const bool&                 a_useFourthOrderSlopes,
                     const bool&                 a_usePrimLimiting,
                     const bool&                 a_useCharLimiting,
                     const bool&                 a_useFlattening,
                     const bool&                 a_useArtificialViscosity,
                     const Real&                 a_artificialViscosity,
                     const RefCellTagger* const  a_refCellTagger,
                     const bool&                 a_useDeltaPhiCorr,
                      const StencilType&          a_stencil);

  /// Check that everything is defined
  bool isDefined() const;

protected:
  // Some default values
  void setDefaultValues();

  //
  int m_verbosity;

  // CFL numbers for hydro
  Real m_cfl;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // max level of refinement at start
  int m_maxInitRefLevel;

  // Tag buffer size
  int  m_tagBufferSize;

  // Initial dt multiplier
  Real m_initialDtMultiplier;

  // Order of the normal predictor (1 -> PLM, 2 -> PPM)
  int m_normalPredOrder;

  // Use 4th other slope computations (otherwise, use 2nd order)
  bool m_useFourthOrderSlopes;

  // Do slope limiting in the primitive or characteristic variables, respect.
  bool m_usePrimLimiting;
  bool m_useCharLimiting;

  // Do slope flattening - MUST BE USING 4th order slopes
  bool m_useFlattening;

  // Apply artificial viscosity of a set value
  bool m_useArtificialViscosity;
  Real m_artificialViscosity;

  // stencil for the derivative (eg to get F from phi)
  StencilType m_stencil;

  // whether or not the deltaPhi correction should be applied
  bool m_useDeltaPhiCorr;

  // Pointer to the object that supplies the physics needed by the integrator
  // (used as a factory)
  const GodunovPhysics* m_godunovPhysics;

  //
  const RefCellTagger* m_refCellTagger;

  // is all defined
  bool m_isDefined;

private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelSelfGravityFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelSelfGravityFactory(const AMRLevelSelfGravityFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
