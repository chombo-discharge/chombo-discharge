#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELCLAWFACTORY_H_
#define _AMRLEVELCLAWFACTORY_H_

#include "AMRLevelFactory.H"
#include "AMRLevelClaw.H"
#include "ClawPatch.H"

#include "UsingNamespace.H"

/// AMR ClawPack factory
/**
 */
class AMRLevelClawFactory : public AMRLevelFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelClawFactory();

  /// Virtual constructor
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  /**
   */
  virtual ~AMRLevelClawFactory();

  /// CFL number
  /**
   */
  virtual void CFL(Real a_cfl);

  /// Physical dimension of the longest side of the domain
  /**
   */
  virtual void domainLength(Real a_domainLength);

//   void IBC(ClawIBC* a_w0)
//   {
//     m_wave_ibc = a_w0->new_physIBC();
//   };
//

  ///
  /**
   */
  virtual void verbosity(const int& verbosity);

  /// Refinement threshold
  /**
   */
  virtual void refinementThreshold(Real a_refineThresh);

  /// Clawpack parameters
  /**
   */

  void set_stateNames(vector<string> stateNames);

  /// Tag buffer size
  /**
   */
  void tagBufferSize(int a_tagBufferSize);

  /// Initial dt multiplier
  /**
   */
  void initialDtMultiplier(Real a_initialDtMultiplier);

  void clawPatch(const ClawPatch& a_clawPatch);

  /// Check that everything is defined
  /**
   */
  bool isDefined() const;



protected:
  // Some default values
  void setDefaultValues();

  int m_verbosity;

  // CFL number
  Real m_cfl;
  bool m_cflSet;
  // ClawIBC* m_wave_ibc;  // Now in ClawPatch

  // Physical dimension of the longest side of the domain
  Real m_domainLength;
  bool m_domainLengthSet;

  // Refinement threshold for gradient
  Real m_refineThresh;
  bool m_refineThreshSet;

  // Tag buffer size
  int  m_tagBufferSize;
  bool m_tagBufferSizeSet;

  // Initial dt multiplier
  Real m_initialDtMultiplier;
  bool m_initialDtMultiplierSet;

  std::vector<string> m_stateNames;

    // Patch integrator (used as a factory)
    ClawPatch m_clawPatch; // Used topy everything to AMRLevelClaw's
    bool                     m_clawPatchSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelClawFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelClawFactory(const AMRLevelClawFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
