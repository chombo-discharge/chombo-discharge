#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELMAPPEDCONSFACTORY_H_
#define _AMRLEVELMAPPEDCONSFACTORY_H_

#include "AMRLevelConsFactory.H"
#include "AMRLevelMappedCons.H"

#include "NamespaceHeader.H"

/// AMR Cons Equation factory
/**
 */
class AMRLevelMappedConsFactory : public AMRLevelConsFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelMappedConsFactory();

  /// Constructor for a single-block coordinate system.
  /**
   */
  AMRLevelMappedConsFactory(
      NewCoordSysFactory *const a_coordSysFact,
      AMRLevelMappedStabilityStrategy* a_stableTimeSteppingStrategy,
      AMRLevelMappedTaggingStrategy* a_taggingStategy);

  /// Constructor for a multi-block coordinate system.
  /**
   */
  AMRLevelMappedConsFactory(
      MultiBlockCoordSysFactory *const a_coordSysFact,
      AMRLevelMappedStabilityStrategy* a_stableTimeSteppingStrategy,
      AMRLevelMappedTaggingStrategy* a_taggingStategy);

  /// Virtual constructor
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  /**
   */
  virtual ~AMRLevelMappedConsFactory();

  /// Set the factory for the coordinate system
  /**
   */
  void coordinateSystemFactory(MultiBlockCoordSysFactory* a_coordSysFact);

  /// Set whether to use source term.
  /**
   */
  virtual void useSourceTerm(bool a_useSourceTerm);

  /// Source term class, if necessary
  /**
   */
  virtual void sourceTerm(const LevelSourceTerm* const a_sourceTerm);

  /// Set whether to modify fluxes by adding artificial dissipation based on 5th derivative.
  /**
   */
  virtual void useArtificialDissipation(bool a_useArtificialDissipation);

  /// Set coefficient of artificial dissipation.
  /**
   */
  virtual void artificialDissipation(Real a_artificialDissipation);

  /**
   */
  void plotPrefix(const std::string& a_plotfile_prefix);

  /// Check that everything is defined
  /**
   */
  bool isDefined() const;

  virtual void writeJ(bool a_writeJ);

  virtual void writeExact(bool a_writeExact);

  virtual void writeError(bool a_writeError);

  virtual void writeMap(bool a_writeMap);

protected:

  // Some default values
  virtual void setDefaultValues();

  virtual void transferSettings(AMRLevelCons* a_newPtr) const;

  // Single block?  Default false.  If true, then m_coordSysFact will be newed and deleted.
  bool m_singleBlock;

  // Use source term?
  bool m_useSourceTerm;

  // Pointer to the object that supplies source term.
  const LevelSourceTerm* m_sourceTerm;

  // Use artificial dissipation?
  bool m_useArtificialDissipation;

  // Coefficient of artificial dissipation
  Real m_artificialDissipation;

  // Factory for coordinate system
  MultiBlockCoordSysFactory* m_coordSysFact;

  // Prefix for plot files
  std::string m_plotPrefix;

  // Whether to write out J.
  bool m_writeJ;

  // Whether to write out exact solution.
  bool m_writeExact;

  // Whether to write out error (difference between calculated and exact solutions).
  bool m_writeError;

  /// Whether to write out mapping file
  bool m_writeMap;

  // Strategy for computing a stable time step.
  AMRLevelMappedStabilityStrategy* m_stabilityStrategy;

  // Strategy for determining which cells to tag for refinement.
  AMRLevelMappedTaggingStrategy* m_taggingStrategy;

private:

  // Copy and assignment not allowed
  void operator=(const AMRLevelMappedConsFactory& a_input);
  AMRLevelMappedConsFactory(const AMRLevelMappedConsFactory& a_input);
};

#include "NamespaceFooter.H"

#endif
