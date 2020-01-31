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

  /// Virtual constructor
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  /**
   */
  virtual ~AMRLevelMappedConsFactory();

  /// Whether to compute timestep from cell-centered velocities (old way)
  /**
   */
  void dtFromCells(bool a_dtFromCells);

  /// Coordinate system factory
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

  /**
   */
  void plotPrefix(const std::string& a_plotfile_prefix);

  /// Whether to write mapped-grid geometry file
  /**
   */
  void writeMap(bool a_writeMap);

  /// Whether to include error in output
  /**
   */
  void writeError(bool a_writeError);

  /// Check that everything is defined
  /**
   */
  bool isDefined() const;

protected:

  // Some default values
  virtual void setDefaultValues();

  virtual void transferSettings(AMRLevelMappedCons* a_newPtr) const;

  // whether to compute timestep from cell data (the old way)
  bool m_dtFromCells;

  // Use source term?
  bool m_useSourceTerm;

  // Pointer to the object that supplies source term.
  const LevelSourceTerm* m_sourceTerm;

  // Pointer to the object that supplies the physics needed by the integrator
  // (used as a factory here)
  //  const MOLPhysics* m_molPhysics;

  // Coordinate system factory class
  MultiBlockCoordSysFactory* m_coordSysFactPtr;
  bool m_coordSysFactSet;

  std::string m_plotfile_prefix;
};

#include "NamespaceFooter.H"

#endif
