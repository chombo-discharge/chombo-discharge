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
  void coordinateSystemFactory(NewCoordSysFactory* a_coordSysFact);

  /**
   */
  void plotPrefix(const std::string& a_plotfile_prefix);

  /// Check that everything is defined
  /**
   */
  bool isDefined() const;

protected:

  // Some default values
  void setDefaultValues();

  virtual void transferSettings(AMRLevelMappedCons* a_newPtr) const;

  // whether to compute timestep from cell data (the old way)
  bool m_dtFromCells;

  // Pointer to the object that supplies the physics needed by the integrator
  // (used as a factory here)
  const MOLPhysics* m_gdnvPhysics;

  // Coordinate system factory class
  NewCoordSysFactory* m_coordSysFactPtr;
  bool m_coordSysFactSet;

  std::string m_plotfile_prefix;
};

#include "NamespaceFooter.H"

#endif
