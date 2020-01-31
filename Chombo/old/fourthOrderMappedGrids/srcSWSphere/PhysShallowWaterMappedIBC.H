#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//
//  PhysShallowWaterMappedIBC.H
// ============
//
// Virtual base class through which a user specifies the initial and boundary
// conditions for advection.
//

#ifndef _PHYSSHALLOWWATERMAPPEDIBC_H_
#define _PHYSSHALLOWWATERMAPPEDIBC_H_

#include "PhysMappedIBC.H"

#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**
   Virtual base class through which a user specifies the initial and
   boundary conditions for a hyperbolic system of PDEs.
 */
class PhysShallowWaterMappedIBC : public PhysMappedIBC
{
public:
  /// Constructor
  /**
   */
  PhysShallowWaterMappedIBC();

  /// Destructor
  /**
   */
  virtual ~PhysShallowWaterMappedIBC();

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  virtual void setFortranCommon(const Real& a_gravity,
                                const Real& a_omega,
                                const Real& a_alpha);

  ///
  /* Return topography height.
   */
  virtual void topoHeight(LevelData<FArrayBox>& a_topoHeight);

protected:
  // True if the Fortran common block has been set
  bool m_isFortranCommonSet;

  // Constant of gravity
  Real m_gravity;

  // Rotation rate of the Earth
  Real m_omega;

  // Rotation angle of the grid
  Real m_alpha;
};

#include "NamespaceFooter.H"
#endif
