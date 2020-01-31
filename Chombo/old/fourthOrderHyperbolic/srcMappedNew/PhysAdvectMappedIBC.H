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
//  PhysAdvectMappedIBC.H
// ============
//
// Virtual base class through which a user specifies the initial and boundary
// conditions for advection.
//

#ifndef _PHYSADVECTMAPPEDIBC_H_
#define _PHYSADVECTMAPPEDIBC_H_

#include "PhysMappedIBC.H"

#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**
   Virtual base class through which a user specifies the initial and
   boundary conditions for a hyperbolic system of PDEs.
 */
class PhysAdvectMappedIBC : public PhysMappedIBC
{
public:
  /// Constructor
  /**
   */
  PhysAdvectMappedIBC();

  /// Destructor
  /**
   */
  virtual ~PhysAdvectMappedIBC();

  ///
  /* Return advection velocity.
   */
  virtual void advVel(LevelData<FArrayBox>& a_advVel) = 0;
};

#include "NamespaceFooter.H"
#endif
