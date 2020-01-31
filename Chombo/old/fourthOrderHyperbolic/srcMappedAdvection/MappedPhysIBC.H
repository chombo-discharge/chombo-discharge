#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  MappedPhysIBC.H
// ============
//
// Interface between AMRTimeDependent PhysIBC class needed for
// CTU advection and the BasicIBC class used in the mapped advection
// example. The basic idea is that this _is_ a PhysIBC which _contains_
// a BasicIBC pointer which does the heavy lifting for the BC.  This is
// a simpler and cleaner solution than trying to derive BasicIBC from PhysIBC
// directly.
//

#ifndef _MAPPEDPHYSIBC_H_
#define _MAPPEDPHYSIBC_H_

#include "PhysIBC.H"
#include "BasicIBC.H"

#include "UsingNamespace.H"

class BasicIBC;

class
MappedPhysIBC : public PhysIBC
{
  friend class BasicIBC;

public:
  /// Constructor
  /**
   */
  MappedPhysIBC();

  /// Destructor
  /**
   */
  virtual ~MappedPhysIBC();

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object. Note that in this context,
     a_dx is really dXi
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// Factory method - this object is its own factory
  /**
     Return a point to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  virtual PhysIBC* new_physIBC();

  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U);

  /// Set boundary fluxes
  /**
   */
  virtual void primBC(FArrayBox&            a_WGdnv,
                      const FArrayBox&      a_Wextrap,
                      const FArrayBox&      a_W,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time);

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  virtual
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  virtual
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time);

  /// set mapped-grid basicIBC member
  void setBasicIBC(BasicIBC* a_basicIBCptr);

  /// set coordinate-system member
  /**
      note that this pointer is not managed by the MappedPhysIBC object.
      (in other words, it's not deleted when this object is destructed)
  */
  void setCoordSys(CoordSys<FArrayBox,FluxBox>* a_coordSysPtr);

protected:

  /// BasicIBC class which does all the mapped-grid IBC heavy lifting
  BasicIBC* m_basicIBCptr;

  /// mapped coordinate system --
  /** note that this is not managed by this object -- it is not freed upon
      deletion of this MappedPhysIBC object.
  */
  CoordSys<FArrayBox,FluxBox>* m_coordSysPtr;

};

#endif
