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
//  PhysMappedIBC.H
// ============
//
// Virtual base class through which a user specifies the initial and boundary
// conditions for a hyperbolic system of PDEs.
//

#ifndef _PHYSMAPPEDIBC_H_
#define _PHYSMAPPEDIBC_H_

#include "PhysIBC.H"
#include "MultiBlockCoordSys.H"
#include "CHArray_fwd.H"

#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**
   Virtual base class through which a user specifies the initial and
   boundary conditions for a hyperbolic system of PDEs.
 */
class PhysMappedIBC : public PhysIBC
{
public:
  /// Constructor
  /**
   */
  PhysMappedIBC();

  /// Destructor
  /**
   */
  virtual ~PhysMappedIBC();

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// Factory method - this object is its own factory
  /**
     Return a point to a new PhysMappedIBC object with m_isDefined = false
     (i.e., its define() must be called before it is used).
   */
  virtual PhysMappedIBC* new_physIBC() = 0;

  virtual void setTime(Real a_time);

  //! Set the coordinate system for this IBC. This object does not assume
  //! of the coordinate system.
  virtual void setCoordSys(MultiBlockCoordSys* a_coordSysPtr);

  //! Set the coordinate system for this IBC. This object does not assume
  //! of the coordinate system. This method is the single-block version
  //! of setCoordSys.
  void setCoordSys(NewCoordSys* a_coordSysPtr);

  virtual void print() {}

  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U) = 0;

  virtual void initializeWithJ(LevelData<FArrayBox>& a_U)
    {
      MayDay::Error("Initialize with J not defined for this derivation from "
                    "PhysMappedIBC");
    }

  /// Old set boundary fluxes for Cartesian methods -- will cause an error
  /**
   */
  virtual void primBC(FArrayBox&            a_WGdnv,
                      const FArrayBox&      a_Wextrap,
                      const FArrayBox&      a_W,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time)
    {
      MayDay::Error("Cartesian interface to primBC not valid for mapped grids");
    }

  /// Set boundary fluxes
  /**
   */
  virtual void primBC(FArrayBox&             a_WGdnv,
                      const FArrayBox&       a_Wextrap,
                      const FArrayBox&       a_W,
                      const FArrayBox *const a_unitNormalBasisPtr,
                      const Interval&        a_velIntv,
                      const int&             a_dir,
                      const Side::LoHiSide&  a_side,
                      const Real&            a_time) = 0;

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
                     const Real&      a_time) = 0;

  /// Old artificial viscosity for Cartesian methods -- will cause an error
  /**
   */
  virtual
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time)
    {
      MayDay::Error("Cartesian interface to artViscBC not valid for mapped "
                    "grids");
    }

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  virtual
  void artViscBC(FArrayBox&                                   a_NtFdir,
                 const CHArray<Real, SpaceDim+1, ArRangeCol>& a_Nctg,
                 const FArrayBox&                             a_U,
                 const FArrayBox&                             a_unitNormalBasis,
                 const FArrayBox&                             a_divVel,
                 const FArrayBox&                             a_csq,
                 const FArrayBox&                             a_dxFace,
                 const Interval&                              a_momIntv,
                 const Real                                   a_alpha,
                 const Real                                   a_beta,
                 const Box&                                   a_loFaceBox,
                 const int                                    a_hasLo,
                 const Box&                                   a_hiFaceBox,
                 const int                                    a_hasHi,
                 const int                                    a_dir) = 0;

  /// Do we have the exact solution
  bool haveExactSoln() const;

protected:

  MultiBlockCoordSys* m_coordSysPtr;

  Real m_time;

  bool m_haveCoordSys;

  bool m_haveTime;

  bool m_haveExactSoln;
};

#include "NamespaceFooter.H"
#endif
