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

  virtual PhysMappedIBC *new_physIBC() = 0;

  virtual void initialize(LevelData<FArrayBox>& a_U) = 0;

  ///
  /* Return advection velocity.
   */
  virtual void advVel(FArrayBox& a_advVel,
                      Real a_time = 0.) = 0;

  ///
  /* These are now pure virtual in PhysMappedIBC.
   */
  virtual void primBC(FArrayBox&             a_WGdnv,
                     const FArrayBox&       a_Wextrap,
                     const FArrayBox&       a_W,
                     const FArrayBox *const a_unitNormalBasisPtr,
                     const Interval&        a_velIntv,
                     const int&             a_dir,
                     const Side::LoHiSide&  a_side,
                     const Real&            a_time)
  { };

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
                 const int                                    a_dir)
  { };

  virtual void setBdrySlopes(FArrayBox&       a_dW,
                             const FArrayBox& a_W,
                             const int&       a_dir,
                             const Real&      a_time);
};

#include "NamespaceFooter.H"
#endif
