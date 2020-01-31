#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _SLIPWALLVISCOUSTENSORDOMAINBC_H_
#define _SLIPWALLVISCOUSTENSORDOMAINBC_H_

#include "RefCountedPtr.H"
#include "ViscousBaseDomainBC.H"
#include "BaseBCFuncEval.H"
#include "BCFunc.H"
#include "NamespaceHeader.H"

///
/**
 */
class SlipWallViscousTensorDomainBC: public ViscousBaseDomainBC
{
public:
  ///
  /**
   */
  SlipWallViscousTensorDomainBC();

  ///
  /**
   */
  virtual ~SlipWallViscousTensorDomainBC();


  ///
  /**
   */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  /**
     Elliptic solver flux.
   */
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);


  virtual void 
  fillVelGhost(FArrayBox&     a_state,
               const Box&     a_valid,
               const Box&     a_domain,
               Real           a_dx,
               bool           a_homogeneous)
  {
    Box valid = a_valid;
    for (int i=0; i<CH_SPACEDIM; ++i)
      {
        Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
        Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
        if (!a_domain.contains(ghostBoxLo))
          {
            ReflectiveVectorBC(a_state,
                               valid,
                               a_dx,
                               i,
                               Side::Lo, 1);
          }
        
        if (!a_domain.contains(ghostBoxHi))
          {
            ReflectiveVectorBC(a_state,
                               valid,
                               a_dx,
                               i,
                               Side::Hi, 1);
          }
        //grow so we hit corners
        valid.grow(i,1);
      }
  }

private:
  void getDivergenceAndGradient(Real&                 a_divergence,
                                Real*                 a_gradient,
                                const int&            a_idir,
                                const FaceIndex&      a_bndryFace,
                                const VolIndex&       a_vof,
                                const EBCellFAB&      a_phi,
                                const RealVect&       a_probLo,
                                const RealVect&       a_dx,
                                const Side::LoHiSide& a_side,
                                const DataIndex&      a_dit,
                                const Real&           a_time,
                                const bool&           a_useHomogeneous);

};

class SlipWallViscousTensorDomainBCFactory: public BaseDomainBCFactory
{
public:
  ///
  /**
   */
  SlipWallViscousTensorDomainBCFactory();

  ///
  /**
   */
  virtual ~SlipWallViscousTensorDomainBCFactory();

  ///
  /**
   */
  virtual SlipWallViscousTensorDomainBC* create(const ProblemDomain& a_domain,
                                                const EBISLayout&    a_layout,
                                                const RealVect&      a_dx);

private:

};

#include "NamespaceFooter.H"
#endif
