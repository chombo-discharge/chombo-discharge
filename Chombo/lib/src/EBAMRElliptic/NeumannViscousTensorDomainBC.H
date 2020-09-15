#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEUMANNVISCOUSTENSORDOMAINBC_H_
#define _NEUMANNVISCOUSTENSORDOMAINBC_H_

#include "RefCountedPtr.H"

#include "ViscousBaseDomainBC.H"
#include "BaseBCValue.H"
#include "BCFunc.H"
#include "NamespaceHeader.H"


class NeumannViscousTensorDomainBC: public ViscousBaseDomainBC
{
public:
  NeumannViscousTensorDomainBC();

  virtual ~NeumannViscousTensorDomainBC();

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
  fillVelGhost(FArrayBox&     a_phi,
               const Box&     a_valid,
               const Box&     a_domain,
               Real           a_dx,
               bool           a_homogeneous)
  {
    Box valid = a_valid;
    Box grownBox = a_valid;
    grownBox.grow(1);
    for (int idir=0; idir<CH_SPACEDIM; ++idir)
      {
        for(SideIterator sit; sit.ok(); ++sit)
          {
            Box choppedBox = grownBox;
            choppedBox.grow(idir,-1);
            Box toRegion = adjCellBox(choppedBox, idir, sit(), 1);

            if(!a_domain.contains(toRegion))
              {
                for (BoxIterator bit(toRegion); bit.ok(); ++bit)
                  {
                    const IntVect& iv = bit();
                    //fake vof just to get the location
                    VolIndex vof(iv, 0);
                    RealVect loc = EBArith::getVoFLocation(vof, a_dx, RealVect::Zero);
                    int isign = sign(sit());
                    IntVect ivneigh = iv - isign*BASISV(idir);
                    RealVect value = bcvaluefunc(loc, idir, sit());
                    if(a_homogeneous) value = RealVect::Zero;
                    for(int comp = 0; comp < SpaceDim; comp++)
                      {
                        //neumann
                        a_phi(iv, comp) = a_phi(ivneigh, comp)  + a_dx*value[comp];
                      }
                  }
              }
          } 
      }//end loop over directions
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

class NeumannViscousTensorDomainBCFactory: public BaseDomainBCFactory
{
public:
  NeumannViscousTensorDomainBCFactory();

  virtual ~NeumannViscousTensorDomainBCFactory();

  ///
  /**
     Sets value of NORMAL gradient of normal component.  Everything
     else in valueland is set to zero.   If you want something more
     elaborate, use the functional form.
  */
  virtual void setValue(Real a_value);

  virtual void setFunction(RefCountedPtr<BaseBCFuncEval> a_flux);

  virtual NeumannViscousTensorDomainBC* create(const ProblemDomain& a_domain,
                                               const EBISLayout&    a_layout,
                                               const RealVect&      a_dx);


private:
  bool m_onlyHomogeneous;
  bool m_isFunction;

  Real m_value;
  RefCountedPtr<BaseBCFuncEval> m_flux;
};

#include "NamespaceFooter.H"
#endif
