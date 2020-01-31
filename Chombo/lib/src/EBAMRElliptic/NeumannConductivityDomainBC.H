#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEUMANNCONDUCTIVITYDOMAINBC_H_
#define _NEUMANNCONDUCTIVITYDOMAINBC_H_

#include "RefCountedPtr.H"
#include "ConductivityBaseDomainBC.H"
#include "NeumannPoissonDomainBC.H"
#include "BaseBCValue.H"
#include "NamespaceHeader.H"

class NeumannConductivityDomainBC: public ConductivityBaseDomainBC
{
public:
  NeumannConductivityDomainBC();

  virtual ~NeumannConductivityDomainBC();

  virtual int whichBC(int                  a_idir,
                      Side::LoHiSide       a_side);


  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

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

  virtual void getFaceGradPhi(Real&                 a_faceFlux,
                              const FaceIndex&      a_face,
                              const int&            a_comp,
                              const EBCellFAB&      a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_dx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useAreaFrac,
                              const RealVect&       a_centroid,
                              const bool&           a_useHomogeneous);

  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_face,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side);

  virtual void 
  fillPhiGhost(FArrayBox&     a_phi,
               const Box&     a_valid,
               const Box&     a_domain,
               Real           a_dx,
               bool           a_homogeneous)
  {
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
                    Real value = bcvaluefunc(loc, idir, sit());
                    if(a_homogeneous) value = 0;
                    a_phi(iv, 0) = a_phi(ivneigh, 0)  + a_dx*value;
                  }
              }
          } 
      }//end loop over directions
  }
private:


};

class NeumannConductivityDomainBCFactory: public BaseDomainBCFactory
{
public:
  NeumannConductivityDomainBCFactory();

  virtual ~NeumannConductivityDomainBCFactory();

  virtual void setValue(Real a_value);

  virtual void setFunction(RefCountedPtr<BaseBCFuncEval> a_flux);

  virtual NeumannConductivityDomainBC* create(const ProblemDomain& a_domain,
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
