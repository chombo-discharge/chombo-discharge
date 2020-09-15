#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VISCOUSBASEDOMAINBC_H_
#define _VISCOUSBASEDOMAINBC_H_

#include "LoHiSide.H"
#include "RealVect.H"
#include "FArrayBox.H"
#include "BaseDomainBC.H"
#include "VolIndex.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "EBFluxFAB.H"
#include "EBISLayout.H"
#include "EBLevelGrid.H"
#include "TensorCFInterp.H"
#include "BaseBCFuncEval.H"
#include "ViscousTensorOp.H"
#include "NamespaceHeader.H"


///
/**
 */
class ViscousBaseDomainBC: public BaseDomainBC
{
public:
  ///
  /**
   */
  ViscousBaseDomainBC()
  {
    m_coefSet  = false;
    m_value = 12345.6789;
    m_func = RefCountedPtr<BaseBCFuncEval>();
    m_isFunction = false;
  }

  virtual ~ViscousBaseDomainBC()
  {
  }
  ///
  /**
   */
  void setCoef(EBLevelGrid                         &  a_eblg,
               Real                                &  a_beta,
               RefCountedPtr<LevelData<EBFluxFAB> >&  a_eta,
               RefCountedPtr<LevelData<EBFluxFAB> >&  a_lambda)
  {
    m_coefSet = true;
    m_eblg    = a_eblg;
    m_beta    = a_beta;
    m_eta     = a_eta;
    m_lambda  = a_lambda;
  }

  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_vof,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side)
  {
    //projection thingy does not apply here
    MayDay::Error();
  }
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
                              const bool&           a_useHomogeneous)
  {
    //projection thingy does not apply here
    MayDay::Error();
  }
  ///
  /**
   */
  virtual void setValue(Real a_value)
  {
    m_value = a_value;
    m_func = RefCountedPtr<BaseBCFuncEval>();

    m_isFunction = false;
  }


  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCFuncEval> a_func)
  {
    m_value = 12345.6789;
    m_func = a_func;

    m_isFunction = true;
  }

  void getFluxFromGrad(BaseFab<Real>&   a_flux,
                       const FArrayBox& a_grad,
                       const DataIndex& a_dit,
                       const int&       a_idir)
  {
    FArrayBox& fluxFAB = (FArrayBox&) a_flux;
    FArrayBox faceDiv(a_flux.box(), 1);
    faceDiv.setVal(0.);
    //compute the derivs as the sum of the appropriate grads
    for (int divDir = 0; divDir < SpaceDim; divDir++)
      {
        int gradComp = TensorCFInterp::gradIndex(divDir,divDir);
        int srccomp = gradComp;
        int dstcomp = 0;
        faceDiv.plus(a_grad, srccomp, dstcomp, 1);
      }

    //need to do this because there is an increment later
    const Box& faceBox = fluxFAB.box();
    fluxFAB.setVal(0.);
    const FArrayBox& lamFace = (const FArrayBox&)((*m_lambda)[a_dit][a_idir].getSingleValuedFAB());
    const FArrayBox& etaFace = (const FArrayBox&)((*m_eta   )[a_dit][a_idir].getSingleValuedFAB());
    ViscousTensorOp::getFluxFromDivAndGrad(fluxFAB, faceDiv, a_grad, etaFace, lamFace, faceBox, a_idir);
  }

  //this makes applyOpRegular faster
  virtual void 
  fillVelGhost(FArrayBox&     a_state,
               const Box&     a_valid,
               const Box&     a_domain,
               Real           a_dx,
               bool           a_homogeneous) = 0;


  RealVect bcvaluefunc(const RealVect      & a_point,
                       const int           & a_dir,
                       const Side::LoHiSide& a_side)
  {
    RealVect retval;
    if(m_isFunction)
      {
        for(int idir = 0; idir < SpaceDim; idir++)
          {
            retval[idir]=  m_func->value(a_point, idir);
          }
      }
    else
      {
        for(int idir = 0; idir < SpaceDim; idir++)
          {
            retval[idir] = m_value;
          }
      }
    return retval;
  }
protected:
  bool m_onlyHomogeneous;
  bool m_isFunction;

  Real m_value;
  RefCountedPtr<BaseBCFuncEval> m_func;


  Real                                  m_beta;
  bool                                  m_coefSet;
  EBLevelGrid                           m_eblg;
  RefCountedPtr<LevelData<EBFluxFAB> >  m_eta;
  RefCountedPtr<LevelData<EBFluxFAB> >  m_lambda;
};

#include "NamespaceFooter.H"
#endif
