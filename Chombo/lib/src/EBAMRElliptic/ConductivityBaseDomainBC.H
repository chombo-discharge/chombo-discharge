#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ConductivityBASEDOMAINBC_H_
#define _ConductivityBASEDOMAINBC_H_

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
#include "NamespaceHeader.H"


///
/**
 */
class ConductivityBaseDomainBC: public BaseDomainBC
{
public:
  ///
  /**
   */
  ConductivityBaseDomainBC()
  {
    m_coefSet  = false;
  }

  void
  setValue(Real a_value)
  {
    m_onlyHomogeneous = false;
    m_isFunction = false;
    m_value = a_value;
    m_func = RefCountedPtr<BaseBCFuncEval>();
  }
  void
  setFunction(RefCountedPtr<BaseBCFuncEval> a_func)
  {
    m_value = 12345.6789;
    m_func = a_func;

    m_onlyHomogeneous = false;
    m_isFunction = true;
  }

  virtual ~ConductivityBaseDomainBC()
  {
  }

  ///
  /**
   */
  void setCoef(EBLevelGrid                         &  a_eblg,
               Real                                &  a_beta,
               RefCountedPtr<LevelData<EBFluxFAB> >&  a_bcoef)
  {
    m_coefSet = true;
    m_eblg    = a_eblg;
    m_beta    = a_beta;
    m_bcoef   = a_bcoef;
  }

  Real bcvaluefunc(const RealVect      & a_point,
                   const int           & a_dir,
                   const Side::LoHiSide& a_side)
 {
    Real retval;
    if(m_onlyHomogeneous)
      {
        retval = 0;
      }
    else if(m_isFunction)
      {
        retval = m_func->value(a_point, 0);
      }
    else
      {
        retval = m_value;
      }
    return retval;
  }

  //this makes applyOpRegular faster
  virtual void 
  fillPhiGhost(FArrayBox&     a_state,
               const Box&     a_valid,
               const Box&     a_domain,
               Real           a_dx,
               bool           a_homogeneous)
  {
    MayDay::Error("This shouldn't be called unless using NWO");
  }

protected:
  bool m_onlyHomogeneous;
  bool m_isFunction;

  Real m_value;
  RefCountedPtr<BaseBCFuncEval> m_func;

  bool                                  m_coefSet;
  EBLevelGrid                           m_eblg;
  Real                                  m_beta;
  RefCountedPtr<LevelData<EBFluxFAB> >  m_bcoef;
};
#include "NamespaceFooter.H"
#endif
