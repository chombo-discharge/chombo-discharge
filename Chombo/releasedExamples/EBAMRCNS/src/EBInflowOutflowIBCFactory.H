#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#ifndef _EBInflowOutflowIBCFACTORY_H
#define _EBInflowOutflowIBCFACTORY_H

#include "EBInflowOutflowIBC.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "EBPhysIBC.H"
#include "EBInflowOutflowIBC.H"
#include "EBPhysIBCFactory.H"

#include "NamespaceHeader.H"

///
/**
 */
class EBInflowOutflowIBCFactory: public EBPhysIBCFactory
{
public:

  ///
  ~EBInflowOutflowIBCFactory();

  ///
  EBInflowOutflowIBCFactory(const Real&     a_gamma,
                            const Real&     a_cv,
                            const Real&     a_tempinflow,
                            const Real&     a_tempinterior,
                            const Real&     a_presinflow,
                            const Real&     a_presinterior,
                            const Real&     a_machinflow,
                            const int&      a_flowdir);

  ///
  EBPhysIBC* create() const;

protected:

  Real     m_gamma         ;
  Real     m_cv            ;
  Real     m_tempinflow    ;
  Real     m_tempinterior  ;
  Real     m_presinflow    ;
  Real     m_presinterior  ;
  Real     m_machinflow    ;
  int      m_flowdir       ;


private:
  EBInflowOutflowIBCFactory()
  {
    MayDay::Error("Invalid operator");
  }

  void operator=(const EBInflowOutflowIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

  EBInflowOutflowIBCFactory(const EBInflowOutflowIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }


};

#include "NamespaceFooter.H"
#endif
