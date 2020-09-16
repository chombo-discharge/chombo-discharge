#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBInflowOutflowIBCFactory.H"
#include "EBInflowOutflowIBC.H"

#include "NamespaceHeader.H"

/******************/
EBInflowOutflowIBCFactory::
EBInflowOutflowIBCFactory(const Real&     a_gamma,
                          const Real&     a_cv,
                          const Real&     a_tempinflow,
                          const Real&     a_tempinterior,
                          const Real&     a_presinflow,
                          const Real&     a_presinterior,
                          const Real&     a_machinflow,
                          const int&      a_flowdir)        
  :EBPhysIBCFactory()
{
  m_gamma         = a_gamma        ;
  m_cv            = a_cv           ;
  m_tempinflow    = a_tempinflow   ;
  m_tempinterior  = a_tempinterior ;
  m_presinflow    = a_presinflow   ;
  m_presinterior  = a_presinterior ;
  m_machinflow    = a_machinflow   ;
  m_flowdir       = a_flowdir      ;
}
/******************/
EBInflowOutflowIBCFactory::
~EBInflowOutflowIBCFactory()
{;}
/******************/
EBPhysIBC*
EBInflowOutflowIBCFactory::
create() const
{
  EBInflowOutflowIBC* retval =
    new EBInflowOutflowIBC(
                           m_gamma        ,
                           m_cv           ,
                           m_tempinflow   ,
                           m_tempinterior ,
                           m_presinflow   ,
                           m_presinterior ,
                           m_machinflow   ,
                           m_flowdir      );

  return static_cast<EBPhysIBC*>(retval);
}
/******************/

#include "NamespaceFooter.H"
