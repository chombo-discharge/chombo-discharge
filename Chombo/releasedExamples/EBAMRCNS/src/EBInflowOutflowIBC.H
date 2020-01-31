#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef  _EB_InflowOutflowIBC_H__
#define  _EB_InflowOutflowIBC_H__

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "EBPhysIBC.H"
#include "EBInflowOutflowF_F.H"
#include "EBSolidF_F.H"

#include "NamespaceHeader.H"

class EBInflowOutflowIBC : public EBPhysIBC
{
public:
  virtual ~EBInflowOutflowIBC();

  EBInflowOutflowIBC(const Real&     a_gamma,
                     const Real&     a_cv,
                     const Real&     a_tempinflow,
                     const Real&     a_tempinterior,
                     const Real&     a_presinflow,
                     const Real&     a_presinterior,
                     const Real&     a_machinflow,
                     const int&      a_flowdir);

  ///
  void define(const ProblemDomain&  a_domain,
              const RealVect&       a_dx);

  ///  For every box in this level, this function is called.
  void fluxBC(EBFluxFAB&            a_flux,
              const EBCellFAB&      a_Wcenter,
              const EBCellFAB&      a_Wextrap,
              const Side::LoHiSide& a_sd,
              const Real&           a_time,
              const EBISBox&        a_ebisBox,
              const DataIndex&      a_dit,
              const Box&            a_box,
              const Box&            a_faceBox,
              const int&            a_dir);

  /// Initialize
  void initialize(LevelData<EBCellFAB>& a_conState,
                  const EBISLayout& a_ebisl) const;

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
  */
  void setBndrySlopes(EBCellFAB&       a_deltaPrim,
                      const EBCellFAB& a_primState,
                      const EBISBox&   a_ebisBox,
                      const Box&       a_box,
                      const int&       a_dir)
  {
    //one sided slopes at boundaries are fine for this
  }
protected:
  bool m_isFortranCommonSet;
  bool m_isDefined;

  Real     m_gamma         ;
  Real     m_cv            ;
  Real     m_tempinflow    ;
  Real     m_tempinterior  ;
  Real     m_presinflow    ;
  Real     m_presinterior  ;
  Real     m_machinflow    ;
  int      m_flowdir       ;

  Real m_dx;
  ProblemDomain  m_domain;

private:
  EBInflowOutflowIBC()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(const EBInflowOutflowIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBInflowOutflowIBC(const EBInflowOutflowIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }


};

#include "NamespaceFooter.H"
#endif
