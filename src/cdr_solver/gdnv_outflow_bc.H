#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _GDNV_OUTFLOW_BC_
#define _GDNV_OUTFLOW_BC_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"
#include "NamespaceHeader.H"

///uses extrapolated value in all cases
class gdnv_outflow_bc : public EBPhysIBC
{
public:
  virtual ~gdnv_outflow_bc();

  ///
  gdnv_outflow_bc();


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

  /// Initialize---no op here
  void initialize(LevelData<EBCellFAB>& a_conState,
                  const EBISLayout& a_ebisl) const;

  /// Set boundary slopes --no op here
  void setBndrySlopes(EBCellFAB&       a_deltaPrim,
                      const EBCellFAB& a_primState,
                      const EBISBox&   a_ebisBox,
                      const Box&       a_box,
                      const int&       a_dir)
  {;}


protected:
  bool     m_isDefined;
  RealVect m_dx;
  ProblemDomain  m_domain;

};


///
/**
 */
class gdnv_outflow_bcfactory: public EBPhysIBCFactory
{
public:

  ///
  ~gdnv_outflow_bcfactory();

  ///
  gdnv_outflow_bcfactory();


  ///
  EBPhysIBC* create() const;


};

#include "NamespaceFooter.H"
#endif
