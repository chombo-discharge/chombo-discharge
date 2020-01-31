#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INFLOWOUTFLOWVISCOUSTENSOREBBC_H_
#define _INFLOWOUTFLOWVISCOUSTENSOREBBC_H_

#include "RefCountedPtr.H"

#include "BaseEBBC.H"
#include "BaseBCFuncEval.H"
#include "LayoutData.H"
#include "IntVectSet.H"
#include "EBStencil.H"
#include "EBAMRCNSParams.H"
#include "DirichletViscousTensorEBBC.H"
#include "NamespaceHeader.H"
class VoFCache;
class VoFStencil;

///
/**
 */
class InflowOutflowViscousTensorEBBC: public DirichletViscousTensorEBBC
{
public:

  ///
  InflowOutflowViscousTensorEBBC();

  InflowOutflowViscousTensorEBBC(const Real & a_dirichletStart,
                                 const int&   a_flowDir,
                                 const ProblemDomain& a_domain,
                                 const EBISLayout&    a_layout,
                                 const RealVect&      a_dx);

  /// 
  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor);
  ///
  /**
   */
  virtual ~InflowOutflowViscousTensorEBBC();

  //always homogeneous so can simplify
  virtual void
  applyEBFlux(EBCellFAB&                    a_lphi,
              const EBCellFAB&              a_phi,
              VoFIterator&                  a_vofit,
              const LayoutData<IntVectSet>& a_cfivs,
              const DataIndex&              a_dit,
              const RealVect&               a_probLo,
              const RealVect&               a_dx,
              const Real&                   a_factor,
              const bool&                   a_useHomogeneous,
              const Real&                   a_time)
  {
    return;
  }

protected:

  Real m_dirichletStart;
  int  m_flowDir;

};

///
/**
 */
class InflowOutflowViscousTensorEBBCFactory: public BaseEBBCFactory
{
public:
  ///
  /**
   */
  InflowOutflowViscousTensorEBBCFactory(const Real & a_dirichletStart,
                                        const int&   a_flowDir)
  {
    m_dirichletStart = a_dirichletStart;
    m_flowDir        = a_flowDir;
  }

  ///
  /**
   */
  virtual ~InflowOutflowViscousTensorEBBCFactory()
  {; }


  ///
  /**
   */
  virtual InflowOutflowViscousTensorEBBC* create(  const ProblemDomain& a_domain,
                                                   const EBISLayout&    a_layout,
                                                   const RealVect&      a_dx,
                                                   const IntVect*       a_ghostCellsPhi,
                                                   const IntVect*       a_ghostCellsRhs)

  {
    CH_TIME("InflowOutflowViscousTensorEBBC::create");
    InflowOutflowViscousTensorEBBC* fresh = new InflowOutflowViscousTensorEBBC(m_dirichletStart, m_flowDir,
                                                                               a_domain,a_layout,a_dx);
    return fresh;
  }
private:

  Real m_dirichletStart;
  int  m_flowDir;

};



#include "NamespaceFooter.H"
#endif
