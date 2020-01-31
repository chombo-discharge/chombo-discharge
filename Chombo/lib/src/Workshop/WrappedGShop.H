#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _WRAPPEDGSHOP_H_
#define _WRAPPEDGSHOP_H_

#include "REAL.H"
#include "RealVect.H"
#include "Box.H"
#include "IntVect.H"
#include "EBISBox.H"
#include "GeometryService.H"
#include "WrappedGShop.H"
#include "LSquares.H"
#include "BaseIF.H"
#include "AMRIO.H"
#include "IrregNode.H"
#include "CutCellMoments.H"
#include "RefCountedPtr.H"
#include "NamespaceHeader.H"
///base class for special refinement
class WGSRefinementCriterion
{
public:
  //default says no special refinement at any point
  virtual bool refineHere(const IrregNode& a_node, const IntVect& a_iv, const Real& a_dx) const
  {
    return false;
  }

  virtual ~WGSRefinementCriterion()
  {
  }

  WGSRefinementCriterion()
  {
  }
};
///
/**
   This is a much simplified version of GeometryShop.
   Refinement has been taken out of the internals of computecutcellmoments and
   the constrained least squares stuff is gone alltogether.   
   This is also done for isotropic dx and there is no more stuff
   about globaldim.   
*/
class WrappedGShop: public GeometryService
{

public:

  int m_phase;

  ///
  /**
     This class will refine a cell if it is between min and max refinements.
     If it violates bounds, that triggers refinement until max refinement.
  */
  WrappedGShop(const RefCountedPtr<BaseIF>  & a_baseIF,
               const RealVect               & a_origin,
               const Real                   & a_dx,
               const ProblemDomain          & a_domain,
               int   minNumberRefines,
               int   maxNumberRefines);
  


  ///
  /**  
     This checks to see if moments are within sane bounds.
     If a_bindMoments = true, it enforces those bounds.  
     Only returns true if the bounds are violated by the relative tolerance.
  */
  bool checkNodeMoments(IrregNode & a_node, 
                        const Real& a_dx,
                        const bool& a_bindMoments,
                        const Real& a_tolerance) const;

  virtual const BaseIF* getBaseIFPtr() const
  {
    //used for exact geometry generation. 
    return &(*m_baseIF);
  }


  ///
  ~WrappedGShop()
  {
  }

  /**
     Return true if every cell in region is regular at the
     refinement described by dx.
  */
  virtual bool isRegular(const Box           & a_region,
                         const ProblemDomain & a_domain,
                         const RealVect      & a_origin,
                         const Real          & a_dx) const;

  ///
  /**
     Return true if every cell in region is covered at the
     refinement described by dx.
  */
  virtual bool isCovered(const Box           & a_region,
                         const ProblemDomain & a_domain,
                         const RealVect      & a_origin,
                         const Real          & a_dx) const;


  virtual bool isIrregular(const Box           & a_region,
                           const ProblemDomain & a_domain,
                           const RealVect      & a_origin,
                           const Real          & a_dx) const;

  virtual bool canGenerateMultiCells() const
  {
    return false;
  }

  bool onBoxBoundary(const IntVect        & a_iv, 
                     const Box            & a_box,
                     const int            & a_dir,
                     const Side::LoHiSide & a_sd) const;

  ///
  /**
     Define the internals of the input ebisRegion.
  */
  virtual void fillGraph(BaseFab<int>        & a_regIrregCovered,
                         Vector<IrregNode>   & a_nodes,
                         const Box           & a_validRegion,
                         const Box           & a_ghostRegion,
                         const ProblemDomain & a_domain,
                         const RealVect      & a_origin,
                         const Real          & a_dx) const;

  /**
   */
  void
  computeVoFInternals(IrregNode                &     a_node,
                      const IntVectSet         &     a_ivsIrreg,
                      const ProblemDomain      &     a_domain,
                      const RealVect           &     a_origin,
                      const Real               &     a_dx,
                      const IntVect            &     a_iv) const;



  //will bound moments and check the number of refinements
  bool needToRefine(IrregNode       & a_node, 
                    const Real      & a_dx,  
                    const int       & a_numRefSoFar)const ;

  void agglomerateMoments(IrregNode              & a_node, 
                          const Vector<IrregNode>& a_refNodes,
                          const Box              & a_refBox,
                          const Real             & a_fineDx,
                          const Real             & a_coarDx) const;

  //converts a RealVect in physical coordinates to a RealVect in coordinates relative to a cell center
  RealVect convert2RelativeCoord(const RealVect& a_rVect)const;

  //converts a IndexTM<Real,SpaceDim> in physical coordinates to a
  //RealVect in coordinates relative to a cell center
  RealVect convert2RelativeCoord(const IndexTM<Real,SpaceDim>& a_rVect)const;

  void
  fillNewNode(IrregNode                &     a_node,
              const IntVectSet         &     a_ivsIrreg,
              const ProblemDomain      &     a_domain,
              const RealVect           &     a_origin,
              const Real               &     a_dx,
              const IntVect            &     a_iv) const;

  void setRefinementCriterion(const RefCountedPtr<WGSRefinementCriterion>& a_refCrit)
  {
    m_refCrit = a_refCrit;
  }

private:

  RefCountedPtr<WGSRefinementCriterion> m_refCrit;
  
  Real                m_threshold;

  //origin
  RealVect            m_origin;

  int                 m_order;
  int                 m_degreeP;
  int                 m_minNumberRefines;
  int                 m_maxNumberRefines;

  ProblemDomain m_domain;

  RefCountedPtr<BaseIF>       m_baseIF;
  
  static Real s_relativeTol;
  /**
     Return true if every cell in region is regular at the
     refinement described by dx.
  */
  bool isRegularEveryPoint(const Box&           a_region,
                           const ProblemDomain& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;

  ///
  /**
     Return true if every cell in region is covered at the
     refinement described by dx.
  */
  bool isCoveredEveryPoint(const Box&           a_region,
                           const ProblemDomain& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;


  void fillArc(Vector<int>                          a_arc[SpaceDim],
               CutCellMoments<SpaceDim>       &     a_cutCellMoments,
               const int                      &     a_hilo,
               const IntVectSet               &     a_ivsIrreg,
               const IntVect                  &     a_curriv) const;


  //stuff disallowed for all the usual reasons.
  WrappedGShop()
  {
    MayDay::Abort("WrappedGShop uses strong construction only");
  }
  WrappedGShop(const WrappedGShop& a_workshopin)
  {
    MayDay::Abort("WrappedGShop disallows copy contruction");
  }
  void operator=(const WrappedGShop& a_workshopin)
  {
    MayDay::Abort("WrappedGShop disallows the assignment operator");
  }

};
#include "NamespaceFooter.H"
#endif
