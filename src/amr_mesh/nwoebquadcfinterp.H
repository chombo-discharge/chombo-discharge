#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _nwoebquadcfinterp_H_
#define _nwoebquadcfinterp_H_

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "BaseFab.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "QuadCFStencil.H"
#include "ProblemDomain.H"
#include "TensorFineStencilSet.H"
#include "REAL.H"
#include "AggStencil.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "Stencils.H"
#include "QuadCFInterp.H"
#include "BaseIVFAB.H"
#include "TensorCFInterp.H"
#include "CornerCopier.H"
#include "EBCFData.H"
#include "EBQuadCFInterp.H"
#include "NamespaceHeader.H"

/// Quadratic coarse-fine interpolation utility 
class nwoebquadcfinterp
{
public:

  ///
  nwoebquadcfinterp(const DisjointBoxLayout& a_gridsFine,
                    const DisjointBoxLayout& a_gridsCoar,
                    const EBISLayout&        a_ebislFine,
                    const EBISLayout&        a_ebislCoar,
                    const ProblemDomain&     a_domainCoar,
                    const int&               a_nref,
                    const int&               a_nvar,
                    const Real&              a_dxFine,
                    const IntVect&           a_ghost,
                    const LayoutData<IntVectSet>&  a_cfivs,
		    const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());

  ///
  ~nwoebquadcfinterp()
  {  }

  ///
  /**
     Coarse / Fine interpolation operator.
  */
  void coarseFineInterp(LevelData<EBCellFAB>&       a_phif,
                        const LevelData<EBCellFAB>& a_phic,
                        int isrc, int idst, int inco);


  ///
  /**
     Coarse-fine interpolation operator with homogeneous BCs.
     Does same interpolation as coarseFineInterp function,
     but with all coarse-level data set to 0.
  */
  void coarseFineInterpH(LevelData<EBCellFAB>& a_phif,
                        int isrc, int idst, int inco);


protected:
  void
  defineInternals(const LayoutData<IntVectSet>& a_cfivs, const EBIndexSpace* const a_ebisPtr);
  void
  defineStencils (const LayoutData<IntVectSet>& a_cfivs);
  
  void
  getStencils(Vector<RefCountedPtr< BaseStencil> >  & a_stencils, 
              Vector<RefCountedPtr< BaseIndex  > >  & a_baseDstVoFs,
              const IntVectSet                      & a_cfivs,
              const EBISBox                         & a_ebisFine,
              const EBISBox                         & a_ebisCoFi,
              const DataIndex                       & a_dit);


  void
  getStencil(VoFStencil           & a_stencil,
             const VolIndex       & a_vofFine,
             const EBISBox        & a_ebisFine,
             const EBISBox        & a_ebisCoFi,
             const DataIndex      & a_dit);

  DisjointBoxLayout m_gridsFine;
  DisjointBoxLayout m_gridsCoar;
  DisjointBoxLayout m_gridsCoFi;
  EBISLayout        m_ebislFine;
  EBISLayout        m_ebislCoar;
  EBISLayout        m_ebislCoFi;
  ProblemDomain     m_domainCoar;
  int               m_nref;
  int               m_nvar;
  Real              m_dxFine;
  IntVect           m_ghost;
  LayoutData<RefCountedPtr< AggStencil<EBCellFAB, EBCellFAB> > > m_stencil;
  LevelData<EBCellFAB> m_bufferCoFi;
  

private:
  //weak construction is bad
  nwoebquadcfinterp()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  nwoebquadcfinterp(nwoebquadcfinterp& a_input)
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(nwoebquadcfinterp& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
