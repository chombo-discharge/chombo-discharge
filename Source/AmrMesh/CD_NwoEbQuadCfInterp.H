#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*!
  @file   CD_NwoEbQuadCfInterp.H
  @brief  This is just a copy of Chombos NWOEBQuadCFInterp which also takes an ebis pointer argument. 
  @author Robert Marskar
  @note   Chombo copyright. 
*/

#ifndef CD_NwoEbQuadCfInterp_H
#define CD_NwoEbQuadCfInterp_H

// Std includes
#include <iostream>
#include <math.h>
#include <stdlib.h>

// Chombo includes
#include <SPACE.H>
#include <REAL.H>
#include <IntVect.H>
#include <Box.H>
#include <BaseFab.H>
#include <DisjointBoxLayout.H>
#include <LevelData.H>
#include <FArrayBox.H>
#include <QuadCFStencil.H>
#include <ProblemDomain.H>
#include <TensorFineStencilSet.H>
#include <REAL.H>
#include <AggStencil.H>
#include <FArrayBox.H>
#include <LevelData.H>
#include <DisjointBoxLayout.H>
#include <EBISLayout.H>
#include <EBCellFAB.H>
#include <Interval.H>
#include <Stencils.H>
#include <QuadCFInterp.H>
#include <BaseIVFAB.H>
#include <TensorCFInterp.H>
#include <CornerCopier.H>
#include <EBCFData.H>
#include <EBQuadCFInterp.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Reimplementation of NWOEBQuadCFInterp which takes an additional argument; the EBIndexSpace
*/
class NwoEbQuadCfInterp {
public:

  /*!
    @brief Constructor
    @param[in] a_gridsFine  Fine grids
    @param[in] a_gridsCoar  Coarse grids
    @param[in] a_ebislFine  Fine EBISLayout
    @param[in] a_ebislCoar  Coarse EBISLayout
    @param[in] a_domainCoar Coarse level domainw
    @param[in] a_nref       Refinement ratio
    @param[in] a_nvar       Number of variables
    @param[in] a_dxFine     Grid resolution on fine level
    @param[in] a_ghost      Ghost cells. 
    @param[in] a_cfivs      Coarse-fine interface
    @param[in] a_ebisPtr    EBIndexSpace pointer. 
   */
  NwoEbQuadCfInterp(const DisjointBoxLayout&           a_gridsFine,
		    const DisjointBoxLayout&           a_gridsCoar,
		    const EBISLayout&                  a_ebislFine,
		    const EBISLayout&                  a_ebislCoar,
		    const ProblemDomain&               a_domainCoar,
		    const int&                         a_nref,
		    const int&                         a_nvar,
		    const Real&                        a_dxFine,
		    const int&                         a_ghost,
		    const LayoutData<IntVectSet>&      a_cfivs,
		    const RefCountedPtr<EBIndexSpace>& a_ebisPtr);

  /*!
    @brief Destructor
  */
  ~NwoEbQuadCfInterp();



  /*!
     @brief Coarse fine interpolation
     @param[out] a_phif Fine level data
     @param[in]  a_phic Coarse level data
     @param[in]  isrc   Source component
     @param[in]  idst   Destination component
     @param[in]  inco   Number of components
  */
  void coarseFineInterp(LevelData<EBCellFAB>&       a_phif,
			const LevelData<EBCellFAB>& a_phic,
			int                         isrc,
			int                         idst,
			int                         inco);


  /*!
    @brief
     Coarse-fine interpolation operator with homogeneous BCs.
     Does same interpolation as coarseFineInterp function,
     but with all coarse-level data set to 0.
     @param[out] a_phif Fine level data
     @param[in]  isrc   Source component
     @param[in]  idst   Destination component
     @param[in]  inco   Number of components
  */
  void coarseFineInterpH(LevelData<EBCellFAB>& a_phif,
			 int                   isrc,
			 int                   idst,
			 int                   inco);


protected:
  void
  defineInternals(const LayoutData<IntVectSet>& a_cfivs, const RefCountedPtr<EBIndexSpace>& a_ebisPtr);
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
  int               m_ghost;
  LayoutData<RefCountedPtr< AggStencil<EBCellFAB, EBCellFAB> > > m_stencil;
  LevelData<EBCellFAB> m_bufferCoFi;
  

private:
  //weak construction is bad
  NwoEbQuadCfInterp() = delete;

  //disallowed for all the usual reasons
  NwoEbQuadCfInterp(NwoEbQuadCfInterp& a_input) = delete;

  //disallowed for all the usual reasons
  void operator=(NwoEbQuadCfInterp& a_input) = delete;
  
};

#include <CD_NamespaceFooter.H>

#endif
