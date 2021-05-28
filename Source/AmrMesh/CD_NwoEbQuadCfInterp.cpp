/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NwoEbQuadCfInterp.cpp
  @brief  Implementation of CD_NwoEbQuadCfInterp.H
  @author Robert Marskar
*/


// Chombo includes
#include <EBCellFactory.H>

// Our includes
#include <CD_NwoEbQuadCfInterp.H>
#include <CD_NamespaceHeader.H>

NwoEbQuadCfInterp::NwoEbQuadCfInterp(const DisjointBoxLayout&           a_gridsFine,
				     const DisjointBoxLayout&           a_gridsCoar,
				     const EBISLayout&                  a_ebislFine,
				     const EBISLayout&                  a_ebislCoar,
				     const ProblemDomain&               a_domainCoar,
				     const int&                         a_nref,
				     const int&                         a_nvar,
				     const Real&                        a_dxFine,
				     const IntVect&                     a_ghost,
				     const LayoutData<IntVectSet>&      a_cfivs,
				     const RefCountedPtr<EBIndexSpace>& a_ebisPtr) {
  m_gridsFine  =  a_gridsFine ;          
  m_gridsCoar  =  a_gridsCoar ;     
  m_ebislFine  =  a_ebislFine ;      
  m_ebislCoar  =  a_ebislCoar ;      
  m_domainCoar =  a_domainCoar;
  m_nref       =  a_nref      ;      
  m_nvar       =  a_nvar      ;      
  m_dxFine     =  a_dxFine    ;      
  m_ghost      =  a_ghost     ;
  
  defineInternals(a_cfivs, a_ebisPtr);
}

NwoEbQuadCfInterp::~NwoEbQuadCfInterp(){

}

void NwoEbQuadCfInterp::defineInternals(const LayoutData<IntVectSet>&  a_cfivs, const RefCountedPtr<EBIndexSpace>& a_ebisPtr) {
  m_gridsCoFi= DisjointBoxLayout();

  coarsen(m_gridsCoFi, m_gridsFine, m_nref);

  a_ebisPtr->fillEBISLayout(m_ebislCoFi,
			    m_gridsCoFi,
			    m_domainCoar, 4); 

  EBCellFactory fact(m_ebislCoFi);
  m_bufferCoFi.define(m_gridsCoFi, m_nvar, 4*IntVect::Unit, fact);

  defineStencils(a_cfivs);
}

#include <CD_NamespaceFooter.H>
