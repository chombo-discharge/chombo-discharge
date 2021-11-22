/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRParticleMesh.cpp
  @brief  Implementation of CD_EBAMRParticleMesh.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBAMRParticleMesh.H>
#include <CD_NamespaceHeader.H>

EBAMRParticleMesh::EBAMRParticleMesh(){
  CH_TIME("EBAMRParticleMesh::EBAMRParticleMesh()");

  m_isDefined = false;
}

EBAMRParticleMesh::EBAMRParticleMesh(const Vector<RefCountedPtr<EBLevelGrid> >& a_eblgs,
				     const Vector<int>&                         a_refRat,
				     const Vector<Real>&                        a_dx,
				     const RealVect&                            a_probLo,
				     const IntVect&                             a_ghost,	      
				     const int                                  a_maxParticleWidth,
				     const int                                  a_finestLevel){
  CH_TIME("EBAMRParticleMesh::EBAMRParticleMesh(full)");

  this->define(a_eblgs, a_refRat, a_dx, a_probLo, a_ghost, a_maxParticleWidth, a_finestLevel);
}

EBAMRParticleMesh::~EBAMRParticleMesh(){
  CH_TIME("EBAMRParticleMesh::~EBAMRParticleMesh()");
}

void EBAMRParticleMesh::define(const Vector<RefCountedPtr<EBLevelGrid> >& a_eblgs,
			       const Vector<int>&                         a_refRat,
			       const Vector<Real>&                        a_dx,
			       const RealVect&                            a_probLo,			       
			       const IntVect&                             a_ghost,
			       const int                                  a_maxParticleWidth,
			       const int                                  a_finestLevel){
  CH_TIME("EBAMRParticleMesh::define");

  m_eblgs            = a_eblgs;
  m_refRat           = a_refRat;
  m_dx               = a_dx;
  m_probLo           = a_probLo;
  m_ghost            = a_ghost;
  m_maxParticleWidth = a_maxParticleWidth;
  m_finestLevel      = a_finestLevel;

  this->defineLevelMotion     ();
  this->defineCoarseFineMotion();
  
  m_isDefined = true;
}

void EBAMRParticleMesh::defineLevelMotion(){
  CH_TIME("EBAMRParticleMesh::defineLevelMotion");


  // TLDR: Define level Copiers. These are defined such that we can move data from valid+ghost -> valid. We need this because when we deposit particles
  //       we will also deposit into ghost cells that overlap with patches on the same level. The data in those ghost cells needs to find its way into
  //       the neighboring patches, which is what these Copiers do. 
  m_levelCopiers.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const EBLevelGrid&       eblg   = *m_eblgs[lvl];
    const ProblemDomain&     domain =  eblg.getDomain();
    const DisjointBoxLayout& dbl    =  eblg.getDBL();

    // Define Copier. Note that the below code is basically the same as ghostDefine(). 
    const bool doExchange = true;           
    
    m_levelCopiers[lvl].define(dbl, dbl, domain, m_ghost, doExchange); // Define Copier as going from valid       -> valid+ghost.
    m_levelCopiers[lvl].reverse();                                     // Define Copier as going from valid+ghost -> valid.
  }
}

void EBAMRParticleMesh::defineCoarseFineMotion(){
  CH_TIME("EBAMRParticleMesh::defineCoarseFineMotion");

  m_coarseFinePM.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){

    const bool hasCoar = (lvl > 0);

    if(hasCoar){
      m_coarseFinePM[lvl] = RefCountedPtr<EBCoarseFineParticleMesh>(new EBCoarseFineParticleMesh(*m_eblgs [lvl-1],
												 *m_eblgs [lvl  ],
												 m_refRat[lvl-1],
												 m_ghost));
    }
    else{
      m_coarseFinePM[lvl] = RefCountedPtr<EBCoarseFineParticleMesh> (nullptr);
    }
  }
}

#include <CD_NamespaceFooter.H>
