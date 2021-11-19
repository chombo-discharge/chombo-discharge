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

EBAMRParticleMesh::~EBAMRParticleMesh(){
  CH_TIME("EBAMRParticleMesh::~EBAMRParticleMesh()");
}

void EBAMRParticleMesh::define(const Vector<RefCountedPtr<EBLevelGrid> >& a_eblgs,
			       const Vector<int>&                         a_refRat,
			       const Vector<Real>&                        a_dx,
			       const IntVect&                             a_ghost,
			       const int                                  a_maxParticleWidth,
			       const int                                  a_finestLevel){
  CH_TIME("EBAMRParticleMesh::define");

  m_eblgs            = a_eblgs;
  m_refRat           = a_refRat;
  m_dx               = a_dx;
  m_ghost            = a_ghost;
  m_maxParticleWidth = a_maxParticleWidth;
  m_finestLevel      = a_finestLevel;

  this->defineLevelCopiers();
  
  m_isDefined = true;
}

void EBAMRParticleMesh::defineLevelCopiers(){
  CH_TIME("EBAMRParticleMesh::defineLevelCopiers");


  // TLDR: Define level Copiers. These are defined such that we can move data from valid+ghost -> valid. We need this because when we deposit particles
  //       we will also deposit into ghost cells that overlap with patches on the same level. The data in those ghost cells needs to find its way into
  //       the neighboring patches, which is what these Copiers do. 
  m_levelCopiers.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const EBLevelGrid&       eblg   = *m_eblgs[lvl];
    const ProblemDomain&     domain =  eblg.getDomain();
    const DisjointBoxLayout& dbl    =  eblg.getDBL();

    // Define Copier. Note that the below code is basically the same as ghostDefine(). 
    const bool doExchange = false;           
    
    m_levelCopiers[lvl].define(dbl, dbl, domain, m_ghost, doExchange); // Define Copier as going from valid       -> valid+ghost.
    m_levelCopiers[lvl].reverse();                                     // Define Copier as going from valid+ghost -> valid.
  }
  
}

#include <CD_NamespaceFooter.H>
