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
			       const int                                  a_maxParticleWidth,
			       const int                                  a_finestLevel){
  CH_TIME("EBAMRParticleMesh::define");

  m_eblgs            = a_eblgs;
  m_refRat           = a_refRat;
  m_dx               = a_dx;
  m_maxParticleWidth = a_maxParticleWidth;
  m_finestLevel      = a_finestLevel;
  
  m_isDefined = true;
}

#include <CD_NamespaceFooter.H>
