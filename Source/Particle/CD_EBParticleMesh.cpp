/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBParticleMesh.cpp
  @brief  Implementationof CD_EBParticleMesh.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBParticleMesh.H>
#include <CD_EBParticleMeshF_F.H>
#include <CD_NamespaceHeader.H>

EBParticleMesh::EBParticleMesh(){
  CH_TIME("EBParticleMesh::EBParticleMesh");
}

EBParticleMesh::EBParticleMesh(const Box&      a_region,
			       const EBISBox&  a_ebisbox,
			       const RealVect& a_dx,
			       const RealVect& a_probLo) {
  CH_TIME("EBParticleMesh::EBParticleMesh");
  
  this->define(a_region, a_ebisbox, a_dx, a_probLo);
}

void EBParticleMesh::define(const Box&      a_region,
			    const EBISBox&  a_ebisbox,
			    const RealVect& a_dx,
			    const RealVect& a_probLo) {
  CH_TIME("EBParticleMesh::define");
  
  m_region        = a_region;
  m_ebisbox       = a_ebisbox;
  m_probLo        = a_probLo;
  m_dx            = a_dx;
}

#include <CD_NamespaceFooter.H>
