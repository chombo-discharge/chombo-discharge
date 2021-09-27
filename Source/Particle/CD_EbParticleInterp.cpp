/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbParticleInterp.cpp
  @brief  Implementationof CD_EbParticleInterp.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EbParticleInterp.H>
#include <CD_EbParticleInterpF_F.H>
#include <CD_NamespaceHeader.H>

EbParticleInterp::EbParticleInterp(){
  CH_TIME("EbParticleInterp::EbParticleInterp");
}

EbParticleInterp::EbParticleInterp(const Box&      a_region,
				   const EBISBox&  a_ebisbox,
				   const RealVect& a_dx,
				   const RealVect& a_probLo,
				   const bool      a_forceIrregNGP){
  CH_TIME("EbParticleInterp::EbParticleInterp");
  
  this->define(a_region, a_ebisbox, a_dx, a_probLo, a_forceIrregNGP);
}

void EbParticleInterp::define(const Box&      a_region,
			      const EBISBox&  a_ebisbox,
			      const RealVect& a_dx,
			      const RealVect& a_probLo,
			      const bool      a_forceIrregNGP){
  CH_TIME("EbParticleInterp::define");
  
  m_region        = a_region;
  m_ebisbox       = a_ebisbox;
  m_probLo        = a_probLo;
  m_dx            = a_dx;
  m_forceIrregNGP = a_forceIrregNGP;
}

#include <CD_NamespaceFooter.H>
