/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_EbParticleInterpImplem.H
  @brief  Implementatoin of CD_EbParticleInterp.H
  @author Robert Marskar
*/

#ifndef CD_EbParticleInterpImplem_H
#define CD_EbParticleInterpImplem_H

#include <CD_NamespaceHeader.H>
  
template <class P>
void EbParticleInterp::deposit(const List<P>& a_particleList, FArrayBox& a_rho, DepositionType::Which a_interpType){
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P& this_particle = lit();
    const RealVect& this_x = this_particle.position();
    const Real& this_m     = this_particle.mass();
    depositParticle(a_rho, 
		    m_probLo, 
		    m_dx, 
		    this_x, 
		    this_m, 
		    a_interpType);
  }
}

template <class P>
void EbParticleInterp::deposit2(const List<P>& a_particleList, FArrayBox& a_rho, DepositionType::Which a_interpType){
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P& this_particle = lit();
    const RealVect& this_x = this_particle.position();
    const Real& this_m     = this_particle.mass();
    depositParticle2(a_rho, 
		     m_probLo, 
		     m_dx, 
		     this_x, 
		     this_m, 
		     a_interpType);
  }
}

template <class P>
void EbParticleInterp::deposit4(const List<P>& a_particleList, FArrayBox& a_rho, DepositionType::Which a_interpType){
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P& this_particle = lit();
    const RealVect& this_x = this_particle.position();
    const Real& this_m     = this_particle.mass();
    depositParticle4(a_rho, 
		     m_probLo, 
		     m_dx, 
		     this_x, 
		     this_m, 
		     a_interpType);
  }
}

template <class P>
void EbParticleInterp::interpolateDiffusion(List<P>&              a_particleList,
					    const FArrayBox&      a_diffusionField,
					    DepositionType::Which a_interpType){
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    Real particleDiffusion(0.0);
    P& p = lit();
    interpolateParticle(particleDiffusion,
			a_diffusionField,
			m_probLo,
			m_dx,
			p.position(),
			a_interpType);
    p.setDiffusion(particleDiffusion);
  }
}

template <class P>
void EbParticleInterp::interpolateMobility(List<P>&              a_particleList,
					   const FArrayBox&      a_mobilityField,
					   DepositionType::Which a_interpType){
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    Real particleMobility(0.0);
    P& p = lit();
    interpolateParticle(particleMobility,
			a_mobilityField,
			m_probLo,
			m_dx,
			p.position(),
			a_interpType);
    p.setMobility(particleMobility);
  }
}
  
template <class P>
void EbParticleInterp::interpolateVelocity(List<P>&              a_particleList,
					   const FArrayBox&      a_velocityField,
					   DepositionType::Which a_interpType){
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    RealVect particleVelocity(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    P& p = lit();
    interpolateParticle(particleVelocity,
			a_velocityField,
			m_probLo,
			m_dx,
			p.position(),
			a_interpType);
    p.setVelocity(particleVelocity);
  }
}

template <class P>
void EbParticleInterp::interpolateAcceleration(List<P>&              a_particleList,
					       const FArrayBox&      a_accelerationField,
					       DepositionType::Which a_interpType){
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    RealVect particleAcceleration(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    P& p = lit();
    interpolateParticle(particleAcceleration,
			a_accelerationField,
			m_probLo,
			m_dx,
			p.position(),
			a_interpType);
    p.setAcceleration(particleAcceleration);
  }
}

#include <CD_NamespaceFooter.H>

#endif
