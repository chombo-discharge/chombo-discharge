/*!
  @file   EBParticleInterpI.H
  @brief  Implementatoin of EBParticleInterp.H
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _EBPARTICLEINTERPI_H_
#define _EBPARTICLEINTERPI_H_

#include "CD_NamespaceHeader.H"
  
template <class P>
void EBParticleInterp::deposit(const List<P>& a_particleList, FArrayBox& a_rho, DepositionType::Which a_interpType){
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
void EBParticleInterp::deposit2(const List<P>& a_particleList, FArrayBox& a_rho, DepositionType::Which a_interpType){
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
void EBParticleInterp::deposit4(const List<P>& a_particleList, FArrayBox& a_rho, DepositionType::Which a_interpType){
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
void EBParticleInterp::interpolateDiffusion(List<P>&              a_particleList,
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
void EBParticleInterp::interpolateMobility(List<P>&              a_particleList,
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
void EBParticleInterp::interpolateVelocity(List<P>&              a_particleList,
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
void EBParticleInterp::interpolateAcceleration(List<P>&              a_particleList,
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
#include "CD_NamespaceFooter.H"

#endif
