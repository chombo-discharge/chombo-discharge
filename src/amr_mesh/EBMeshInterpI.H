/*!
  @file   EBMeshInterpI.H
  @brief  Implementatoin of EBMeshInterp.H
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _EBMESHINTERPI_H_
#define _EBMESHINTERPI_H_


template <class P>
void EBMeshInterp::deposit(const List<P>&        a_particleList,
			   FArrayBox&            a_rho,
			   DepositionType::Which a_interpType){
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P& this_particle = lit();
    const RealVect& this_x = this_particle.position();
    const Real& this_m = this_particle.mass();
    depositParticle(a_rho, 
		    m_prob_lo, 
		    m_dx, 
		    this_x, 
		    this_m, 
		    a_interpType);
  }
}

template <class P>
void EBMeshInterp::interpolateDiffusion(List<P>&              a_particleList,
					const FArrayBox&      a_diffusionField,
					DepositionType::Which a_interpType){
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    Real particleDiffusion(0.0);
    P& p = lit();
    interpolateParticle(particleDiffusion,
			a_diffusionField,
			m_prob_lo,
			m_dx,
			p.position(),
			a_interpType);
    p.setDiffusion(particleDiffusion);
  }
}
  
template <class P>
void EBMeshInterp::interpolateVelocity(List<P>&              a_particleList,
				       const FArrayBox&      a_velocityField,
				       DepositionType::Which a_interpType){
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    RealVect particleVelocity(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    P& p = lit();
    interpolateParticle(particleVelocity,
			a_velocityField,
			m_prob_lo,
			m_dx,
			p.position(),
			a_interpType);
    p.setVelocity(particleVelocity);
  }
}

template <class P>
void EBMeshInterp::interpolateAcceleration(List<P>&              a_particleList,
					   const FArrayBox&      a_accelerationField,
					   DepositionType::Which a_interpType){
  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    RealVect particleAcceleration(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    P& p = lit();
    interpolateParticle(particleAcceleration,
			a_accelerationField,
			m_prob_lo,
			m_dx,
			p.position(),
			a_interpType);
    p.setAcceleration(particleAcceleration);
  }
}

#endif
