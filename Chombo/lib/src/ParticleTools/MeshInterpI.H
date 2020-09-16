#ifndef _MESHINTERPI_H_
#define _MESHINTERPI_H_

#include "NamespaceHeader.H"

template <class P>
void MeshInterp::deposit(const List<P>& a_particleList,
			 FArrayBox& a_rho,
			 InterpType& a_interpType)
{
  for (ListIterator<P> lit(a_particleList); lit; ++lit)
  {
    const P& this_particle = lit();
    const RealVect& this_x = this_particle.position();
    const Real& this_m = this_particle.mass();
    depositParticle(a_rho, 
		    m_domainLeftEdge, 
		    m_dx, 
		    this_x, 
		    this_m, 
		    a_interpType);
  }
}

template <class P>
void MeshInterp::interpolate(List<P>& a_particleList,
			     const FArrayBox& a_field,
			     InterpType& a_interpType)
{
  for (ListIterator<P> lit(a_particleList); lit; ++lit)
    {
      RealVect particleField(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
      P& this_particle = lit();
      RealVect& this_x = this_particle.position();
      interpolateParticle(particleField,
			  a_field,
			  m_domainLeftEdge,
			  m_dx,
			  this_x,
			  a_interpType);
      this_particle.setAcceleration(particleField);
    }
}

#include "NamespaceFooter.H"

#endif // include guard
