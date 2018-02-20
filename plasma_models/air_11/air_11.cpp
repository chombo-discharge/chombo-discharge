/*!
  @file   air_11.cpp
  @brief  Implementation of air_11.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air_11.H"
#include "air_11_species.H"

air_11::air_11(){

}

air_11::~air_11(){

}

Vector<Real> air_11::compute_cdr_diffusion_coefficients(const Real&         a_time,
							const RealVect&     a_pos,
							const RealVect&     a_E,
							const Vector<Real>& a_cdr_densities) const {

}

Vector<RealVect> air_11::compute_cdr_velocities(const Real&         a_time,
						const RealVect&     a_pos,
						const RealVect&     a_E,
						const Vector<Real>& a_cdr_densities) const {

}
  
Vector<Real> air_11::compute_cdr_source_terms(const Real              a_time,
					      const RealVect&         a_pos,
					      const RealVect&         a_E,
					      const RealVect&         a_gradE,
					      const Vector<Real>&     a_cdr_densities,
					      const Vector<Real>&     a_rte_densities,
					      const Vector<RealVect>& a_grad_cdr) const {

}

Vector<Real> air_11::compute_cdr_electrode_fluxes(const Real&         a_time,
						  const RealVect&     a_pos,
						  const RealVect&     a_normal,
						  const RealVect&     a_E,
						  const Vector<Real>& a_cdr_densities,
						  const Vector<Real>& a_cdr_velocities,
						  const Vector<Real>& a_rte_fluxes,
						  const Vector<Real>& a_extrap_cdr_fluxes) const {

}

Vector<Real> air_11::compute_cdr_dielectric_fluxes(const Real&         a_time,
						   const RealVect&     a_pos,
						   const RealVect&     a_normal,
						   const RealVect&     a_E,
						   const Vector<Real>& a_cdr_densities,
						   const Vector<Real>& a_cdr_velocities,
						   const Vector<Real>& a_rte_fluxes,
						   const Vector<Real>& a_extrap_cdr_fluxes) const {

}

Vector<Real> air_11::compute_rte_source_terms(const Real&         a_time,
					      const RealVect&     a_pos,
					      const RealVect&     a_E,
					      const Vector<Real>& a_cdr_densities) const {

}

Real air_11::initial_sigma(const Real a_time, const RealVect& a_pos) const {

}
