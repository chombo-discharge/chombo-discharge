/*!
  @file   air_11eed.cpp
  @brief  Implementation of air_11eed.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air_11eed.H"
#include "air_11eed_species.H"
#include "units.H"

#include <ParmParse.H>

air_11eed::air_11eed(){

}

air_11eed::~air_11eed(){

}

void air_11eed::get_gas_parameters(Real& a_Tg, Real& a_p, Real& a_N, Real& a_O2frac, Real& a_N2frac){
    ParmParse pp("air_11eed");
    pp.get("gas_temperature", a_Tg);
    pp.get("pressure", a_p);
    pp.get("O2_frac", a_O2frac);
    pp.get("N2_frac", a_N2frac);

    const Real tot_frac = a_O2frac + a_N2frac; 
    a_p      = a_p*units::s_atm2pascal;
    a_O2frac = a_O2frac/tot_frac; // Normalize to one
    a_N2frac = a_N2frac/tot_frac;
    a_N      = a_p*units::s_Na/(a_Tg*units::s_R);
}

Vector<Real> air_11eed::compute_cdr_diffusion_coefficients(const Real&         a_time,
							const RealVect&     a_pos,
							const RealVect&     a_E,
							const Vector<Real>& a_cdr_densities) const {

}

Vector<RealVect> air_11eed::compute_cdr_velocities(const Real&         a_time,
						const RealVect&     a_pos,
						const RealVect&     a_E,
						const Vector<Real>& a_cdr_densities) const {

}
  
Vector<Real> air_11eed::compute_cdr_source_terms(const Real              a_time,
					      const RealVect&         a_pos,
					      const RealVect&         a_E,
					      const RealVect&         a_gradE,
					      const Vector<Real>&     a_cdr_densities,
					      const Vector<Real>&     a_rte_densities,
					      const Vector<RealVect>& a_grad_cdr) const {

}

Vector<Real> air_11eed::compute_cdr_electrode_fluxes(const Real&         a_time,
						  const RealVect&     a_pos,
						  const RealVect&     a_normal,
						  const RealVect&     a_E,
						  const Vector<Real>& a_cdr_densities,
						  const Vector<Real>& a_cdr_velocities,
						  const Vector<Real>& a_rte_fluxes,
						  const Vector<Real>& a_extrap_cdr_fluxes) const {

}

Vector<Real> air_11eed::compute_cdr_dielectric_fluxes(const Real&         a_time,
						   const RealVect&     a_pos,
						   const RealVect&     a_normal,
						   const RealVect&     a_E,
						   const Vector<Real>& a_cdr_densities,
						   const Vector<Real>& a_cdr_velocities,
						   const Vector<Real>& a_rte_fluxes,
						   const Vector<Real>& a_extrap_cdr_fluxes) const {

}

Vector<Real> air_11eed::compute_rte_source_terms(const Real&         a_time,
					      const RealVect&     a_pos,
					      const RealVect&     a_E,
					      const Vector<Real>& a_cdr_densities) const {

}

Real air_11eed::initial_sigma(const Real a_time, const RealVect& a_pos) const {

}
