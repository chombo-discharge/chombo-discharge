/*!
  @file   gradient_tagger.cpp
  @brief  Implementation of gradient_tagger.H
  @author Robert Marskar
  @date   July 2018
*/

#include "gradient_tagger.H"

#include <ParmParse.H>

gradient_tagger::gradient_tagger(){
  m_num_tracers = 0; // This is set by overriding the define function

  m_coar_curv = 0.1;
  m_refi_curv = 0.2;
  m_cutoff    = 1.E-6;

  { // Get input parameters
    ParmParse pp("gradient_tagger");
    pp.query("coarsen_curvature", m_coar_curv);
    pp.query("refine_curvature",  m_refi_curv);
    pp.query("realtive_cutoff",   m_cutoff);
  }
}

gradient_tagger::~gradient_tagger(){

}

void gradient_tagger::define(const RefCountedPtr<plasma_kinetics>&        a_plaskin,
			     const RefCountedPtr<time_stepper>&           a_timestepper,
			     const RefCountedPtr<amr_mesh>&               a_amr,
			     const RefCountedPtr<computational_geometry>& a_compgeom,
			     const RefCountedPtr<physical_domain>&        a_physdom){
  cell_tagger::define(a_plaskin, a_timestepper, a_amr, a_compgeom, a_physdom);

  m_num_tracers = a_plaskin->get_num_species();
}

Vector<Real> gradient_tagger::tracer(const RealVect&         a_pos,
				     const Real&             a_time,
				     const Real&             a_dx,
				     const RealVect&         a_E,
				     const Real&             a_min_E,
				     const Real&             a_max_E,
				     const RealVect&         a_grad_E,
				     const Real&             a_min_grad_E,
				     const Real&             a_max_grad_E,
				     const Real&             a_rho,
				     const Real&             a_min_rho,
				     const Real&             a_max_rho,
				     const RealVect&         a_grad_rho,
				     const Real&             a_min_grad_rho,
				     const Real&             a_max_grad_rho,
				     const Vector<Real>&     a_ion_densities,
				     const Vector<Real>&     a_min_ion_densities,
				     const Vector<Real>&     a_max_ion_densities,
				     const Vector<RealVect>& a_ion_gradients,
				     const Vector<Real>&     a_min_ion_gradients,
				     const Vector<Real>&     a_max_ion_gradients,
				     const Vector<Real>&     a_photon_densities,
				     const Vector<Real>&     a_min_photon_densities,
				     const Vector<Real>&     a_max_photon_densities){
  Vector<Real> tracers(m_num_tracers);
  for (int i = 0; i < tracers.size(); i++){
    const Real cur_phi = a_ion_densities[i];
    const Real cur_max = a_max_ion_densities[i];
    tracers[i] = (cur_phi > m_cutoff*cur_max) ? cur_phi : 0.0;
  }
  return tracers;
}


bool gradient_tagger::coarsen_cell(const RealVect&         a_pos,
				   const Real&             a_time,
				   const Real&             a_dx,
				   const int&              a_lvl,
				   const Vector<Real>&     a_tracer,
				   const Vector<RealVect>& a_grad_tracer){

  bool coarsen = false;
  for (int i = 0; i < a_tracer.size(); i++){
    const bool test = a_grad_tracer[i].vectorLength()*a_dx/a_tracer[i] < m_coar_curv || a_tracer[i] <= 0.0;

    if(test){
      coarsen = true;
    }
  }

  return coarsen;  
}

bool gradient_tagger::refine_cell(const RealVect&         a_pos,
				  const Real&             a_time,
				  const Real&             a_dx,
				  const int&              a_lvl,
				  const Vector<Real>&     a_tracer,
				  const Vector<RealVect>& a_grad_tracer){

  bool refine = false;
  for (int i = 0; i < a_tracer.size(); i++){
    const bool test = a_grad_tracer[i].vectorLength()*a_dx/a_tracer[i] > m_refi_curv && a_tracer[i] > 0.0;

    if(test){
      refine = true;
    }
  }

  return refine;
}
