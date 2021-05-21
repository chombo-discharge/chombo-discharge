/*!
  @file   gradient_tagger.cpp
  @brief  Implementation of gradient_tagger.H
  @author Robert Marskar
  @date   July 2018
*/

#include "gradient_tagger.H"

#include <ParmParse.H>

gradient_tagger::gradient_tagger(){
  m_name = "gradient_tagger";
}

gradient_tagger::~gradient_tagger(){

}

void gradient_tagger::define(const RefCountedPtr<plasma_kinetics>&        a_plaskin,
			     const RefCountedPtr<TimeStepper>&           a_timeStepper,
			     const RefCountedPtr<AmrMesh>&               a_amr,
			     const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry,
			     const RefCountedPtr<physical_domain>&        a_physdom){
  CellTagger::define(a_plaskin, a_timeStepper, a_amr, a_computationalGeometry, a_physdom);

  m_num_tracers = a_plaskin->get_num_species();
}

void gradient_tagger::parseOptions(){
  parseVerbosity();
  parseBoxes();
  parseBuffer();

  ParmParse pp(m_name.c_str());
  pp.get("coarsen_curvature", m_coar_curv);
  pp.get("refine_curvature",  m_refi_curv);
  pp.get("relative_cutoff",   m_cutoff);
}

Vector<Real> gradient_tagger::tracer(const RealVect         a_pos,
				     const Real             a_time,
				     const Real             a_dx,
				     const RealVect         a_E,
				     const Real             a_min_E,
				     const Real             a_max_E,
				     const RealVect         a_grad_E,
				     const Real             a_min_grad_E,
				     const Real             a_max_grad_E,
				     const Real             a_rho,
				     const Real             a_min_rho,
				     const Real             a_max_rho,
				     const RealVect         a_grad_rho,
				     const Real             a_min_grad_rho,
				     const Real             a_max_grad_rho,
				     const Vector<Real>     a_ion_densities,
				     const Vector<Real>     a_min_ion_densities,
				     const Vector<Real>     a_max_ion_densities,
				     const Vector<RealVect> a_ion_gradients,
				     const Vector<Real>     a_min_ion_gradients,
				     const Vector<Real>     a_max_ion_gradients,
				     const Vector<Real>     a_Photon_densities,
				     const Vector<Real>     a_min_Photon_densities,
				     const Vector<Real>     a_max_Photon_densities){
  Vector<Real> tracers(m_num_tracers);
  for (int i = 0; i < tracers.size(); i++){
    const Real cur_phi = a_ion_densities[i];
    const Real cur_max = a_max_ion_densities[i];
    tracers[i] = (cur_phi > m_cutoff*cur_max) ? cur_phi : 0.0;
  }
  return tracers;
}


bool gradient_tagger::coarsenCell(const RealVect         a_pos,
				   const Real             a_time,
				   const Real             a_dx,
				   const int              a_lvl,
				   const Vector<Real>     a_tracer,
				   const Vector<RealVect> a_grad_tracer){

  bool coarsen = false;
  for (int i = 0; i < a_tracer.size(); i++){
    const bool test = a_grad_tracer[i].vectorLength()*a_dx/a_tracer[i] < m_coar_curv || a_tracer[i] <= 0.0;

    if(test){
      coarsen = true;
    }
  }

  return coarsen;  
}

bool gradient_tagger::refineCell(const RealVect         a_pos,
				  const Real             a_time,
				  const Real             a_dx,
				  const int              a_lvl,
				  const Vector<Real>     a_tracer,
				  const Vector<RealVect> a_grad_tracer){

  bool refine = false;
  for (int i = 0; i < a_tracer.size(); i++){
    const bool test = a_grad_tracer[i].vectorLength()*a_dx/a_tracer[i] > m_refi_curv && a_tracer[i] > 0.0;

    if(test){
      refine = true;
    }
  }

  return refine;
#include "CD_NamespaceFooter.H"
