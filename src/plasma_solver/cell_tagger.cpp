/*!
  @file cell_tagger.cpp
  @brief Implementation of cell_tagger.H
  @author Robert marskar
  @date Nov. 2017
*/

#include "cell_tagger.H"

cell_tagger::cell_tagger(const int a_num_tracers){
  m_num_tracers = a_num_tracers;
}

cell_tagger::~cell_tagger(){

}

void cell_tagger::define(const RefCountedPtr<plasma_kinetics>&        a_plaskin,
			 const RefCountedPtr<time_stepper>&           a_timestepper,
			 const RefCountedPtr<amr_mesh>&               a_amr,
			 const RefCountedPtr<computational_geometry>& a_compgeom,
			 const RefCountedPtr<physical_domain>&        a_physdom){

  m_plaskin     = a_plaskin;
  m_timestepper = a_timestepper;
  m_amr         = a_amr;
  m_compgeom    = a_compgeom;
  m_physdom     = a_physdom;
}



void cell_tagger::regrid(){

}

void cell_tagger::compute_tracers(){

}

void cell_tagger::compute_tracer_gradient(){

}

void cell_tagger::tag_cells(Vector<IntVectSet>& a_tags,
			    const Vector<RefCountedPtr<LayoutData<IntVectSet> > >& a_layout_tags,
			    const int a_finestLevel){

}

Vector<EBAMRCellData>& cell_tagger::get_tracer_fields() {
  return m_tracer;
}

Vector<Real> cell_tagger::tracer(const RealVect&         a_pos,
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
  Vector<Real> ret(m_num_tracers, 0.);
  return ret;
}

bool cell_tagger::coarsen_cell(const RealVect&         a_pos,
			       const Real&             a_time,
			       const Real&             a_dx,
			       const int&              a_lvl,
			       const Vector<Real>&     a_tracer,
			       const Vector<RealVect>& a_grad_tracer){
  return false;
}

bool cell_tagger::refine_cell(const RealVect&         a_pos,
			      const Real&             a_time,
			      const Real&             a_dx,
			      const int&              a_lvl,
			      const Vector<Real>&     a_tracer,
			      const Vector<RealVect>& a_grad_tracer){

}


