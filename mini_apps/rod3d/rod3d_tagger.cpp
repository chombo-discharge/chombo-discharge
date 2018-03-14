/*!
  @file   rod3d_tagger.cpp
  @brief  Implementation rod3d_tagger.cpp
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "rod3d_tagger.H"

#include <ParmParse.H>
#include <PolyGeom.H>

rod3d_tagger::rod3d_tagger(){
  m_num_tracers = 1;

  m_coar_curv = 0.1;
  m_coar_mag  = 0.1;
  m_refi_curv = 0.2;
  m_refi_mag  = 0.75;

  // Get input parameters if we have them
  {
    ParmParse pp("rod3d_tagger");

    pp.query("coarsen_curvature", m_coar_curv);
    pp.query("coarsen_magnitude", m_coar_mag);
    pp.query("refine_curvature",  m_refi_curv);
    pp.query("refine_magnitude",  m_refi_mag);
  }

  { // Geometrical inputs
    Vector<Real> vec1, vec2;
    ParmParse pp("rod_slab");
    pp.get("electrode_radius", m_rod_radius);
    pp.getarr("electrode_center1", vec1, 0, SpaceDim);
    pp.getarr("electrode_center2", vec2, 0, SpaceDim);
    
    const RealVect c1   = RealVect(D_DECL(vec1[0], vec1[1], vec1[2]));
    const RealVect c2   = RealVect(D_DECL(vec2[0], vec2[1], vec2[2]));
    
    m_axis       = (c1 - c2)/(c1-c2).vectorLength();
    m_rod_center = c2 + m_axis*m_rod_radius;
  }
}


rod3d_tagger::~rod3d_tagger(){

}


Vector<Real> rod3d_tagger::tracer(const RealVect&         a_pos,
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
  tracers[0] = a_E.vectorLength()/a_max_E;

  return tracers;
}


bool rod3d_tagger::coarsen_cell(const RealVect&         a_pos,
				const Real&             a_time,
				const Real&             a_dx,
				const int&              a_lvl,
				const Vector<Real>&     a_tracer,
				const Vector<RealVect>& a_grad_tracer){

  const bool coarsen = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] < m_coar_curv && a_tracer[0] < m_coar_mag;
  
  return coarsen;
}


bool rod3d_tagger::refine_cell(const RealVect&         a_pos,
			       const Real&             a_time,
			       const Real&             a_dx,
			       const int&              a_lvl,
			       const Vector<Real>&     a_tracer,
			       const Vector<RealVect>& a_grad_tracer){
  const bool refine_tag      = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] > m_refi_curv || a_tracer[0] > m_refi_mag;
  const bool refine_geom     = (a_pos - m_rod_center).vectorLength() < 1.5*m_rod_radius;
  const bool refine_restrict = PolyGeom::dot(a_pos - m_rod_center, m_axis) < 1.2*m_rod_radius;

  return (refine_tag || refine_geom) && refine_restrict;
}
