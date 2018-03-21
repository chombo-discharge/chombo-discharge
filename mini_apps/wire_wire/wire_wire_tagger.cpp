/*!
  @file   wire_wire_tagger.cpp
  @brief  Implementation wire_wire_tagger.cpp
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "wire_wire_tagger.H"

#include <ParmParse.H>
#include <PolyGeom.H>

wire_wire_tagger::wire_wire_tagger(){
  m_num_tracers = 1;

  m_coar_curv  = 0.1;
  m_coar_mag   = 0.1;
  m_refi_curv  = 0.2;
  m_refi_mag   = 0.75;
  m_radius_fac = 1.5;
  m_geom_refine_depth = 0;

  // Get input parameters if we have them
  {
    ParmParse pp("wire_wire_tagger");

    pp.query("coarsen_curvature", m_coar_curv);
    pp.query("coarsen_magnitude", m_coar_mag);
    pp.query("refine_curvature",  m_refi_curv);
    pp.query("refine_magnitude",  m_refi_mag);
    pp.query("refine_radius",     m_radius_fac);
    pp.query("refine_geom_depth", m_geom_refine_depth);
  }

  { // Geometrical inputs
    Real r1, r2;
    Vector<Real> v1, v2, v3, v4;
    ParmParse pp("wire_wire_geometry");
    pp.get("wire1_radius", r1);
    pp.get("wire2_radius", r2);
    pp.getarr("wire1_center1", v1, 0, SpaceDim);
    pp.getarr("wire2_center1", v2, 0, SpaceDim);
#if CH_SPACEDIM == 3
    pp.getarr("wire1_center2", v3, 0, SpaceDim);
    pp.getarr("wire2_center2", v4, 0, SpaceDim);
#endif

    const RealVect c1 = RealVect(D_DECL(v1[0], v1[1], v1[2]));
    const RealVect c2 = RealVect(D_DECL(v2[0], v2[1], v2[2]));
#if CH_SPACEDIM==3
    const RealVect c3 = RealVect(D_DECL(v3[0], v3[1], v3[2]));
    const RealVect c4 = RealVect(D_DECL(v4[0], v4[1], v4[2]));
#endif

    const Real r = Max(r1, r2);
#if CH_SPACEDIM == 2

    if(c1 <= c2){
      m_lo = c1 - m_radius_fac*r*RealVect::Unit;
      m_hi = c2 + m_radius_fac*r*RealVect::Unit;
    }
    else {
      m_lo = c2 - m_radius_fac*r*RealVect::Unit;
      m_hi = c1 + m_radius_fac*r*RealVect::Unit;
    }
#elif CH_SPACEDIM==3
    Vector<RealVect> c(4);
    c[0] = c1;
    c[1] = c2;
    c[2] = c3;
    c[3] = c4;

    m_lo =  1.E99*RealVect::Unit;
    m_hi = -1.E99*RealVect::Unit;
    for (int i = 0; i < 4; i++){
      if(c[i] < m_lo) m_lo = c[i];
      if(c[i] > m_hi) m_hi = c[i];
    }

    m_lo = m_lo - m_radius_fac*r*RealVect::Unit;
    m_hi = m_hi - m_radius_fac*r*RealVect::Unit;
#endif
  }
  
}


wire_wire_tagger::~wire_wire_tagger(){

}


Vector<Real> wire_wire_tagger::tracer(const RealVect&         a_pos,
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


bool wire_wire_tagger::coarsen_cell(const RealVect&         a_pos,
				    const Real&             a_time,
				    const Real&             a_dx,
				    const int&              a_lvl,
				    const Vector<Real>&     a_tracer,
				    const Vector<RealVect>& a_grad_tracer){

  const bool coarsen = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] < m_coar_curv && a_tracer[0] < m_coar_mag;
  
  return coarsen;
}


bool wire_wire_tagger::refine_cell(const RealVect&         a_pos,
				   const Real&             a_time,
				   const Real&             a_dx,
				   const int&              a_lvl,
				   const Vector<Real>&     a_tracer,
				   const Vector<RealVect>& a_grad_tracer){

  const bool refine_curv     = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] > m_refi_curv && a_tracer[0] > 1.E-3;
  const bool refine_magn     = a_tracer[0] > m_refi_mag;
  const bool refine_geom     = (a_pos > m_lo) && (a_pos < m_hi) && (a_lvl < m_geom_refine_depth);
  return refine_curv || refine_magn || refine_geom;
}
