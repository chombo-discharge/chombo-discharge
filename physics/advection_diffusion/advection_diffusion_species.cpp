/*!
  @file   advection_diffusion_species.cpp
  @brief  Implementation of advection_diffusion_species.H
  @author Robert Marskar
  @date   March 2020
*/

#include "advection_diffusion_species.H"

#include <ParmParse.H>
#include <PolyGeom.H>

using namespace physics::advection_diffusion;

advection_diffusion_species::advection_diffusion_species(){

  ParmParse pp("advection_diffusion");

  m_charge = 0;
  m_name   = "scalar species";

  Vector<Real> v;
  pp.get   ("diffusion",      m_diffusive);
  pp.get   ("advection",      m_mobile);
  pp.get   ("blob_amplitude", m_blob_amplitude);
  pp.get   ("blob_radius",    m_blob_radius);
  pp.getarr("blob_center",    v, 0, SpaceDim); m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));
}

advection_diffusion_species::~advection_diffusion_species(){

}

Real advection_diffusion_species::initial_data(const RealVect a_pos, const Real a_time) const{
  const RealVect d = a_pos - m_blob_center;
  const Real d2 = PolyGeom::dot(d,d);
  const Real r2 = m_blob_radius*m_blob_radius;

  return m_blob_amplitude*exp(-0.5*d2*d2/(r2*r2));
}
