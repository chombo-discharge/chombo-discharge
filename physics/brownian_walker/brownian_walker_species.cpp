/*!
  @file   brownian_walker_species.cpp
  @brief  Implementation of brownian_walker_species.H
  @author Robert Marskar
  @date   March 2020
*/

#include "brownian_walker_species.H"

#include <ParmParse.H>

using namespace physics::brownian_walker;

brownian_walker_species::brownian_walker_species(){
  m_name   = "scalar species";
  m_charge = 0;

  ParmParse pp("brownian_walker");
  Vector<Real> v;
  
  pp.get   ("diffusion",      m_diffusive);
  pp.get   ("advection",      m_mobile);
  pp.get   ("blob_amplitude", m_blob_amplitude);
  pp.get   ("blob_radius",    m_blob_radius);
  pp.getarr("blob_center",    v, 0, SpaceDim); m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.get   ("num_particles",  m_num_particles);
}

brownian_walker_species::~brownian_walker_species(){

}
