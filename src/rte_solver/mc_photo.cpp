/*!
  @file   mc_photo.cpp
  @brief  Implementation of mc_photo.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "mc_photo.H"
#include "data_ops.H"

#include <BoxIterator.H>
#include <Particle.H>
#include <ParticleData.H>
#include <BinItem.H>
#include <MeshInterp.H>

mc_photo::mc_photo(){

}

mc_photo::~mc_photo(){

}

bool mc_photo::advance(const Real a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source, const bool a_zerophi){
  data_ops::set_value(a_state, 0.0);

  // Deposit on particle on each patch on the coarsest level
  const RealVect origin = m_physdom->get_prob_lo();
  const Real dx = m_amr->get_dx()[0];
  const DisjointBoxLayout& dbl = m_amr->get_grids()[0];
  const ProblemDomain& dom = m_amr->get_domains()[0];
  const int boxsize = m_amr->get_max_box_size();

  ParticleData<Particle> m_particles(dbl, dom, boxsize, dx*RealVect::Unit, origin);

  // Create initial particles
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box = dbl.get(dit());

    List<Particle> particles;

    FArrayBox& state  = (*a_state[0])[dit()].getFArrayBox();
    FArrayBox& source = (*a_source[0])[dit()].getFArrayBox();
    
    for (BoxIterator bit(box); bit.ok(); ++bit){
      const IntVect iv = bit();
      const RealVect pos = origin + ((RealVect)iv+0.5)*dx*RealVect::Unit;

      Particle part(source(iv,0)/1.E20, pos);
      particles.append(part);
    }
    m_particles[dit()].addItemsDestructive(particles);
  }

  // Move all the particles a little bit
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box = dbl.get(dit());

    List<Particle>& particles = m_particles[dit()].listItems();

    for (ListIterator<Particle> lit(particles); lit.ok(); ++lit){

      Particle& p = lit();

      p.position() += -1.E-3*RealVect(0.,1.);
    }
  }

  m_particles.gatherOutcast();
  m_particles.remapOutcast();


  // Interpolate to mesh
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());
    
    FArrayBox& state  = (*a_state[0])[dit()].getFArrayBox();

    List<Particle>& particles = m_particles[dit()].listItems();
    
    // Interp
    MeshInterp interp(box, dx*RealVect::Unit, origin);
    InterpType type = InterpType::CIC;
    interp.deposit(particles, state, type);
  }

  
}
  
void mc_photo::allocate_internals(){
  const int ncomp = 1;
  m_amr->allocate(m_state,  m_phase, ncomp); // This is the deposited 
  m_amr->allocate(m_source, m_phase, ncomp);
}
  
void mc_photo::cache_state(){

}

void mc_photo::deallocate_internals(){
  m_amr->deallocate(m_state);
  m_amr->deallocate(m_source);
}

void mc_photo::regrid(const int a_old_finest_level, const int a_new_finest_level){

}

void mc_photo::compute_boundary_flux(EBAMRIVData& a_ebflux, const EBAMRCellData& a_state){
  data_ops::set_value(a_ebflux, 0.0);
}


void mc_photo::compute_domain_flux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_state){
  data_ops::set_value(a_domainflux, 0.0);
}

void mc_photo::compute_flux(EBAMRCellData& a_flux, const EBAMRCellData& a_state){
  MayDay::Abort("mc_photo::compute_flux - I don't think this should ever be called.");
}

void mc_photo::compute_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_state){
  MayDay::Abort("mc_photo::compute_density - I don't think this should ever be called.");
}


void mc_photo::write_plot_file(){

}

int mc_photo::query_ghost() const {
  return 3;
}
