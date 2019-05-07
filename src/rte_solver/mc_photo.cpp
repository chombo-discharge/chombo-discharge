/*!
  @file   mc_photo.cpp
  @brief  Implementation of mc_photo.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "mc_photo.H"
#include "data_ops.H"

#include <BoxIterator.H>

mc_photo::mc_photo(){
  this->set_verbosity(-1);
  this->set_stationary(false);
}

mc_photo::~mc_photo(){

}

bool mc_photo::advance(const Real a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source, const bool a_zerophi){
  data_ops::set_value(a_state, 0.0);

  const RealVect origin  = m_physdom->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();
  const int boxsize      = m_amr->get_max_box_size();
  if(boxsize != m_amr->get_blocking_factor()){
    MayDay::Abort("mc_photo::advance - only constant box sizes supported for particle methods");
  }
  
  Vector<ParticleData<Particle>* > n_particles(1 + finest_level);
  Vector<ParticleValidRegion* > n_pvr(1 + finest_level);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& dom = m_amr->get_domains()[lvl];
    const Real vol = pow(dx, SpaceDim);
    const int ref_ratio = m_amr->get_ref_rat()[lvl];

    const bool has_coar = lvl > 0;

    n_particles[lvl] = new ParticleData<Particle>(dbl, dom, boxsize, dx*RealVect::Unit, origin);
    if(!has_coar){
      n_pvr[lvl] = new ParticleValidRegion(dbl, NULL, 1, 0);
    }
    if(has_coar) { // Get particles from coarser region and put onto this one
      collectValidParticles(n_particles[lvl]->outcast(),
      			    *n_particles[lvl-1],
      			    m_pvr[lvl]->mask(),
      			    dx*RealVect::Unit,
      			    ref_ratio,
			    false,
			    origin);
      m_particles[lvl]->outcast().clear();
    }

    // Create new particles on this level using fine-resolution data
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box = dbl.get(dit());
      const EBISBox& ebisbox = (*a_state[lvl])[dit()].getEBISBox();
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = IntVectSet(box);


      FArrayBox& source = (*a_source[lvl])[dit()].getFArrayBox();
      FArrayBox& state  = (*a_state[lvl])[dit()].getFArrayBox();


      List<Particle> particles;
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const IntVect iv = vof.gridIndex();
	const RealVect pos = origin + ((RealVect)iv + 0.5*RealVect::Unit)*dx;

	Particle part(source(iv,0)*vol*a_dt, pos);
	particles.append(part);
      }
      (*n_particles[lvl])[dit()].addItemsDestructive(particles);

      // Interp
      MeshInterp interp(box, dx*RealVect::Unit, origin);
      InterpType type = InterpType::CIC;
      interp.deposit((*n_particles[lvl])[dit()].listItems(), state, type);
    }

    n_particles[lvl]->gatherOutcast();
    n_particles[lvl]->remapOutcast();
  }

  //   // Deposit on particle on each patch on the coarsest level
  // const RealVect origin = m_physdom->get_prob_lo();
  // const Real dx = m_amr->get_dx()[0];
  // const DisjointBoxLayout& dbl = m_amr->get_grids()[0];
  // const ProblemDomain& dom = m_amr->get_domains()[0];
  // const int boxsize = m_amr->get_max_box_size();

  // ParticleData<Particle> m_particles(dbl, dom, boxsize, dx*RealVect::Unit, origin);

  // // Create initial particles
  // for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
  //   const Box& box = dbl.get(dit());

  //   List<Particle> particles;

  //   FArrayBox& state  = (*a_state[0])[dit()].getFArrayBox();
  //   FArrayBox& source = (*a_source[0])[dit()].getFArrayBox();
    
  //   for (BoxIterator bit(box); bit.ok(); ++bit){
  //     const IntVect iv = bit();
  //     const RealVect pos = origin + ((RealVect)iv+0.5)*dx*RealVect::Unit;

  //     Particle part(source(iv,0)/1.E20, pos);
  //     particles.append(part);
  //   }
  //   m_particles[dit()].addItemsDestructive(particles);
  // }

  // // Move all the particles a little bit
  // for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
  //   const Box& box = dbl.get(dit());

  //   List<Particle>& particles = m_particles[dit()].listItems();

  //   for (ListIterator<Particle> lit(particles); lit.ok(); ++lit){

  //     Particle& p = lit();

  //     p.position() += -1.E-3*RealVect(0.,1.);
  //   }
  // }

  // m_particles.gatherOutcast();
  // m_particles.remapOutcast();


  // // Interpolate to mesh
  // for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
  //   const Box box = dbl.get(dit());
    
  //   FArrayBox& state  = (*a_state[0])[dit()].getFArrayBox();

  //   List<Particle>& particles = m_particles[dit()].listItems();
    
  //   // Interp
  //   MeshInterp interp(box, dx*RealVect::Unit, origin);
  //   InterpType type = InterpType::CIC;
  //   interp.deposit(particles, state, type);
  // }
}
  
void mc_photo::allocate_internals(){
  CH_TIME("mc_photo::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  const int buffer = 1;
  const int ncomp  = 1;
  m_amr->allocate(m_state,  m_phase, ncomp); // This is the deposited 
  m_amr->allocate(m_source, m_phase, ncomp);
  
  m_amr->allocate(m_particles);
  m_amr->allocate(m_pvr, buffer);
}
  
void mc_photo::cache_state(){
  CH_TIME("mc_photo::cache_state");
  if(m_verbosity > 5){
    pout() << m_name + "::cache_state" << endl;
  }
}

void mc_photo::deallocate_internals(){
  m_amr->deallocate(m_state);
  m_amr->deallocate(m_source);
}

void mc_photo::regrid(const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("mc_photo::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  this->allocate_internals();
}

void mc_photo::compute_boundary_flux(EBAMRIVData& a_ebflux, const EBAMRCellData& a_state){
  CH_TIME("mc_photo::compute_boundary_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_boundary_flux" << endl;
  }
  data_ops::set_value(a_ebflux, 0.0);
}


void mc_photo::compute_domain_flux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_state){
  CH_TIME("mc_photo::compute_domain_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_domain_flux" << endl;
  }
  data_ops::set_value(a_domainflux, 0.0);
}

void mc_photo::compute_flux(EBAMRCellData& a_flux, const EBAMRCellData& a_state){
  MayDay::Abort("mc_photo::compute_flux - I don't think this should ever be called.");
}

void mc_photo::compute_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_state){
  MayDay::Abort("mc_photo::compute_density - I don't think this should ever be called.");
}


void mc_photo::write_plot_file(){
  CH_TIME("mc_photo::write_plot_file");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_file" << endl;
  }
}

int mc_photo::query_ghost() const {
  return 3;
}
