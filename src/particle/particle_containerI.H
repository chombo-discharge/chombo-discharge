/*!
  @file   particle_containerI.H
  @brief  Implementation of particle_container.H
  @author Robert Marskar
  @date   April 2020
*/

template <class P>
particle_container<P>::particle_container(const Vector<DisjointBoxLayout>& a_grids,
					  const Vector<ProblemDomain>&     a_domains,
					  const Vector<Real>&              a_dx,
					  const Vector<int>&               a_ref_rat,
					  const RealVect&                  a_prob_lo,
					  const int                        a_blocking_factor,
					  const int                        a_finest_level,
					  const int                        a_pvr_buffer){

  this->define(a_grids, a_domains, a_dx, a_ref_rat, a_prob_lo, a_blocking_factor, a_finest_level, a_pvr_buffer);
}

template <class P>
particle_container<P>::~particle_container(){

}

template <class P>
void particle_container<P>::define(const Vector<DisjointBoxLayout>& a_grids,
				   const Vector<ProblemDomain>&     a_domains,
				   const Vector<Real>&              a_dx,
				   const Vector<int>&               a_ref_rat,
				   const RealVect&                  a_prob_lo,
				   const int                        a_blocking_factor,
				   const int                        a_finest_level,
				   const int                        a_pvr_buffer){
  m_grids        = a_grids;
  m_domains      = a_domains;
  m_ref_rat      = a_ref_rat;
  m_prob_lo      = a_prob_lo;
  m_finest_level = a_finest_level;
  m_buffer       = a_pvr_buffer;

  m_dx.resize(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_dx[lvl] = a_dx[lvl]*RealVect::Unit;
  }

  this->setup_pvr();
  this->setup_particle_data();

  m_is_defined = true;
}

template <class P>
void particle_container<P>::setup_pvr(){
  m_pvr.resize(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_grids[lvl];
    const ProblemDomain& domain  = m_domains[lvl];

    const bool has_coar = lvl > 0;

    if(!has_coar){
      m_pvr[lvl] = RefCountedPtr<ParticleValidRegion> (new ParticleValidRegion(dbl, NULL, 1, 0));
    }
    else{
      const int ref_ratio = m_ref_rat[lvl-1];
      m_pvr[lvl] = RefCountedPtr<ParticleValidRegion> (new ParticleValidRegion(dbl, m_pvr[lvl-1]->mask(), ref_ratio, m_buffer));
    }
  }
}

template <class P>
void particle_container<P>::setup_particle_data(){
  m_particles.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_particles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
									    m_domains[lvl],
									    m_blocking_factor,
									    m_dx[lvl]*RealVect::Unit,
									    m_prob_lo));
  }
}

template <class P>
AMRParticles<P>& particle_container<P>::get_particles(){
  return m_particles;
}

template <class P>
ParticleData<P>& particle_container<P>::operator[](const int a_lvl){
  return *m_particles[a_lvl];
}

template <class P>
void particle_container<P>::add_particles(const List<P>& a_particles){

  // Put the initial particles on the coarsest grid level
  List<P>& outcastBase = m_particles[0]->outcast();
  outcastBase.join(a_particles);
  m_particles[0]->remapOutcast(); 

  // Move particles to finer levels if they belong there. This piece of code moves particles from lvl-1
  // and onto the outcast list on level lvl. Then, we remap the outcast list
  for (int lvl = 1; lvl <= m_finest_level; lvl++){
    collectValidParticles(m_particles[lvl]->outcast(),
			  *m_particles[lvl-1],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  m_ref_rat[lvl-1],
			  false,
			  m_prob_lo);
    m_particles[lvl]->remapOutcast();
  }
}

template <class P>
void particle_container<P>::add_particles_destructive(List<P>& a_particles){
  MayDay::Abort("particle_container::add_particles_destructive - not implemented");

  // Put the initial particles on the coarsest grid level
  List<P>& outcastBase = m_particles[0]->outcast();
  outcastBase.catenate(a_particles);
  m_particles[0]->remapOutcast(); 

  // Move particles to finer levels if they belong there. This piece of code moves particles from lvl-1
  // and onto the outcast list on level lvl. Then, we remap the outcast list
  for (int lvl = 1; lvl <= m_finest_level; lvl++){
    collectValidParticles(m_particles[lvl]->outcast(),
			  *m_particles[lvl-1],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  m_ref_rat[lvl-1],
			  false,
			  m_prob_lo);
    m_particles[lvl]->remapOutcast();
  }
}

template <class P>
void particle_container<P>::cache_particles(particle_container<P>& a_otherContainer){
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    collectValidParticles(a_otherContainer[lvl],
			  *m_particles[lvl],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  1,
			  false, 
			  m_prob_lo);
    a_otherContainer[lvl]->remapOutcast();
  }
}

template <class P>
void particle_container<P>::remap(){

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const bool has_coar = lvl > 0;

    // Collect coarser level particles onto this levels outcast list if they fit in this levels PVR
    if(has_coar){
      collectValidParticles(m_particles[lvl]->outcast(),
			    *m_particles[lvl-1],
			    m_pvr[lvl]->mask(),
			    m_dx[lvl],
			    m_ref_rat[lvl],
			    false, 
			    m_prob_lo);
    }

    // Rebin this level
    m_particles[lvl]->gatherOutcast();
    m_particles[lvl]->remapOutcast();

    // Transfer finer-level particles
    if(has_coar){
      collectValidParticles(m_particles[lvl-1]->outcast(),
			    *m_particles[lvl],
			    m_pvr[lvl]->mask(),
			    m_dx[lvl],
			    1,
			    true,
			    m_prob_lo);
      m_particles[lvl-1]->remapOutcast();
    }
  }

  // If outcasts lists are not empty (particles could have moved entirely out of the grid through the PVR), we need
  // do collect those, too
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    if(m_particles->numOutcast() > 0){
      this->remap_lost_particles();
      break;
    }
  }
}

template <class P>
void particle_container<P>::remap_lost_particles(){
  // TLDR: This routine is called because there is a chance that particles might hop across more than one refinement boundary,
  //       and the remap_amr_particles only performs two-level operations. This piece of code is responsible for gathering
  //       "lost particles" which were not successfully remapped through two-level operations

  // Gather all "lost" outcasts on the coarsest level
  List<P>& coarsest_outcast = m_particles[0]->outcast();
  for (int lvl = 1; lvl <= m_finest_level; lvl++){
    coarsest_outcast.catenate(m_particles[lvl]->outcast());
    m_particles[lvl]->outcast().clear();
  }

  // Remap the coarsest level and discard particles that don't fit in the coarsest level (i.e. ones that have left the domain)
  m_particles[0]->remapOutcast();
  m_particles[0]->outcast().clear();

  // Collect the particles 
  for (int lvl = 1; lvl <= m_finest_level; lvl++){
    collectValidParticles(m_particles[lvl]->outcast(),
			  *m_particles[lvl-1],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  m_ref_rat[lvl-1],
			  false, 
			  m_prob_lo);
  }
}
