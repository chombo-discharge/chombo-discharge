/*!
  @file   particle_containerI.H
  @brief  Implementation of particle_container.H
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _PARTICLE_CONTAINERI_H_
#define _PARTICLE_CONTAINERI_H_

#include "particle_containerI.H"
#include "NeighborIterator.H"

template <class P>
particle_container<P>::particle_container(){
  m_is_defined  = false;
  m_cell_sorted = false;
}

template <class P>
particle_container<P>::particle_container(const Vector<DisjointBoxLayout>& a_grids,
					  const Vector<ProblemDomain>&     a_domains,
					  const Vector<Real>&              a_dx,
					  const Vector<int>&               a_ref_rat,
					  const RealVect&                  a_prob_lo,
					  const int                        a_blocking_factor,
					  const int                        a_finest_level,
					  const int                        a_pvr_buffer,
					  const int                        a_halo,
					  const std::string                a_realm){

  this->define(a_grids, a_domains, a_dx, a_ref_rat, a_prob_lo, a_blocking_factor, a_finest_level, a_pvr_buffer, a_halo, a_realm);


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
				   const int                        a_pvr_buffer,
				   const int                        a_halo_buffer,
				   const std::string                a_realm){
  m_grids           = a_grids;
  m_domains         = a_domains;
  m_ref_rat         = a_ref_rat;
  m_prob_lo         = a_prob_lo;
  m_finest_level    = a_finest_level;
  m_buffer          = a_pvr_buffer;
  m_realm           = a_realm;
  m_blocking_factor = a_blocking_factor;
  m_halo_buffer     = a_halo_buffer;


  m_dx.resize(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_dx[lvl] = a_dx[lvl]*RealVect::Unit;
  }

  const int base = 0;

  // Do the define stuff. 
  this->setup_pvr(base, m_finest_level);
  this->define_halo_masks();
  this->setup_particle_data(base, m_finest_level);

  m_is_defined = true;
}

template <class P>
void particle_container<P>::setup_pvr(const int a_base, const int a_finest_level){
  m_pvr.resize(1 + a_finest_level);
  for (int lvl = a_base; lvl <= a_finest_level; lvl++){
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
void particle_container<P>::setup_particle_data(const int a_base, const int a_finest_level){
  m_particles.resize(1 + a_finest_level);
  m_halo_particles.resize(1 + a_finest_level);
  m_non_halo_particles.resize(1 + a_finest_level);
  m_cell_particles.resize(1 + a_finest_level);
  
  for (int lvl = a_base; lvl <= a_finest_level; lvl++){
    int halo_blocking_factor;
    if(lvl > 0){
      halo_blocking_factor = m_blocking_factor + m_ref_rat[lvl-1]*m_halo_buffer;
    }
    else{
      halo_blocking_factor = m_blocking_factor;
    }
    
    m_particles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
									    m_domains[lvl],
									    m_blocking_factor,
									    m_dx[lvl]*RealVect::Unit,
									    m_prob_lo));

    m_halo_particles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
										 m_domains[lvl],
										 m_blocking_factor,
										 m_dx[lvl]*RealVect::Unit,
										 m_prob_lo));

    m_non_halo_particles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
										     m_domains[lvl],
										     m_blocking_factor,
										     m_dx[lvl]*RealVect::Unit,
										     m_prob_lo));


    m_cell_particles[lvl] = RefCountedPtr<LayoutData<BinFab<P> > > (new LayoutData<BinFab<P> >(m_grids[lvl]));
  }
}

template <class P>
const std::string particle_container<P>::get_realm() const {
  return m_realm;
}

template <class P>
const AMRPVR& particle_container<P>::get_pvr() const {
  return m_pvr;
}

template <class P>
AMRParticles<P>& particle_container<P>::get_particles(){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::get_particles - particles are sorted by cell!");
  return m_particles;
}

template <class P>
const AMRParticles<P>& particle_container<P>::get_particles() const {
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::get_particles - particles are sorted by cell!");
  return m_particles;
}

template <class P>
AMRParticles<P>& particle_container<P>::get_halo_particles(){
  CH_assert(m_is_defined);
  return m_halo_particles;
}

template <class P>
const AMRParticles<P>& particle_container<P>::get_halo_particles() const {
  CH_assert(m_is_defined);
  return m_halo_particles;
}

template <class P>
AMRParticles<P>& particle_container<P>::get_non_halo_particles(){
  CH_assert(m_is_defined);
  return m_non_halo_particles;
}

template <class P>
const AMRParticles<P>& particle_container<P>::get_non_halo_particles() const {
  CH_assert(m_is_defined);
  return m_non_halo_particles;
}

template <class P>
ParticleData<P>& particle_container<P>::operator[](const int a_lvl){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::operator[](a_lvl) - particles are sorted by cell!");
  return *m_particles[a_lvl];
}

template <class P>
const ParticleData<P>& particle_container<P>::operator[](const int a_level) const {
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::operator[](a_lvl) - particles are sorted by cell!");
  return *m_particles[a_level];
}

template <class P>
AMRCellParticles<P>& particle_container<P>::get_cell_particles(){
  if(!m_cell_sorted) MayDay::Abort("particle_container<P>::get_cell_particles()- particles are not sorted by cell!");
  return m_cell_particles;
}

template <class P>
const AMRCellParticles<P>& particle_container<P>::get_cell_particles() const{
  if(!m_cell_sorted) MayDay::Abort("particle_container<P>::get_cell_particles()- particles are not sorted by cell!");
  return m_cell_particles;
}
  
template <class P>
LayoutData<BinFab<P> >& particle_container<P>::get_cell_particles(const int a_level){
  if(!m_cell_sorted) MayDay::Abort("particle_container<P>::get_cell_particles(level)- particles are not sorted by cell!");
  return *m_cell_particles[a_level];
}

template <class P>
const LayoutData<BinFab<P> >& particle_container<P>::get_cell_particles(const int a_level) const{
  if(!m_cell_sorted) MayDay::Abort("particle_container<P>::get_cell_particles(level)- particles are not sorted by cell!");
  return *m_cell_particles[a_level];
}

template <class P>
void particle_container<P>::get_cell_particles(BinFab<P>& cellParticles, const int a_lvl, const DataIndex a_dit) const {
  CH_assert(m_is_defined);
  if(!m_cell_sorted) MayDay::Abort("particle_container<P>::get_cell_particles - particles are not sorted by cell!");
  
  cellParticles.define(m_grids[a_lvl].get(a_dit),  m_dx[a_lvl], m_prob_lo);
  cellParticles.addItems((*m_particles[a_lvl])[a_dit].listItems());
}

template <class P>
void particle_container<P>::get_cell_particles_destructive(BinFab<P>& cellParticles, const int a_lvl, const DataIndex a_dit){
  CH_assert(m_is_defined);
  if(!m_cell_sorted) MayDay::Abort("particle_container<P>::get_cell_particles_destructive - particles are not sorted by cell!");
  
  cellParticles.define(m_grids[a_lvl].get(a_dit),  m_dx[a_lvl], m_prob_lo);
  cellParticles.addItemsDestructive((*m_particles[a_lvl])[a_dit].listItems());
}

template <class P>
BinFab<P>& particle_container<P>::get_cell_particles(const int a_level, const DataIndex a_dit){
  return (*m_cell_particles[a_level])[a_dit];
}

template <class P>
const BinFab<P>& particle_container<P>::get_cell_particles(const int a_level, const DataIndex a_dit) const{
  return (*m_cell_particles[a_level])[a_dit];
}

template <class P>
void particle_container<P>::sort_particles_by_cell(){
  if(!m_cell_sorted){

    for (int lvl = 0; lvl <= m_finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_grids[lvl];
      
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	BinFab<P>& cellParticles = (*m_cell_particles[lvl])[dit()];
	cellParticles.define(dbl.get(dit()), m_dx[lvl], m_prob_lo);

	cellParticles.addItemsDestructive((*m_particles[lvl])[dit()].listItems());
      }
    }
    
    m_cell_sorted = true;
  }
}

template <class P>
void particle_container<P>::sort_particles_by_patch(){
  if(m_cell_sorted){

    const int comp = 0;
    
    for (int lvl = 0; lvl <= m_finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_grids[lvl];
      
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	ListBox<P>& patchParticles = (*m_particles[lvl])[dit()];
	BinFab<P>& cellParticles   = (*m_cell_particles[lvl])[dit()];
	
	const Box box = dbl.get(dit());
	for (BoxIterator bit(box); bit.ok(); ++bit){
	  const IntVect iv = bit();
	  
	  patchParticles.addItemsDestructive(cellParticles(iv, comp));

	}

	cellParticles.clear();
      }
    }
    
    m_cell_sorted = false;
  }
}

template <class P>
void particle_container<P>::add_particles(const List<P>& a_particles){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(List<P>) - particles are sorted by cell!");
  
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
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(List<P>) - particles are sorted by cell!");

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
void particle_container<P>::add_particles(const BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit){
  CH_assert(m_is_defined);
  CH_assert(m_grids[a_lvl].get(a_dit) == a_particles.getRegion());
  if(!m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(BinFab<P>) - particles are not sorted by cell!");

  const int comp = 0;
  
  BinFab<P>& boxParticles = (*m_cell_particles[a_lvl])[a_dit];

  const Box box = m_grids[a_lvl].get(a_dit);
  for (BoxIterator bit(box); bit.ok(); ++bit){
    const IntVect iv = bit();

    List<P>& myParts       = boxParticles(iv, comp);
    const List<P>& inParts = a_particles(iv, comp);

    myParts.join(inParts);
  }
}

template <class P>
void particle_container<P>::add_particles_destructive(BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit){
  CH_assert(m_is_defined);
  CH_assert(m_grids[a_lvl].get(a_dit) == a_particles.getRegion());
  if(!m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(BinFab<P>) - particles are not sorted by cell!");

  const int comp = 0;

  BinFab<P>& boxParticles = (*m_cell_particles[a_lvl])[a_dit];

  const Box box = m_grids[a_lvl].get(a_dit);
  for (BoxIterator bit(box); bit.ok(); ++bit){
    const IntVect iv = bit();

    List<P>& myParts       = boxParticles(iv, comp);
    const List<P>& inParts = a_particles(iv, comp);

    myParts.catenate(inParts);
  }
}

template <class P>
void particle_container<P>::add_particles(const particle_container<P>& a_otherContainer){
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(particle_container<P>) - particles are sorted by cell!");
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_grids[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

      ListBox<P>& myParticles          = (*m_particles[lvl])[dit()];
      const ListBox<P>& otherParticles = a_otherContainer[lvl][dit()];

      myParticles.addItems(otherParticles.listItems());
    }
  }
}

template <class P>
void particle_container<P>::add_particles_destructive(particle_container<P>& a_otherContainer){
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(particle_container<P>) - particles are sorted by cell!");
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_grids[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

      ListBox<P>& myParticles    = (*m_particles[lvl])[dit()];
      ListBox<P>& otherParticles = a_otherContainer[lvl][dit()];

      myParticles.addItemsDestructive(otherParticles.listItems());
    }
  }
}

template <class P>
void particle_container<P>::cache_particles(particle_container<P>& a_otherContainer){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(particle_container<P>) - particles are sorted by cell!");
  
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    collectValidParticles(a_otherContainer[lvl],
			  *m_particles[lvl],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  1,
			  false, 
			  m_prob_lo);
    a_otherContainer[lvl].remapOutcast();
  }
}

template <class P>
void particle_container<P>::remap(){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(particle_container<P>) - particles are sorted by cell!");
  
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const bool has_coar = lvl > 0;

    // Rebin this level
    m_particles[lvl]->gatherOutcast();
    m_particles[lvl]->remapOutcast();

    // Collect coarser level particles onto this levels outcast list if they fit in this levels PVR
    if(has_coar){
      collectValidParticles(m_particles[lvl]->outcast(),
			    *m_particles[lvl-1],
			    m_pvr[lvl]->mask(),
			    m_dx[lvl],
			    m_ref_rat[lvl],
			    false, 
			    m_prob_lo);

      // Remap new outcasts
      m_particles[lvl]->remapOutcast();
    }


    // Transfer finer-level particles
    if(has_coar){
      collectValidParticles(m_particles[lvl-1]->outcast(),
			    *m_particles[lvl],
			    m_pvr[lvl]->mask(),
			    m_dx[lvl],
			    1,
			    true,
			    m_prob_lo);
      m_particles[lvl-1]->gatherOutcast();
      m_particles[lvl-1]->remapOutcast();
    }
  }

  // If outcasts lists are not empty (particles could have moved entirely out of the grid through the PVR), we need
  // do collect those, too
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    if(m_particles[lvl]->numOutcast() > 0){
      //      if(procID() == 0) std::cout << "particle_container<P>::remap - outcast list not empty" << std::endl;
      this->remap_lost_particles();
      break;
    }
  }
}

template <class P>
void particle_container<P>::remap_lost_particles(){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(particle_container<P>) - particles are sorted by cell!");
  
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
    m_particles[lvl]->remapOutcast();
  }
}

template <class P>
void particle_container<P>::levelRemap(){
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(particle_container<P>) - particles are sorted by cell!");
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    this->levelRemap(lvl);
  }
}

template <class P>
void particle_container<P>::levelRemap(const int a_lvl){
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(particle_container<P>) - particles are sorted by cell!");
  m_particles[a_lvl]->gatherOutcast();
  m_particles[a_lvl]->remapOutcast();
}

template <class P>
void particle_container<P>::pre_regrid(const int a_lmin){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(particle_container<P>) - particles are sorted by cell!");

  // Fill cache particles on each level
  m_cache_particles.resize(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_cache_particles[lvl].clear();
    
    for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit){
      List<P>& p = (*m_particles[lvl])[dit()].listItems();

      m_cache_particles[lvl].catenate(p);
      p.clear();
    }
  }
}

template <class P>
void particle_container<P>::regrid(const Vector<DisjointBoxLayout>& a_grids,
				   const Vector<ProblemDomain>&     a_domains,
				   const Vector<Real>&              a_dx,
				   const Vector<int>&               a_ref_rat,
				   const int                        a_lmin,
				   const int                        a_newFinestLevel){
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::regrid(...) - particles are sorted by cell!");



  // Update this stuff
  m_grids        = a_grids;
  m_domains      = a_domains;
  m_ref_rat      = a_ref_rat;
  m_finest_level = a_newFinestLevel;

  m_dx.resize(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_dx[lvl] = a_dx[lvl]*RealVect::Unit;
  }

  this->define_halo_masks();
  this->setup_pvr(a_lmin, m_finest_level);
  this->setup_particle_data(a_lmin, m_finest_level);
  
  const int oldFinestLevel = m_cache_particles.size() - 1;

  // If a level was removed, put the particles on the new finest level
  if(m_finest_level < oldFinestLevel){
    List<P>& cacheFinest = m_cache_particles[m_finest_level];
    for (int lvl = m_finest_level + 1; lvl <= oldFinestLevel; lvl++){
      List<P>& p = m_cache_particles[lvl];

      cacheFinest.catenate(p);
      p.clear();
    }
  }

  for (int lvl = 0; lvl <= a_newFinestLevel; lvl++){
    //  for (int lvl = 0; lvl <= oldFinestLevel; lvl++){
    List<P>& outcast = m_particles[lvl]->outcast();

    outcast.catenate(m_cache_particles[lvl]);
    m_cache_particles[lvl].clear();

    m_particles[lvl]->remapOutcast();
  }

  // If grids were removed the fine particle outcast list is definitely not empty. Transfer fine particles
  // to the outcast list below and remap the coarse level again. Some of these particles may be outside the
  // coar level PVR, but that is taken care of later. 
  for (int lvl = a_newFinestLevel; lvl > 0; lvl--){
    List<P>& fineParticles = m_particles[lvl]->outcast();
    List<P>& coarParticles = m_particles[lvl-1]->outcast();

    coarParticles.catenate(fineParticles);
    fineParticles.clear();

    m_particles[lvl-1]->remapOutcast();
  }

  // Now do remap so that the PVR is respected. 
  this->remap();

  m_cache_particles.resize(0);
}

template <class P>
void particle_container<P>::discard_invalid_particles(){
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::add_particles(particle_container<P>) - particles are sorted by cell!");  
  // Move particles to the outcast list on the fine level. Then discard that list. 
  for (int lvl = 1; lvl <= m_finest_level; lvl++){
    collectValidParticles(m_particles[lvl]->outcast(),
			  *m_particles[lvl-1],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  m_ref_rat[lvl-1],
			  false,
			  m_prob_lo);
    (m_particles[lvl]->outcast()).clear();
  }
}

template <class P>
size_t particle_container<P>::get_num_valid_local() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    n += m_particles[lvl]->numValidLocal();
  }

  return n;
}

template <class P>
size_t particle_container<P>::get_num_valid_global() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    n += m_particles[lvl]->numValid();
  }

  return n;
}

template <class P>
size_t particle_container<P>::get_num_outcast_local() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    n += m_particles[lvl]->numOutcastLocal();
  }

  return n;
}

template <class P>
size_t particle_container<P>::get_num_outcast_global() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    n += m_particles[lvl]->numOutcast();
  }

  return n;
}

template <class P>
size_t particle_container<P>::get_num_halo_local() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    n += m_halo_particles[lvl]->numValidLocal();
  }

  return n;
}

template <class P>
size_t particle_container<P>::get_num_halo_global() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    n += m_halo_particles[lvl]->numValid();
  }

  return n;
}

template <class P>
void particle_container<P>::define_halo_masks(){
  
  const int ncomp = 1;
  
  m_halo.resize(1 + m_finest_level);

  for (int lvl = 0; lvl < m_finest_level; lvl++){

    m_halo[lvl] = RefCountedPtr<LevelData<BaseFab<bool> > > (new LevelData<BaseFab<bool> > (m_grids[lvl], ncomp, IntVect::Zero));

    if(lvl < m_finest_level){

      const DisjointBoxLayout& gridsCoar = m_grids[lvl];
      const DisjointBoxLayout& gridsFine = m_grids[lvl+1];

      const ProblemDomain& domainCoar = m_domains[lvl];
      const ProblemDomain& domainFine = m_domains[lvl+1];

      this->define_halo_mask(*m_halo[lvl],
			     domainCoar,
			     domainFine,
			     gridsCoar,
			     gridsFine,
			     m_halo_buffer,
			     m_ref_rat[lvl]);
    }
  }
}

template <class P>
void particle_container<P>::define_halo_mask(LevelData<BaseFab<bool> >& a_coarMask,
					     const ProblemDomain&       a_domainCoar,
					     const ProblemDomain&       a_domainFine,
					     const DisjointBoxLayout&   a_gridsCoar,
					     const DisjointBoxLayout&   a_gridsFine,
					     const int                  a_buffer,
					     const int                  a_ref_rat){

  // TLDR: This function creates a GLOBAL view of the CFIVS, a_buffer cells wide on the coarse side of the refinement boundary. If this is not performant,
  //       there are probably better ways to do this. 
  
  const int comp  = 0;
  const int ncomp = 1;

  // Coarsen the fine grid and make a mask on the coarsened fine grid
  DisjointBoxLayout gridsCoFi;
  coarsen(gridsCoFi, a_gridsFine, a_ref_rat);

  IntVectSet halo;                 // Global CFIVS - should find a better way to do this.
  NeighborIterator nit(gridsCoFi); // Neighbor iterator

  // Go through the cofi grid and set the halo to true
  for (DataIterator dit(gridsCoFi); dit.ok(); ++dit){
    const Box coFiBox = gridsCoFi[dit()];

    // Make IntVect set consisting of only ghost cells (a_buffer) for each box. 
    Box grownBox = grow(coFiBox, a_buffer);
    grownBox    &= a_domainCoar;

    // Subtract non-ghosted box. 
    IntVectSet myHalo(grownBox);
    myHalo -= coFiBox;

    // Subtract non-ghosted neighboring boxes. 
    for (nit.begin(dit()); nit.ok(); ++nit){
      const Box neighborBox = gridsCoFi[nit()];
      myHalo -= neighborBox;
    }

    // Add to halo. 
    halo |= myHalo;
  }

  // Get the global cfivs
  gather_and_broadcast(halo);

  // Set the mask
  for (DataIterator dit = a_gridsCoar.dataIterator(); dit.ok(); ++dit){
    BaseFab<bool>& mask = a_coarMask[dit()];

    mask.setVal(false);

    const IntVectSet curHalo = halo & a_gridsCoar[dit()];
    for (IVSIterator ivsIt(curHalo); ivsIt.ok(); ++ivsIt){
      mask(ivsIt(), comp) = true;
    }
  }
}

template <class P>
void particle_container<P>::gather_and_broadcast(IntVectSet& a_set){
#ifdef CH_MPI
  Vector<IntVectSet> procSet;
  const int destProc = uniqueProc(SerialTask::compute);
  gather(procSet, a_set, destProc);
  a_set.makeEmpty();
  if(procID() == destProc){
    for (int i = 0; i < procSet.size(); i++){
      a_set |= procSet[i];
    }
  }
  broadcast(a_set, destProc);
#endif
}

template <class P>
void particle_container<P>::copy_halo_particles() const {
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    if(lvl < m_finest_level){
      this->copy_halo_particles(lvl);
    }
  }
}

template <class P>
void particle_container<P>::copy_halo_particles(const int a_level) const {
  
  m_halo_particles[a_level]->clear();

  const RealVect dx = m_dx[a_level];

  // 1. Copy particles from m_particles to m_halo_particles if they lie in the halo mask. 
  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit){
    const BaseFab<bool>& mask = (*m_halo[a_level])[dit()];
    const List<P>& particles  = (*m_particles[a_level])[dit()].listItems();
    List<P>& halo_particles   = (*m_halo_particles[a_level])[dit()].listItems();

    for (ListIterator<P> lit(particles); lit.ok(); ++lit){
      const RealVect x = lit().position();
      const IntVect iv(D_DECL((int)floor((x[0] - m_prob_lo[0])/dx[0]),
			      (int)floor((x[1] - m_prob_lo[1])/dx[1]),
			      (int)floor((x[2] - m_prob_lo[2])/dx[2])));

      if(mask(iv)) halo_particles.add(lit());
    }
  }
}

template <class P>
void particle_container<P>::clear_halo_particles() const {
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_halo_particles[lvl]->clear();
  }
}

template <class P>
void particle_container<P>::copy_non_halo_particles() const {
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    if(lvl < m_finest_level){
      this->copy_non_halo_particles(lvl);
    }
  }
}

template <class P>
void particle_container<P>::copy_non_halo_particles(const int a_level) const {
  
  m_non_halo_particles[a_level]->clear();

  const RealVect dx = m_dx[a_level];

  // 1. Copy particles from m_particles to m_non_halo_particles if they don't lie in the halo mask
  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit){
    const BaseFab<bool>& mask   = (*m_halo[a_level])[dit()];
    const List<P>& particles    = (*m_particles[a_level])[dit()].listItems();
    List<P>& non_halo_particles = (*m_non_halo_particles[a_level])[dit()].listItems();

    for (ListIterator<P> lit(particles); lit.ok(); ++lit){
      const RealVect x = lit().position();
      const IntVect iv(D_DECL((int)floor((x[0] - m_prob_lo[0])/dx[0]),
			      (int)floor((x[1] - m_prob_lo[1])/dx[1]),
			      (int)floor((x[2] - m_prob_lo[2])/dx[2])));

      if(!mask(iv)) non_halo_particles.add(lit());
    }
  }
}

template <class P>
void particle_container<P>::clear_non_halo_particles() const {
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_non_halo_particles[lvl]->clear();
  }
}

#endif
