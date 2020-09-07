/*!
  @file   particle_containerI.H
  @brief  Implementation of particle_container.H
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _PARTICLE_CONTAINERI_H_
#define _PARTICLE_CONTAINERI_H_

#include "particle_containerI.H"

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
					  const std::string                a_realm){

  this->define(a_grids, a_domains, a_dx, a_ref_rat, a_prob_lo, a_blocking_factor, a_finest_level, a_pvr_buffer, a_realm);
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
				   const std::string                a_realm){
  m_grids           = a_grids;
  m_domains         = a_domains;
  m_ref_rat         = a_ref_rat;
  m_prob_lo         = a_prob_lo;
  m_finest_level    = a_finest_level;
  m_buffer          = a_pvr_buffer;
  m_realm           = a_realm;
  m_blocking_factor = a_blocking_factor;


  m_dx.resize(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_dx[lvl] = a_dx[lvl]*RealVect::Unit;
  }

  const int base = 0;
  
  this->setup_pvr(base, m_finest_level);
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
  m_cell_particles.resize(1 + a_finest_level);
  for (int lvl = a_base; lvl <= a_finest_level; lvl++){
    m_particles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
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
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::levelRemap(particle_container<P>) - particles are sorted by cell!");
  m_particles[a_lvl]->gatherOutcast();
  m_particles[a_lvl]->remapOutcast();
}

template <class P>
void particle_container<P>::pre_regrid(const int a_lmin){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("particle_container<P>::pre_regrid(particle_container<P>) - particles are sorted by cell!");
  
  m_cacheParticles.resize(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_cacheParticles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
										 m_domains[lvl],
										 m_blocking_factor,
										 m_dx[lvl]*RealVect::Unit,
										 m_prob_lo));

    collectValidParticles(*m_cacheParticles[lvl],
			  *m_particles[lvl],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  1,
			  false, 
			  m_prob_lo);
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

  this->setup_pvr(a_lmin, m_finest_level);
  this->setup_particle_data(a_lmin, m_finest_level);

  const int oldFinestLevel = m_cache_particles.size() - 1;

  for (int lvl = 0; lvl <= a_newFinestLevel; lvl++){

    collectValidParticles(*m_particles[lvl],
			  *m_cacheParticles[lvl],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  1,
			  false, 
			  m_prob_lo);
    // const BoxLayout& oldGrids = m_cacheParticles[lvl]->getBoxes();

    // List<P>& outcast = m_particles[lvl]->outcast();
    
    // for (DataIterator dit(oldGrids); dit.ok(); ++dit){
    //   List<P>& cacheParticles = (*m_cacheParticles[lvl])[dit()].listItems();
    //   outcast.catenate(cacheParticles);
    //   cacheParticles.clear();
    // }

    // m_particles[lvl]->remapOutcast();
  }

  // Now do remap
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

#endif
