/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ParticleContainerImplem.H
  @brief  Implementation of CD_ParticleContainer.H
  @author Robert Marskar
*/

#ifndef CD_ParticleContainerImplem_H
#define CD_ParticleContainerImplem_H

// Chombo includes
#include <NeighborIterator.H>

// Our includes
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>
  
template <class P>
ParticleContainer<P>::ParticleContainer(){
  m_is_defined  = false;
  m_cell_sorted = false;
}

template <class P>
ParticleContainer<P>::ParticleContainer(const Vector<DisjointBoxLayout>& a_grids,
					const Vector<ProblemDomain>&     a_domains,
					const Vector<Real>&              a_dx,
					const Vector<int>&               a_ref_rat,
					const RealVect&                  a_probLo,
					const int                        a_blocking_factor,
					const int                        a_finestLevel,
					const int                        a_pvr_buffer,
					const std::string                a_realm){

  this->define(a_grids, a_domains, a_dx, a_ref_rat, a_probLo, a_blocking_factor, a_finestLevel, a_pvr_buffer, a_realm);


}

template <class P>
ParticleContainer<P>::~ParticleContainer(){

}

template <class P>
void ParticleContainer<P>::define(const Vector<DisjointBoxLayout>& a_grids,
				  const Vector<ProblemDomain>&     a_domains,
				  const Vector<Real>&              a_dx,
				  const Vector<int>&               a_ref_rat,
				  const RealVect&                  a_probLo,
				  const int                        a_blocking_factor,
				  const int                        a_finestLevel,
				  const int                        a_pvr_buffer,
				  const std::string                a_realm){
  m_grids           = a_grids;
  m_domains         = a_domains;
  m_ref_rat         = a_ref_rat;
  m_probLo         = a_probLo;
  m_finestLevel    = a_finestLevel;
  m_buffer          = a_pvr_buffer;
  m_realm           = a_realm;
  m_blockingFactor = a_blocking_factor;

  m_dx.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    m_dx[lvl] = a_dx[lvl]*RealVect::Unit;
  }

  const int base = 0;

  // Do the define stuff. 
  this->setupPVR(base, m_finestLevel);
  this->setupParticleData(base, m_finestLevel);

  m_is_defined  = true;
  m_cell_sorted = false;
}

template <class P>
void ParticleContainer<P>::setupPVR(const int a_base, const int a_finestLevel){
  m_pvr.resize(1 + a_finestLevel);
  for (int lvl = a_base; lvl <= a_finestLevel; lvl++){
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
void ParticleContainer<P>::setupParticleData(const int a_base, const int a_finestLevel){
  m_particles.resize(1 + a_finestLevel);
  m_halo_particles.resize(1 + a_finestLevel);
  m_non_halo_particles.resize(1 + a_finestLevel);
  m_cell_particles.resize(1 + a_finestLevel);
  
  for (int lvl = a_base; lvl <= a_finestLevel; lvl++){
    m_particles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
									    m_domains[lvl],
									    m_blockingFactor,
									    m_dx[lvl]*RealVect::Unit,
									    m_probLo));

    m_halo_particles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
										 m_domains[lvl],
										 m_blockingFactor,
										 m_dx[lvl]*RealVect::Unit,
										 m_probLo));

    m_non_halo_particles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
										     m_domains[lvl],
										     m_blockingFactor,
										     m_dx[lvl]*RealVect::Unit,
										     m_probLo));


    m_cell_particles[lvl] = RefCountedPtr<LayoutData<BinFab<P> > > (new LayoutData<BinFab<P> >(m_grids[lvl]));
  }
}

template <class P>
int ParticleContainer<P>::getFinestLevel() const {
  return m_finestLevel;
}

template <class P>
const std::string ParticleContainer<P>::getRealm() const {
  return m_realm;
}

template <class P>
const AMRPVR& ParticleContainer<P>::getPVR() const {
  return m_pvr;
}

template <class P>
const Vector<DisjointBoxLayout>& ParticleContainer<P>::getGrids() const {
  return m_grids;
}

template <class P>
AMRParticles<P>& ParticleContainer<P>::getParticles(){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::getParticles - particles are sorted by cell!");
  return m_particles;
}

template <class P>
const AMRParticles<P>& ParticleContainer<P>::getParticles() const {
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::getParticles - particles are sorted by cell!");
  return m_particles;
}

template <class P>
AMRParticles<P>& ParticleContainer<P>::getMaskParticles(){
  CH_assert(m_is_defined);
  return m_halo_particles;
}

template <class P>
const AMRParticles<P>& ParticleContainer<P>::getMaskParticles() const {
  CH_assert(m_is_defined);
  return m_halo_particles;
}

template <class P>
AMRParticles<P>& ParticleContainer<P>::getNonMaskParticles(){
  CH_assert(m_is_defined);
  return m_non_halo_particles;
}

template <class P>
const AMRParticles<P>& ParticleContainer<P>::getNonMaskParticles() const {
  CH_assert(m_is_defined);
  return m_non_halo_particles;
}

template <class P>
Vector<List<P> >& ParticleContainer<P>::getCacheParticles() {
  return m_cacheParticles;
}

template <class P>
const Vector<List<P> >& ParticleContainer<P>::getCacheParticles() const {
  return m_cacheParticles;
}

template <class P>
ParticleData<P>& ParticleContainer<P>::operator[](const int a_lvl){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::operator[](a_lvl) - particles are sorted by cell!");
  return *m_particles[a_lvl];
}

template <class P>
const ParticleData<P>& ParticleContainer<P>::operator[](const int a_level) const {
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::operator[](a_lvl) - particles are sorted by cell!");
  return *m_particles[a_level];
}

template <class P>
AMRCellParticles<P>& ParticleContainer<P>::getCellParticles(){
  if(!m_cell_sorted) MayDay::Abort("ParticleContainer<P>::getCellParticles()- particles are not sorted by cell!");
  return m_cell_particles;
}

template <class P>
const AMRCellParticles<P>& ParticleContainer<P>::getCellParticles() const{
  if(!m_cell_sorted) MayDay::Abort("ParticleContainer<P>::getCellParticles()- particles are not sorted by cell!");
  return m_cell_particles;
}
  
template <class P>
LayoutData<BinFab<P> >& ParticleContainer<P>::getCellParticles(const int a_level){
  if(!m_cell_sorted) MayDay::Abort("ParticleContainer<P>::getCellParticles(level)- particles are not sorted by cell!");
  return *m_cell_particles[a_level];
}

template <class P>
const LayoutData<BinFab<P> >& ParticleContainer<P>::getCellParticles(const int a_level) const{
  if(!m_cell_sorted) MayDay::Abort("ParticleContainer<P>::getCellParticles(level)- particles are not sorted by cell!");
  return *m_cell_particles[a_level];
}

template <class P>
void ParticleContainer<P>::getCellParticles(BinFab<P>& cellParticles, const int a_lvl, const DataIndex a_dit) const {
  CH_assert(m_is_defined);
  if(!m_cell_sorted) MayDay::Abort("ParticleContainer<P>::getCellParticles - particles are not sorted by cell!");
  
  cellParticles.define(m_grids[a_lvl].get(a_dit),  m_dx[a_lvl], m_probLo);
  cellParticles.addItems((*m_particles[a_lvl])[a_dit].listItems());
}

template <class P>
void ParticleContainer<P>::getCellParticlesDestructive(BinFab<P>& cellParticles, const int a_lvl, const DataIndex a_dit){
  CH_assert(m_is_defined);
  if(!m_cell_sorted) MayDay::Abort("ParticleContainer<P>::getCellParticlesDestructive - particles are not sorted by cell!");
  
  cellParticles.define(m_grids[a_lvl].get(a_dit),  m_dx[a_lvl], m_probLo);
  cellParticles.addItemsDestructive((*m_particles[a_lvl])[a_dit].listItems());
}

template <class P>
BinFab<P>& ParticleContainer<P>::getCellParticles(const int a_level, const DataIndex a_dit){
  return (*m_cell_particles[a_level])[a_dit];
}

template <class P>
const BinFab<P>& ParticleContainer<P>::getCellParticles(const int a_level, const DataIndex a_dit) const{
  return (*m_cell_particles[a_level])[a_dit];
}

template <class P>
void ParticleContainer<P>::sortParticlesByCell(){
  if(!m_cell_sorted){

    for (int lvl = 0; lvl <= m_finestLevel; lvl++){
      const DisjointBoxLayout& dbl = m_grids[lvl];
      
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	BinFab<P>& cellParticles = (*m_cell_particles[lvl])[dit()];
	cellParticles.define(dbl.get(dit()), m_dx[lvl], m_probLo);

	cellParticles.addItemsDestructive((*m_particles[lvl])[dit()].listItems());
      }
    }
    
    m_cell_sorted = true;
  }
}

template <class P>
void ParticleContainer<P>::sortParticlesByPatch(){
  if(m_cell_sorted){

    const int comp = 0;
    
    for (int lvl = 0; lvl <= m_finestLevel; lvl++){
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
void ParticleContainer<P>::addParticles(const List<P>& a_particles){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(List<P>) - particles are sorted by cell!");
  
  // Put the initial particles on the coarsest grid level
  List<P>& outcastBase = m_particles[0]->outcast();
  outcastBase.join(a_particles);
  m_particles[0]->remapOutcast(); 

  // Move particles to finer levels if they belong there. This piece of code moves particles from lvl-1
  // and onto the outcast list on level lvl. Then, we remap the outcast list
  for (int lvl = 1; lvl <= m_finestLevel; lvl++){
    collectValidParticles(m_particles[lvl]->outcast(),
			  *m_particles[lvl-1],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  m_ref_rat[lvl-1],
			  false,
			  m_probLo);
    m_particles[lvl]->remapOutcast();
  }
}

template <class P>
void ParticleContainer<P>::addParticlesDestructive(List<P>& a_particles){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(List<P>) - particles are sorted by cell!");

  // Put the initial particles on the coarsest grid level
  List<P>& outcastBase = m_particles[0]->outcast();
  outcastBase.catenate(a_particles);
  m_particles[0]->remapOutcast(); 

  // Move particles to finer levels if they belong there. This piece of code moves particles from lvl-1
  // and onto the outcast list on level lvl. Then, we remap the outcast list
  for (int lvl = 1; lvl <= m_finestLevel; lvl++){
    collectValidParticles(m_particles[lvl]->outcast(),
			  *m_particles[lvl-1],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  m_ref_rat[lvl-1],
			  false,
			  m_probLo);
    m_particles[lvl]->remapOutcast();
  }
}

template <class P>
void ParticleContainer<P>::addParticles(const BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit){
  CH_assert(m_is_defined);
  CH_assert(m_grids[a_lvl].get(a_dit) == a_particles.getRegion());
  if(!m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(BinFab<P>) - particles are not sorted by cell!");

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
void ParticleContainer<P>::addParticlesDestructive(BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit){
  CH_assert(m_is_defined);
  CH_assert(m_grids[a_lvl].get(a_dit) == a_particles.getRegion());
  if(!m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(BinFab<P>) - particles are not sorted by cell!");

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
void ParticleContainer<P>::addParticles(const ParticleContainer<P>& a_otherContainer){
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");

  if(m_realm == a_otherContainer.getRealm()){
    for (int lvl = 0; lvl <= m_finestLevel; lvl++){
      const DisjointBoxLayout& dbl = m_grids[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	ListBox<P>& myParticles          = (*m_particles[lvl])[dit()];
	const ListBox<P>& otherParticles = a_otherContainer[lvl][dit()];

	myParticles.addItems(otherParticles.listItems());
      }
    }
  }
  else{
    for (int lvl = 0; lvl <= m_finestLevel; lvl++){
      List<P>& outcast = m_particles[lvl]->outcast();
      outcast.clear();

      const DisjointBoxLayout& dbl = a_otherContainer.getGrids()[lvl];
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const List<P>& otherParticles = a_otherContainer[lvl][dit()].listItems();

	outcast.join(otherParticles);
      }
    }

    this->remap();
  }
}

template <class P>
void ParticleContainer<P>::addParticlesDestructive(ParticleContainer<P>& a_otherContainer){
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");

  if(m_realm == a_otherContainer.getRealm()){
    for (int lvl = 0; lvl <= m_finestLevel; lvl++){
      const DisjointBoxLayout& dbl = m_grids[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	ListBox<P>& myParticles    = (*m_particles[lvl])[dit()];
	ListBox<P>& otherParticles = a_otherContainer[lvl][dit()];

	myParticles.addItemsDestructive(otherParticles.listItems());
      }
    }
  }
  else{
    for (int lvl = 0; lvl <= m_finestLevel; lvl++){
      List<P>& outcast = m_particles[lvl]->outcast();
      outcast.clear();

      const DisjointBoxLayout& dbl = a_otherContainer.getGrids()[lvl];
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	List<P>& otherParticles = a_otherContainer[lvl][dit()].listItems();

	outcast.catenate(otherParticles);

	otherParticles.clear();
      }
    }

    this->remap();
  }
}

template <class P>
void ParticleContainer<P>::cacheParticles(ParticleContainer<P>& a_otherContainer){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    collectValidParticles(a_otherContainer[lvl],
			  *m_particles[lvl],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  1,
			  false, 
			  m_probLo);
    a_otherContainer[lvl].remapOutcast();
  }
}

template <class P>
void ParticleContainer<P>::remap(){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
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
			    m_probLo);

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
			    m_probLo);
      m_particles[lvl-1]->gatherOutcast();
      m_particles[lvl-1]->remapOutcast();
    }
  }

  // If outcasts lists are not empty (particles could have moved entirely out of the grid through the PVR), we need
  // do collect those, too
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    if(m_particles[lvl]->numOutcast() > 0){
      //      if(procID() == 0) std::cout << "ParticleContainer<P>::remap - outcast list not empty" << std::endl;
      this->remapLostParticles();
      break;
    }
  }
}

template <class P>
void ParticleContainer<P>::remapLostParticles(){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  
  // TLDR: This routine is called because there is a chance that particles might hop across more than one refinement boundary,
  //       and the remap_amr_particles only performs two-level operations. This piece of code is responsible for gathering
  //       "lost particles" which were not successfully remapped through two-level operations

  // Gather all "lost" outcasts on the coarsest level
  List<P>& coarsest_outcast = m_particles[0]->outcast();
  for (int lvl = 1; lvl <= m_finestLevel; lvl++){
    coarsest_outcast.catenate(m_particles[lvl]->outcast());
    m_particles[lvl]->outcast().clear();
  }

  // Remap the coarsest level and discard particles that don't fit in the coarsest level (i.e. ones that have left the domain)
  m_particles[0]->remapOutcast();
  m_particles[0]->outcast().clear();

  // Collect the particles 
  for (int lvl = 1; lvl <= m_finestLevel; lvl++){
    collectValidParticles(m_particles[lvl]->outcast(),
			  *m_particles[lvl-1],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  m_ref_rat[lvl-1],
			  false, 
			  m_probLo);
    m_particles[lvl]->remapOutcast();
  }
}

template <class P>
void ParticleContainer<P>::levelRemap(){
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    this->levelRemap(lvl);
  }
}

template <class P>
void ParticleContainer<P>::levelRemap(const int a_lvl){
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  m_particles[a_lvl]->gatherOutcast();
  m_particles[a_lvl]->remapOutcast();
}

template <class P>
void ParticleContainer<P>::preRegrid(const int a_lmin){
  CH_assert(m_is_defined);
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");

  // Fill cache particles on each level
  m_cacheParticles.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    m_cacheParticles[lvl].clear();
    
    for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit){
      List<P>& p = (*m_particles[lvl])[dit()].listItems();

      m_cacheParticles[lvl].catenate(p);
      p.clear();
    }
  }
}

template <class P>
void ParticleContainer<P>::regrid(const Vector<DisjointBoxLayout>& a_grids,
				  const Vector<ProblemDomain>&     a_domains,
				  const Vector<Real>&              a_dx,
				  const Vector<int>&               a_ref_rat,
				  const int                        a_lmin,
				  const int                        a_newFinestLevel){
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::regrid(...) - particles are sorted by cell!");

  // Update this stuff
  m_grids        = a_grids;
  m_domains      = a_domains;
  m_ref_rat      = a_ref_rat;
  m_finestLevel = a_newFinestLevel;

  m_dx.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    m_dx[lvl] = a_dx[lvl]*RealVect::Unit;
  }
  
  this->setupPVR(a_lmin, m_finestLevel);
  this->setupParticleData(a_lmin, m_finestLevel);
  
  const int oldFinestLevel = m_cacheParticles.size() - 1;

  // If a level was removed, put the particles on the new finest level
  if(m_finestLevel < oldFinestLevel){
    List<P>& cacheFinest = m_cacheParticles[m_finestLevel];
    for (int lvl = m_finestLevel + 1; lvl <= oldFinestLevel; lvl++){
      List<P>& p = m_cacheParticles[lvl];

      cacheFinest.catenate(p);
      p.clear();
    }
  }
  
  for (int lvl = 0; lvl <= a_newFinestLevel; lvl++){
    List<P>& outcast = m_particles[lvl]->outcast();

    // Add particles from cache 
    if(lvl <= oldFinestLevel){
      outcast.catenate(m_cacheParticles[lvl]);
      m_cacheParticles[lvl].clear();
    }
    
    m_particles[lvl]->remapOutcast();
  }
  
  // If grids were removed the fine particle outcast list is definitely not empty. Transfer fine particles
  // to the outcast list below and remap the coarse level again. Some of these particles may be outside the
  // coar level PVR, but that is taken care of later. After this, only the outcast list on lvl = 0 has particles. 
  for (int lvl = a_newFinestLevel; lvl > 0; lvl--){
    List<P>& fineParticles = m_particles[lvl]->outcast();
    List<P>& coarParticles = m_particles[lvl-1]->outcast();

    coarParticles.catenate(fineParticles);
    fineParticles.clear();

    m_particles[lvl-1]->remapOutcast();
  }
  
  // Now do remap so that the PVR is respected.
  this->remap();

  m_cacheParticles.resize(0);
}

template <class P>
void ParticleContainer<P>::discardInvalidParticles(){
  if(m_cell_sorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");  
  // Move particles to the outcast list on the fine level. Then discard that list. 
  for (int lvl = 1; lvl <= m_finestLevel; lvl++){
    collectValidParticles(m_particles[lvl]->outcast(),
			  *m_particles[lvl-1],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  m_ref_rat[lvl-1],
			  false,
			  m_probLo);
    (m_particles[lvl]->outcast()).clear();
  }
}

template <class P>
size_t ParticleContainer<P>::getNumberOfValidParticesLocal() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    n += m_particles[lvl]->numValidLocal();
  }

  return n;
}

template <class P>
size_t ParticleContainer<P>::getNumberOfValidParticesGlobal() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    n += m_particles[lvl]->numValid();
  }

  return n;
}

template <class P>
size_t ParticleContainer<P>::getNumberOfOutcastParticesLocal() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    n += m_particles[lvl]->numOutcastLocal();
  }

  return n;
}

template <class P>
size_t ParticleContainer<P>::getNumberOfOutcastParticesGlobal() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    n += m_particles[lvl]->numOutcast();
  }

  return n;
}

template <class P>
size_t ParticleContainer<P>::getNumberOfMaskParticesLocal() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    n += m_halo_particles[lvl]->numValidLocal();
  }

  return n;
}

template <class P>
size_t ParticleContainer<P>::getNumberOfMaskParticesGlobal() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    n += m_halo_particles[lvl]->numValid();
  }

  return n;
}

template <class P>
Real ParticleContainer<P>::getWeightValidLocal() const {

  Real weight = 0.0;
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit){
      const List<P>& particles = (*m_particles[lvl])[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit){
	weight += lit().mass();
      }
    }
  }
  
  return weight;
}

template <class P>
Real ParticleContainer<P>::getWeightValidGlobal() const{

  Real weight = this->getWeightValidLocal();

#ifdef CH_MPI
  Real tmp = weight;
  int result = MPI_Allreduce(&tmp, &weight, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("particle_Container<P>::getWeightValidGlobal - communication error!");
  }
#endif

  return weight;
}

template <class P>
void ParticleContainer<P>::copyMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool> > > >& a_halo) const{
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    if(lvl < m_finestLevel){
      this->copyMaskParticles(lvl, *a_halo[lvl]);
    }
  }
}

template <class P>
void ParticleContainer<P>::copyNonMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool> > > >& a_halo) const{
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    if(lvl < m_finestLevel){
      this->copyNonMaskParticles(lvl, *a_halo[lvl]);
    }
  }
}

template <class P>
void ParticleContainer<P>::copyMaskParticles(const int a_level, const LevelData<BaseFab<bool> >& a_halo) const {
  
  m_halo_particles[a_level]->clear();

  const RealVect dx = m_dx[a_level];

  // 1. Copy particles from m_particles to m_halo_particles if they lie in the halo mask.
  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit){
    const BaseFab<bool>& mask = a_halo[dit()];
    const List<P>& particles  = (*m_particles[a_level])[dit()].listItems();
    List<P>& halo_particles   = (*m_halo_particles[a_level])[dit()].listItems();
    
    for (ListIterator<P> lit(particles); lit.ok(); ++lit){
      const RealVect x = lit().position();
      const IntVect iv(D_DECL((int)floor((x[0] - m_probLo[0])/dx[0]),
			      (int)floor((x[1] - m_probLo[1])/dx[1]),
			      (int)floor((x[2] - m_probLo[2])/dx[2])));
	
      if(mask(iv)) halo_particles.add(lit());
    }
  }
}

template <class P>
void ParticleContainer<P>::copyNonMaskParticles(const int a_level, const LevelData<BaseFab<bool> >& a_halo) const {
  
  m_non_halo_particles[a_level]->clear();

  const RealVect dx = m_dx[a_level];

  // 1. Copy particles from m_particles to m_non_halo_particles if they don't lie in the halo mask
  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit){
    const BaseFab<bool>& mask   = a_halo[dit()];
    const List<P>& particles    = (*m_particles[a_level])[dit()].listItems();
    List<P>& non_halo_particles = (*m_non_halo_particles[a_level])[dit()].listItems();

    for (ListIterator<P> lit(particles); lit.ok(); ++lit){
      const RealVect x = lit().position();
      const IntVect iv(D_DECL((int)floor((x[0] - m_probLo[0])/dx[0]),
			      (int)floor((x[1] - m_probLo[1])/dx[1]),
			      (int)floor((x[2] - m_probLo[2])/dx[2])));

      if(!mask(iv)) non_halo_particles.add(lit());
    }
  }
}

template <class P>
void ParticleContainer<P>::clearParticles() {
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    m_particles[lvl]->clear();
  }
}

template <class P>
void ParticleContainer<P>::clearMaskParticles() const {
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    m_halo_particles[lvl]->clear();
  }
}

template <class P>
void ParticleContainer<P>::clearNonMaskParticles() const {
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    m_non_halo_particles[lvl]->clear();
  }
}

#include <CD_NamespaceFooter.H>

#endif
