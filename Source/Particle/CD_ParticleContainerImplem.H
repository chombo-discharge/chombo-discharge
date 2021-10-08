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
#include <CH_Timer.H>

// Our includes
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>
  
template <class P>
ParticleContainer<P>::ParticleContainer(){
  m_isDefined    = false;
  m_isCellSorted = false;
}

template <class P>
ParticleContainer<P>::ParticleContainer(const Vector<DisjointBoxLayout>& a_grids,
					const Vector<ProblemDomain>&     a_domains,
					const Vector<Real>&              a_dx,
					const Vector<int>&               a_refRat,
					const RealVect&                  a_probLo,
					const int                        a_blockingFactor,
					const int                        a_finestLevel,
					const int                        a_pvrBuffer,
					const std::string                a_realm){
  CH_TIME("ParticleContainer<P>::ParticleContainer");

  this->define(a_grids, a_domains, a_dx, a_refRat, a_probLo, a_blockingFactor, a_finestLevel, a_pvrBuffer, a_realm);
}

template <class P>
ParticleContainer<P>::~ParticleContainer(){
  CH_TIME("ParticleContainer<P>::~ParticleContainer");
}

template <class P>
void ParticleContainer<P>::define(const Vector<DisjointBoxLayout>& a_grids,
				  const Vector<ProblemDomain>&     a_domains,
				  const Vector<Real>&              a_dx,
				  const Vector<int>&               a_refRat,
				  const RealVect&                  a_probLo,
				  const int                        a_blockingFactor,
				  const int                        a_finestLevel,
				  const int                        a_pvrBuffer,
				  const std::string                a_realm){
  CH_TIME("ParticleContainer<P>::~ParticleContainer");

  CH_assert(a_pvrBuffer >= 0);
  
  m_grids           = a_grids;
  m_domains         = a_domains;
  m_refRat          = a_refRat;
  m_probLo          = a_probLo;
  m_finestLevel     = a_finestLevel;
  m_pvrBuffer       = a_pvrBuffer;
  m_realm           = a_realm;
  m_blockingFactor  = a_blockingFactor;

  m_dx.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    m_dx[lvl] = a_dx[lvl]*RealVect::Unit;
  }

  constexpr int base = 0;

  // Do the define stuff.
  this->setupGrownGrids  (base, m_finestLevel);
  this->setupPVR         (base, m_finestLevel);
  this->setupParticleData(base, m_finestLevel);

  m_isDefined    = true;
  m_isCellSorted = false;
}

template <class P>
void ParticleContainer<P>::setupGrownGrids(const int a_base, const int a_finestLevel){
  CH_TIME("ParticleContainer<P>::setupGrownGrids");

  // TLDR: This routine sets up the buffer grids which can be used when we need to fetch particles that fall off each levels' grid. This is very useful
  //       when we want to fetch particles that lie on the coarse grid of a refinement boundary.
  
  m_grownGrids.resize(1 + a_finestLevel);

  for (int lvl = 0; lvl <= a_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_grids  [lvl];
    const ProblemDomain& domain  = m_domains[lvl];

    // Grow boxes by refinement factor on the finer levels. 
    Vector<Box> boxes = dbl.boxArray();    
    if(lvl > 0){
      for (auto& box : boxes.stdVector()){
	box.grow(m_refRat[lvl-1]);
	box &= domain;
      }
    }

    m_grownGrids[lvl].define(boxes, dbl.procIDs());
  }
}

template <class P>
void ParticleContainer<P>::setupPVR(const int a_base, const int a_finestLevel){
  CH_TIME("ParticleContainer<P>::setupGrownGrids");

  // TLDR: This routine sets up the ParticleValidRegion. The PVR is used when we want to ensure that particles that live on a specific grid
  //       level are a specified distance away from the refinement boundary. 
  m_pvr.resize(1 + a_finestLevel);
  for (int lvl = a_base; lvl <= a_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_grids  [lvl];
    const ProblemDomain& domain  = m_domains[lvl];

    const bool hasCoar = lvl > 0;

    if(!hasCoar){
      m_pvr[lvl] = RefCountedPtr<ParticleValidRegion> (new ParticleValidRegion(dbl, nullptr, 1, 0));
    }
    else{
      const int ref_ratio = m_refRat[lvl-1];
      m_pvr[lvl] = RefCountedPtr<ParticleValidRegion> (new ParticleValidRegion(dbl, m_pvr[lvl-1]->mask(), m_refRat[lvl-1], m_pvrBuffer));
    }
  }
}

template <class P>
void ParticleContainer<P>::setupParticleData(const int a_base, const int a_finestLevel){
  CH_TIME("ParticleContainer<P>::setupParticleData");

  // TLDR: This sets up the most commonly used particle data holders for this particle AMR container. This means that
  //       we allocate:
  //
  //          1. The "normal" particle container data holder m_particles. This is defined on the DisjointBoxLayout
  //             which is the natural place for the particles to live.  
  //
  //          2. A buffer particle data holder which is defined on a grown DisjointBoxLayout. This is useful when particles
  //             on the coarse level need to deposit to the fine level.
  //
  //          3. A data holder for "masked particles", providing an opportunity to copy/transfer some of the particles on a grid
  //             level to a separate data holder if they lie within a "mask". Typically used for extracting coarse-level that live
  //             just outside the fine grid (i.e. on the coarse side of the refinement boundary).
  //
  //          4. A data holder for "non-masked particles". This is the same as for the "masked particles", except that it is used for
  //             the inverse mask. Together, m_maskParticles and m_nonMaskParticles provide an opportunity for segregating particles
  //             that live on cells where the mask is true/false. 
  //
  //          5. A data holder for storing cell-sorted particles. Very useful when particles need to be sorted by cell rather than patch. 

  
  m_particles.          resize(1 + a_finestLevel);
  m_bufferParticles.    resize(1 + a_finestLevel);  
  m_maskParticles.      resize(1 + a_finestLevel);
  m_nonMaskParticles.   resize(1 + a_finestLevel);

  m_cellSortedParticles.resize(1 + a_finestLevel);
  
  for (int lvl = a_base; lvl <= a_finestLevel; lvl++){
    m_particles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
									    m_domains[lvl],
									    m_blockingFactor,
									    m_dx[lvl],
									    m_probLo));

    m_bufferParticles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grownGrids[lvl],
										  m_domains[lvl],
										  m_blockingFactor,
										  m_dx[lvl],
										  m_probLo));    

    m_maskParticles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
										m_domains[lvl],
										m_blockingFactor,
										m_dx[lvl],
										m_probLo));

    m_nonMaskParticles[lvl] = RefCountedPtr<ParticleData<P> > (new ParticleData<P>(m_grids[lvl],
										   m_domains[lvl],
										   m_blockingFactor,
										   m_dx[lvl],
										   m_probLo));


    m_cellSortedParticles[lvl] = RefCountedPtr<LayoutData<BinFab<P> > > (new LayoutData<BinFab<P> >(m_grids[lvl]));
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
  CH_assert(m_isDefined);
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::getParticles - particles are sorted by cell!");
  return m_particles;
}

template <class P>
const AMRParticles<P>& ParticleContainer<P>::getParticles() const {
  CH_assert(m_isDefined);
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::getParticles - particles are sorted by cell!");
  return m_particles;
}

template <class P>
AMRParticles<P>& ParticleContainer<P>::getBufferParticles(){
  CH_assert(m_isDefined);

  return m_particles;
}

template <class P>
const AMRParticles<P>& ParticleContainer<P>::getBufferParticles() const {
  CH_assert(m_isDefined);

  return m_particles;
}

template <class P>
AMRParticles<P>& ParticleContainer<P>::getMaskParticles(){
  CH_assert(m_isDefined);
  return m_maskParticles;
}

template <class P>
const AMRParticles<P>& ParticleContainer<P>::getMaskParticles() const {
  CH_assert(m_isDefined);
  return m_maskParticles;
}

template <class P>
AMRParticles<P>& ParticleContainer<P>::getNonMaskParticles(){
  CH_assert(m_isDefined);
  return m_nonMaskParticles;
}

template <class P>
const AMRParticles<P>& ParticleContainer<P>::getNonMaskParticles() const {
  CH_assert(m_isDefined);
  return m_nonMaskParticles;
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
  CH_assert(m_isDefined);
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::operator[](a_lvl) - particles are sorted by cell!");
  return *m_particles[a_lvl];
}

template <class P>
const ParticleData<P>& ParticleContainer<P>::operator[](const int a_level) const {
  CH_assert(m_isDefined);
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::operator[](a_lvl) - particles are sorted by cell!");
  return *m_particles[a_level];
}

template <class P>
AMRCellParticles<P>& ParticleContainer<P>::getCellParticles(){
  if(!m_isCellSorted) MayDay::Abort("ParticleContainer<P>::getCellParticles()- particles are not sorted by cell!");
  return m_cellSortedParticles;
}

template <class P>
const AMRCellParticles<P>& ParticleContainer<P>::getCellParticles() const{
  if(!m_isCellSorted) MayDay::Abort("ParticleContainer<P>::getCellParticles()- particles are not sorted by cell!");
  return m_cellSortedParticles;
}
  
template <class P>
LayoutData<BinFab<P> >& ParticleContainer<P>::getCellParticles(const int a_level){
  if(!m_isCellSorted) MayDay::Abort("ParticleContainer<P>::getCellParticles(level)- particles are not sorted by cell!");
  return *m_cellSortedParticles[a_level];
}

template <class P>
const LayoutData<BinFab<P> >& ParticleContainer<P>::getCellParticles(const int a_level) const{
  if(!m_isCellSorted) MayDay::Abort("ParticleContainer<P>::getCellParticles(level)- particles are not sorted by cell!");
  return *m_cellSortedParticles[a_level];
}

template <class P>
void ParticleContainer<P>::getCellParticles(BinFab<P>& cellParticles, const int a_lvl, const DataIndex a_dit) const {
  CH_assert(m_isDefined);
  if(!m_isCellSorted) MayDay::Abort("ParticleContainer<P>::getCellParticles - particles are not sorted by cell!");
  
  cellParticles.define(m_grids[a_lvl].get(a_dit),  m_dx[a_lvl], m_probLo);
  cellParticles.addItems((*m_particles[a_lvl])[a_dit].listItems());
}

template <class P>
void ParticleContainer<P>::getCellParticlesDestructive(BinFab<P>& cellParticles, const int a_lvl, const DataIndex a_dit){
  CH_assert(m_isDefined);
  if(!m_isCellSorted) MayDay::Abort("ParticleContainer<P>::getCellParticlesDestructive - particles are not sorted by cell!");
  
  cellParticles.define(m_grids[a_lvl].get(a_dit),  m_dx[a_lvl], m_probLo);
  cellParticles.addItemsDestructive((*m_particles[a_lvl])[a_dit].listItems());
}

template <class P>
BinFab<P>& ParticleContainer<P>::getCellParticles(const int a_level, const DataIndex a_dit){
  return (*m_cellSortedParticles[a_level])[a_dit];
}

template <class P>
const BinFab<P>& ParticleContainer<P>::getCellParticles(const int a_level, const DataIndex a_dit) const{
  return (*m_cellSortedParticles[a_level])[a_dit];
}

template <class P>
void ParticleContainer<P>::sortParticlesByCell(){
  if(!m_isCellSorted){

    for (int lvl = 0; lvl <= m_finestLevel; lvl++){
      const DisjointBoxLayout& dbl = m_grids[lvl];
      
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	BinFab<P>& cellParticles = (*m_cellSortedParticles[lvl])[dit()];
	cellParticles.define(dbl.get(dit()), m_dx[lvl], m_probLo);

	cellParticles.addItemsDestructive((*m_particles[lvl])[dit()].listItems());
      }
    }
    
    m_isCellSorted = true;
  }
}

template <class P>
void ParticleContainer<P>::sortParticlesByPatch(){
  if(m_isCellSorted){

    const int comp = 0;
    
    for (int lvl = 0; lvl <= m_finestLevel; lvl++){
      const DisjointBoxLayout& dbl = m_grids[lvl];
      
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	ListBox<P>& patchParticles = (*m_particles[lvl])[dit()];
	BinFab<P>& cellParticles   = (*m_cellSortedParticles[lvl])[dit()];
	
	const Box box = dbl.get(dit());
	for (BoxIterator bit(box); bit.ok(); ++bit){
	  const IntVect iv = bit();
	  
	  patchParticles.addItemsDestructive(cellParticles(iv, comp));

	}

	cellParticles.clear();
      }
    }
    
    m_isCellSorted = false;
  }
}

template <class P>
void ParticleContainer<P>::addParticles(const List<P>& a_particles){
  CH_assert(m_isDefined);
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(List<P>) - particles are sorted by cell!");
  
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
			  m_refRat[lvl-1],
			  false,
			  m_probLo);
    m_particles[lvl]->remapOutcast();
  }
}

template <class P>
void ParticleContainer<P>::addParticlesDestructive(List<P>& a_particles){
  CH_assert(m_isDefined);
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(List<P>) - particles are sorted by cell!");

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
			  m_refRat[lvl-1],
			  false,
			  m_probLo);
    m_particles[lvl]->remapOutcast();
  }
}

template <class P>
void ParticleContainer<P>::addParticles(const BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit){
  CH_assert(m_isDefined);
  CH_assert(m_grids[a_lvl].get(a_dit) == a_particles.getRegion());
  if(!m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(BinFab<P>) - particles are not sorted by cell!");

  const int comp = 0;
  
  BinFab<P>& boxParticles = (*m_cellSortedParticles[a_lvl])[a_dit];

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
  CH_assert(m_isDefined);
  CH_assert(m_grids[a_lvl].get(a_dit) == a_particles.getRegion());
  if(!m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(BinFab<P>) - particles are not sorted by cell!");

  const int comp = 0;

  BinFab<P>& boxParticles = (*m_cellSortedParticles[a_lvl])[a_dit];

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
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");

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
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");

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
  CH_assert(m_isDefined);
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  
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
  CH_assert(m_isDefined);
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const bool hasCoar = lvl > 0;

    // Rebin this level
    m_particles[lvl]->gatherOutcast();
    m_particles[lvl]->remapOutcast();

    // Collect coarser level particles onto this levels outcast list if they fit in this levels PVR
    if(hasCoar){
      collectValidParticles(m_particles[lvl]->outcast(),
			    *m_particles[lvl-1],
			    m_pvr[lvl]->mask(),
			    m_dx[lvl],
			    m_refRat[lvl],
			    false, 
			    m_probLo);

      // Remap new outcasts
      m_particles[lvl]->remapOutcast();
    }


    // Transfer finer-level particles
    if(hasCoar){
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
  CH_assert(m_isDefined);
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  
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
			  m_refRat[lvl-1],
			  false, 
			  m_probLo);
    m_particles[lvl]->remapOutcast();
  }
}

template <class P>
void ParticleContainer<P>::levelRemap(){
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    this->levelRemap(lvl);
  }
}

template <class P>
void ParticleContainer<P>::levelRemap(const int a_lvl){
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  m_particles[a_lvl]->gatherOutcast();
  m_particles[a_lvl]->remapOutcast();
}

template <class P>
void ParticleContainer<P>::preRegrid(const int a_lmin){
  CH_assert(m_isDefined);
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");

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
				  const Vector<int>&               a_refRat,
				  const int                        a_lmin,
				  const int                        a_newFinestLevel){
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::regrid(...) - particles are sorted by cell!");

  // Update this stuff
  m_grids        = a_grids;
  m_domains      = a_domains;
  m_refRat      = a_refRat;
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
  if(m_isCellSorted) MayDay::Abort("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");  
  // Move particles to the outcast list on the fine level. Then discard that list. 
  for (int lvl = 1; lvl <= m_finestLevel; lvl++){
    collectValidParticles(m_particles[lvl]->outcast(),
			  *m_particles[lvl-1],
			  m_pvr[lvl]->mask(),
			  m_dx[lvl],
			  m_refRat[lvl-1],
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
    n += m_maskParticles[lvl]->numValidLocal();
  }

  return n;
}

template <class P>
size_t ParticleContainer<P>::getNumberOfMaskParticesGlobal() const {
  size_t n = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    n += m_maskParticles[lvl]->numValid();
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
  
  m_maskParticles[a_level]->clear();

  const RealVect dx = m_dx[a_level];

  // 1. Copy particles from m_particles to m_maskParticles if they lie in the halo mask.
  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit){
    const BaseFab<bool>& mask = a_halo[dit()];
    const List<P>& particles  = (*m_particles[a_level])[dit()].listItems();
    List<P>& halo_particles   = (*m_maskParticles[a_level])[dit()].listItems();
    
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
  
  m_nonMaskParticles[a_level]->clear();

  const RealVect dx = m_dx[a_level];

  // 1. Copy particles from m_particles to m_nonMaskParticles if they don't lie in the halo mask
  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit){
    const BaseFab<bool>& mask   = a_halo[dit()];
    const List<P>& particles    = (*m_particles[a_level])[dit()].listItems();
    List<P>& non_halo_particles = (*m_nonMaskParticles[a_level])[dit()].listItems();

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
    m_maskParticles[lvl]->clear();
  }
}

template <class P>
void ParticleContainer<P>::clearNonMaskParticles() const {
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    m_nonMaskParticles[lvl]->clear();
  }
}

#include <CD_NamespaceFooter.H>

#endif
