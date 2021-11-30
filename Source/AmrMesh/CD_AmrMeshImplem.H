/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_AmrMeshImplem.H
  @brief  Implementation of CD_AmrMesh.H
  @author Robert Marskar
*/

#ifndef CD_AmrMeshImplem_H
#define CD_AmrMeshImplem_H

// Our includes
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>
  
template<typename T>
void AmrMesh::deallocate(Vector<T*>& a_data) const {
  CH_TIME("AmrMesh::deallocate(Vector<T*>)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::deallocate(Vector<T*>)" << endl;
  }
  
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    delete a_data[lvl];
  }
}

template<typename T>
void AmrMesh::deallocate(Vector<RefCountedPtr<T> >& a_data) const {
  CH_TIME("AmrMesh::deallocate(Vector<RefCountedPtr<T>)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::deallocate(Vector<RefCountedPtr<T>)" << endl;
  }
  
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    //    delete &a_data[lvl];
    a_data[lvl] = RefCountedPtr<T>(0);
#if 0
    if(!a_data[lvl].isNull()){
      delete &(*a_data[lvl]);
      a_data[lvl] = RefCountedPtr<T> (NULL);
    }

#endif
  }
}

template<typename T>
void AmrMesh::deallocate(EBAMRData<T>& a_data) const {
  CH_TIME("AmrMesh::deallocate(EBAMRData<T>)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::deallocate(EBAMRData<T>)" << endl;
  }
  
  return this->deallocate(a_data.getData());
}

template<typename T>
void AmrMesh::alias(Vector<T*>& a_alias, const Vector<RefCountedPtr<T> >& a_data) const {
  CH_TIME("AmrMesh::alias(Vector<T*>, Vector<RefCountedPtr<T>)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::alias(Vector<T*>, Vector<RefCountedPtr<T>)" << endl;
  }

  a_alias.resize(a_data.size());
  
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    a_alias[lvl] = &(*a_data[lvl]);
  }
}

template<typename T, typename S>
void AmrMesh::alias(Vector<T*>& a_alias, const EBAMRData<S>& a_data) const {
  CH_TIME("AmrMesh::alias(Vector<T*>, EBAMRData<S>)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::alias(Vector<T*>, EBAMRData<S>" << endl;
  }
  
  return this->alias(a_alias, a_data.getData());
}

template <typename T>
void AmrMesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >& a_particles, const std::string a_realm) const {
  CH_TIME("AmrMesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >, string)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >, string)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >, string) - could not find Realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  if(m_maxBoxSize != m_blockingFactor){
    MayDay::Abort("AmrMesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >, string) - only constant box sizes supported for particle methods");
  }
  
  a_particles.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()  [lvl];
    const ProblemDomain& domain  = m_realms[a_realm]->getDomains()[lvl];
    
    const Real dx = m_dx[lvl];
    
    a_particles[lvl] = RefCountedPtr<ParticleData<T> > (new ParticleData<T>(dbl,
									    domain,
									    m_blockingFactor,
									    dx*RealVect::Unit,
									    m_probLo));
  }
}

template <typename T>
void AmrMesh::allocate(ParticleContainer<T>& a_container,
		       const int              a_buffer,
		       const std::string      a_realm) const {
  CH_TIME("AmrMesh::allocate(ParticleContainer<T>, int, string)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(ParticleContainer<T>, int, string)" << endl;
  }
  
  if(!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::allocate(ParticleContainer<T>, int, string) - could not find Realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  if(m_maxBoxSize != m_blockingFactor){
    MayDay::Abort("AmrMesh::allocate(ParticleContainer<T>, int, string) - only constant box sizes are supported for particle methods");
  }
  
  a_container.define(m_realms[a_realm]->getGrids(),
		     m_realms[a_realm]->getDomains(),
		     m_realms[a_realm]->getDx(),
		     m_realms[a_realm]->getRefinementRatios(),
		     m_probLo,
		     m_blockingFactor,
		     m_finestLevel,
		     a_buffer,
		     a_realm);
}

template <typename T>
void AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >& a_data) const {
  CH_TIME("AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >)" << endl;
  }

  a_data.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    a_data[lvl] = RefCountedPtr<T> (new T());
  }
}

template <typename T>
void AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >& a_data, const int a_finestLevel) const {
  CH_TIME("AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >, int)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >, int)" << endl;
  }

  a_data.resize(1 + a_finestLevel);
  
  for (int lvl = 0; lvl <= a_finestLevel; lvl++){
    a_data[lvl] = RefCountedPtr<T> (new T());
  }
}

template <typename T>
void AmrMesh::allocatePointer(EBAMRData<T>& a_data) const {
  CH_TIME("AmrMesh::allocatePointer(EBAMRData<T>)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocatePointer(EBAMRData<T>)" << endl;
  }
  
  return this->allocatePointer(a_data.getData());
}

template <typename T>
void AmrMesh::allocatePointer(EBAMRData<T>& a_data, const int a_finestLevel) const {
  CH_TIME("AmrMesh::allocatePointer(EBAMRData<T>, int)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocatePointer(EBAMRData<T>, int)" << endl;
  }
  
  return this->allocatePointer(a_data.getData(), a_finestLevel);
}

template <class P, const Real&(P::*particleScalarField)() const>
void AmrMesh::depositParticles(EBAMRCellData&              a_meshData,
			       const std::string&          a_realm,
			       const phase::which_phase&   a_phase,	       
			       const ParticleContainer<P>& a_particles,
			       const DepositionType        a_depositionType,
			       const CoarseFineDeposition  a_coarseFineDeposition,
			       const bool                  a_forceIrregNGP){
  CH_TIME("AmrMesh::depositParticles");

  EBAMRParticleMesh& particleMesh = m_realms[a_realm]->getParticleMesh(a_phase);

  particleMesh.deposit<P, particleScalarField>(a_meshData, a_particles, a_depositionType, a_coarseFineDeposition, a_forceIrregNGP);
}

template <class P, Real(P::*particleScalarField)() const>
void AmrMesh::depositParticles(EBAMRCellData&              a_meshData,
			       const std::string&          a_realm,
			       const phase::which_phase&   a_phase,	       
			       const ParticleContainer<P>& a_particles,
			       const DepositionType        a_depositionType,
			       const CoarseFineDeposition  a_coarseFineDeposition,
			       const bool                  a_forceIrregNGP){
  CH_TIME("AmrMesh::depositParticles");

  EBAMRParticleMesh& particleMesh = m_realms[a_realm]->getParticleMesh(a_phase);

  particleMesh.deposit<P, particleScalarField>(a_meshData, a_particles, a_depositionType, a_coarseFineDeposition, a_forceIrregNGP);
}

template <class P, const RealVect&(P::*particleVectorField)() const>
void AmrMesh::depositParticles(EBAMRCellData&              a_meshData,
			       const std::string&          a_realm,
			       const phase::which_phase&   a_phase,	       
			       const ParticleContainer<P>& a_particles,
			       const DepositionType        a_depositionType,
			       const CoarseFineDeposition  a_coarseFineDeposition,
			       const bool                  a_forceIrregNGP){
  CH_TIME("AmrMesh::depositParticles");

  EBAMRParticleMesh& particleMesh = m_realms[a_realm]->getParticleMesh(a_phase);

  particleMesh.deposit<P, particleVectorField>(a_meshData, a_particles, a_depositionType, a_coarseFineDeposition, a_forceIrregNGP);
}

template <class P, Real&(P::*particleScalarField)()>
void AmrMesh::interpolateParticles(ParticleContainer<P>&      a_particles,
				   const std::string&         a_realm,
				   const phase::which_phase&  a_phase,	       			    				   
				   const EBAMRCellData&       a_meshScalarField,
				   const DepositionType       a_interpType,
				   const bool                 a_forceIrregNGP) const{
  CH_TIME("AmrMesh::interpolateParticles(scalar)");

  EBAMRParticleMesh& particleMesh = m_realms[a_realm]->getParticleMesh(a_phase);

  particleMesh.interpolate<P, particleScalarField>(a_particles, a_meshScalarField, a_interpType, a_forceIrregNGP);
}

template <class P, RealVect&(P::*particleVectorField)()>
void AmrMesh::interpolateParticles(ParticleContainer<P>&      a_particles,
				   const std::string&         a_realm,
				   const phase::which_phase&  a_phase,	       			    				   
				   const EBAMRCellData&       a_meshVectorField,
				   const DepositionType       a_interpType,
				   const bool                 a_forceIrregNGP) const{
  CH_TIME("AmrMesh::interpolateParticles(vector)");

  EBAMRParticleMesh& particleMesh = m_realms[a_realm]->getParticleMesh(a_phase);

  particleMesh.interpolate<P, particleVectorField>(a_particles, a_meshVectorField, a_interpType, a_forceIrregNGP);  
}

template <class P>
void AmrMesh::removeCoveredParticlesIF(ParticleContainer<P>&     a_particles,
				       const phase::which_phase& a_phase,
				       const Real                a_tolerance) const {
  CH_TIME("AmrMesh::removeCoveredParticlesIF");

  CH_assert(!a_particles.isCellSorted());

  // Figure out the implicit function
  RefCountedPtr<BaseIF> implicitFunction;

  switch(a_phase){
  case phase::gas:
    implicitFunction = m_baseif.at(phase::gas);
    break;
  case phase::solid:
    implicitFunction = m_baseif.at(phase::solid);    
    break;
  default:
    MayDay::Error("AmrMesh::removeCoveredParticlesIF - logic bust");
  }

  // Get the realm where the particles live. 
  const std::string whichRealm = a_particles.getRealm();

  // Go through all particles and remove them if they are less than dx*a_tolerance away from the EB. 
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl   = this->getGrids(whichRealm)[lvl];
    const Real               dx    = this->getDx()[lvl];
    const Real               tol   = a_tolerance*dx;

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      List<P>& particles = a_particles[lvl][dit()].listItems();

      // Check if particles are outside the implicit function. 
      for (ListIterator<P> lit(particles); lit.ok(); ++lit){
	const RealVect& pos = lit().position();
	
	const Real f = implicitFunction->value(pos);
	
	if(f > tol) {
	  particles.remove(lit);
	}
      }
    }
  }
}

template <class P>
void AmrMesh::removeCoveredParticlesDiscrete(ParticleContainer<P>&     a_particles,
					     const phase::which_phase& a_phase,
					     const Real                a_tolerance) const {
  CH_TIME("AmrMesh::removeCoveredParticlesIF");

  CH_assert(!a_particles.isCellSorted());

  // Get the realm where the particles live. 
  const std::string whichRealm = a_particles.getRealm();  

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl   = this->getGrids(whichRealm)[lvl];
    const EBISLayout&        ebisl = this->getEBISLayout(whichRealm, a_phase)[lvl];
    const Real               dx    = this->getDx()[lvl];
    const Real               tol   = a_tolerance*dx;

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const EBISBox& ebisBox = ebisl[dit()];
      const Box      region  = ebisBox.getRegion();

      const bool isRegular   = ebisBox.isAllRegular();
      const bool isCovered   = ebisBox.isAllCovered();
      const bool isIrregular = !isRegular && !isCovered;

      List<P>& particles = a_particles[lvl][dit()].listItems();

      if(isCovered){
	particles.clear();
      }
      else if(isIrregular){
	for (ListIterator<P> lit(particles); lit.ok(); ++lit){
	  const P&        particle = lit();
	  const RealVect& pos      = particle.position();

	  // Get the cell where the particle lives. 
	  const RealVect rv = (pos - m_probLo)/dx;
	  const IntVect  iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

	  CH_assert(region.contains(iv));

	  if(ebisBox.isCovered(iv)){
	    particles.remove(lit());
	  }
	  else if(ebisBox.isIrregular(iv)){
	    // There could potentially be multiple degrees of freedom -- we check if the particle is inside at least one of the VoFs in the cell. 

	    bool insideAtLeastOneVoF = false;

	    // For each VoF, compute the projection of the particle position onth the EB face. Beacuse the EB normal points outwards, we have a positive
	    // value if the particle lies in the valid region. 
	    const std::vector<VolIndex> vofs = ebisBox.getVoFs(iv).stdVector();
	    for (const auto& vof : vofs){
	      const RealVect ebCentroid     = ebisBox.bndryCentroid(vof);
	      const RealVect ebNormal       = ebisBox.normal       (vof);
	      const Real     faceProjection = ebNormal.dotProduct(pos - ebCentroid);

	      if(faceProjection >= tol){
		insideAtLeastOneVoF = true;
		break;
	      }
	    }

	    if(!insideAtLeastOneVoF){
	      particles.remove(lit());
	    }
	  }
	}
      }
    }
  }
}

template <class P>
void AmrMesh::removeCoveredParticlesVoxels(ParticleContainer<P>&     a_particles,
					   const phase::which_phase& a_phase) const {
  CH_TIME("AmrMesh::removeCoveredParticlesVoxels");

  CH_assert(!a_particles.isCellSorted());

  // Get the realm where the particles live. 
  const std::string whichRealm = a_particles.getRealm();  

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl   = this->getGrids(whichRealm)[lvl];
    const EBISLayout&        ebisl = this->getEBISLayout(whichRealm, a_phase)[lvl];
    const Real               dx    = this->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const EBISBox& ebisBox = ebisl[dit()];
      const Box      region  = ebisBox.getRegion();

      const bool isRegular   = ebisBox.isAllRegular();
      const bool isCovered   = ebisBox.isAllCovered();
      const bool isIrregular = !isRegular && !isCovered;

      List<P>& particles = a_particles[lvl][dit()].listItems();

      if(isCovered){
	particles.clear();
      }
      else if(isIrregular){
	for (ListIterator<P> lit(particles); lit.ok(); ++lit){
	  const P&        particle = lit();
	  const RealVect& pos      = particle.position();

	  // Get the cell where the particle lives. 
	  const RealVect rv = (pos - m_probLo)/dx;
	  const IntVect  iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

	  CH_assert(region.contains(iv));

	  if(ebisBox.isCovered(iv)){
	    particles.remove(lit());
	  }
	}
      }
    }
  }  
}

#include <CD_NamespaceFooter.H>

#endif
