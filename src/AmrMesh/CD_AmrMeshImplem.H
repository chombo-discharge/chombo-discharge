/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   AmrMeshImplem.H
  @brief  Implementation of AmrMesh.H
  @author Robert Marskar
*/

#ifndef CD_AmrMeshImplem_H
#define CD_AmrMeshImplem_H

// Our includes
#include <CD_NamespaceHeader.H>
  
template<typename T>
void AmrMesh::deallocate(Vector<T*>& a_data){
  CH_TIME("AmrMesh::deallocate");
  if(m_verbosity > 5){
    pout() << "AmrMesh::deallocate" << endl;
  }
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    delete a_data[lvl];
  }
}

template<typename T>
void AmrMesh::deallocate(Vector<RefCountedPtr<T> >& a_data){
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
void AmrMesh::deallocate(EBAMRData<T>& a_data){
  return this->deallocate(a_data.get_data());
}

template<typename T>
void AmrMesh::alias(Vector<T*>& a_alias, const Vector<RefCountedPtr<T> >& a_data){
  CH_TIME("AmrMesh::alias");
  if(m_verbosity > 5){
    pout() << "AmrMesh::alias" << endl;
  }

  a_alias.resize(a_data.size());
  
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    a_alias[lvl] = &(*a_data[lvl]);
  }
}

template<typename T, typename S>
void AmrMesh::alias(Vector<T*>& a_alias, const EBAMRData<S>& a_data){
  return this->alias(a_alias, a_data.get_data());
}

template <typename T>
void AmrMesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >& a_particles, const std::string a_Realm){
  CH_TIME("AmrMesh::allocate(particle_data, Realm)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(particle_data, Realm)" << endl;
  }

  if(!this->queryRealm(a_Realm)) {
    std::string str = "AmrMesh::allocate(particle_data, Realm) - could not find Realm '" + a_Realm + "'";
    MayDay::Abort(str.c_str());
  }

  if(m_maxBoxSize != m_blockingFactor){
    MayDay::Abort("AmrMesh::allocate(particle_data, Realm) - only constant box sizes supported for particle methods");
  }
  
  a_particles.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_Realms[a_Realm]->getGrids()[lvl];
    const ProblemDomain& domain  = m_Realms[a_Realm]->getDomains()[lvl];
    const Real dx = m_dx[lvl];
    a_particles[lvl] = RefCountedPtr<ParticleData<T> > (new ParticleData<T>(dbl,
									    domain,
									    m_blockingFactor,
									    dx*RealVect::Unit,
									    m_probLo));
  }
}

template <typename T>
void AmrMesh::allocate(particle_container<T>& a_container,
		       const int              a_buffer,
		       const std::string      a_Realm){
  CH_TIME("AmrMesh::allocate(container, buffer, Realm)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(container, buffer, Realm)" << endl;
  }
  
  if(!this->queryRealm(a_Realm)) {
    std::string str = "AmrMesh::allocate(container, buffer, Realm) - could not find Realm '" + a_Realm + "'";
    MayDay::Abort(str.c_str());
  }

  if(m_maxBoxSize != m_blockingFactor){
    MayDay::Abort("AmrMesh::allocate(container, buffer, Realm) - only constant box sizes are supported for particle methods");
  }
  
  a_container.define(m_Realms[a_Realm]->getGrids(),
		     m_Realms[a_Realm]->getDomains(),
		     m_Realms[a_Realm]->getDx(),
		     m_Realms[a_Realm]->getRefinementRatios(),
		     m_probLo,
		     m_blockingFactor,
		     m_finestLevel,
		     a_buffer,
		     a_Realm);
}

template <typename T>
void AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >& a_data){
  CH_TIME("AmrMesh::allocatePointer");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocatePointer" << endl;
  }

  a_data.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    a_data[lvl] = RefCountedPtr<T> (new T());
  }
}

template <typename T>
void AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >& a_data, const int a_finestLevel){
  CH_TIME("AmrMesh::allocatePointer");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocatePointer" << endl;
  }

  a_data.resize(1 + a_finestLevel);
  for (int lvl = 0; lvl <= a_finestLevel; lvl++){
    a_data[lvl] = RefCountedPtr<T> (new T());
  }
}

template <typename T>
void AmrMesh::allocatePointer(EBAMRData<T>& a_data){
  return this->allocatePointer(a_data.get_data());
}

template <typename T>
void AmrMesh::allocatePointer(EBAMRData<T>& a_data, const int a_finestLevel){
  return this->allocatePointer(a_data.get_data(), a_finestLevel);
}

#include <CD_NamespaceFooter.H>

#endif
