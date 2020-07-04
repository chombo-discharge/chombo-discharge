/*!
  @file amr_meshI.H
  @brief Implementation of template functions in amr_mesh:H
  @author Robert Marskar
  @date Jan. 2018
*/

template<typename T> void amr_mesh::deallocate(Vector<T*>& a_data){
  CH_TIME("amr_mesh::deallocate");
  if(m_verbosity > 5){
    pout() << "amr_mesh::deallocate" << endl;
  }
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    delete a_data[lvl];
  }
}

template<typename T> void amr_mesh::deallocate(Vector<RefCountedPtr<T> >& a_data){
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

template<typename T> void amr_mesh::alias(Vector<T*>& a_alias, const Vector<RefCountedPtr<T> >& a_data){
  CH_TIME("amr_mesh::alias");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias" << endl;
  }

  a_alias.resize(a_data.size());
  
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    a_alias[lvl] = &(*a_data[lvl]);
  }
}

template<typename T> void amr_mesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >& a_particles){
  CH_TIME("amr_mesh::allocate(AMR ParticleData)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(AMR ParticleData)" << endl;
  }

  if(m_max_box_size != m_blocking_factor){
    MayDay::Abort("amr_mesh::allocate(photons) - only constant box sizes supported for particle methods");
  }
  
  a_particles.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realm->get_grids()[lvl];
    const ProblemDomain& domain  = m_realm->get_domains()[lvl];
    const Real dx = m_realm->get_dx()[lvl];
    a_particles[lvl] = RefCountedPtr<ParticleData<T> > (new ParticleData<T>(dbl,
									    domain,
									    m_blocking_factor,
									    dx*RealVect::Unit,
									    m_prob_lo));
  }
}

template <typename T> void amr_mesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >& a_particles, const std::string a_realm){
  CH_TIME("amr_mesh::allocate(particle_data, realm)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(particle_data, realm)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::allocate(particle_data, realm) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  if(m_max_box_size != m_blocking_factor){
    MayDay::Abort("amr_mesh::allocate(particle_data, realm) - only constant box sizes supported for particle methods");
  }
  
  a_particles.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    const ProblemDomain& domain  = m_realms[a_realm]->get_domains()[lvl];
    const Real dx = m_dx[lvl];
    a_particles[lvl] = RefCountedPtr<ParticleData<T> > (new ParticleData<T>(dbl,
									    domain,
									    m_blocking_factor,
									    dx*RealVect::Unit,
									    m_prob_lo));
  }
}
																

template <typename T> void amr_mesh::allocate(particle_container<T>& a_container, const int a_pvr_buffer){
  CH_TIME("amr_mesh::allocate(container, buffer)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(container, buffer)" << endl;
  }
  
  if(m_max_box_size != m_blocking_factor){
    MayDay::Abort("amr_mesh::allocate(particle_container<T>) - only constant box sizes are supported for particle methods");
  }
  
  a_container.define(m_realm->get_grids(),
		     m_realm->get_domains(),
		     m_realm->get_dx(),
		     m_realm->get_ref_rat(),
		     m_prob_lo,
		     m_blocking_factor,
		     m_finest_level,
		     a_pvr_buffer);
}

template <typename T> void amr_mesh::allocate(particle_container<T>& a_container, const int a_buffer, const std::string a_realm){
  CH_TIME("amr_mesh::allocate(container, buffer, realm)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(container, buffer, realm)" << endl;
  }
  
  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::allocate(container, buffer, realm) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  if(m_max_box_size != m_blocking_factor){
    MayDay::Abort("amr_mesh::allocate(container, buffer, realm) - only constant box sizes are supported for particle methods");
  }
  
  a_container.define(m_realms[a_realm]->get_grids(),
		     m_realms[a_realm]->get_domains(),
		     m_realms[a_realm]->get_dx(),
		     m_realms[a_realm]->get_ref_rat(),
		     m_prob_lo,
		     m_blocking_factor,
		     m_finest_level,
		     a_buffer);
}

template<typename T > void amr_mesh::allocate_ptr(Vector<RefCountedPtr<T> >& a_data){
  CH_TIME("amr_mesh::allocate_ptr");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate_ptr" << endl;
  }

  a_data.resize(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl] = RefCountedPtr<T> (new T());
  }
}
