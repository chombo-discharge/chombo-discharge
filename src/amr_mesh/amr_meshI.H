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

template<typename T > void amr_mesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >& a_particles){
  CH_TIME("amr_mesh::allocate(AMR ParticleData)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(AMR ParticleData)" << endl;
  }

  if(m_max_box_size != m_blocking_factor){
    MayDay::Abort("amr_mesh::allocate(photons) - only constant box sizes supported for particle methods");
  }
  
  a_particles.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_particles[lvl] = RefCountedPtr<ParticleData<T> > (new ParticleData<T>(m_grids[lvl],
									    m_domains[lvl],
									    m_blocking_factor,
									    m_dx[lvl]*RealVect::Unit,
									    m_prob_lo));
  }
}


template <typename T>
void amr_mesh::allocate(particle_container<T>& a_container, const int a_pvr_buffer){
  if(m_max_box_size != m_blocking_factor){
    MayDay::Abort("amr_mesh::allocate(particle_container<T>) - only constant box sizes are supported for particle methods");
  }
  a_container.define(m_grids,
		     m_domains,
		     m_dx,
		     m_ref_ratios,
		     m_prob_lo,
		     m_blocking_factor,
		     m_finest_level,
		     a_pvr_buffer);
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
