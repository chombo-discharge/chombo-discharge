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
