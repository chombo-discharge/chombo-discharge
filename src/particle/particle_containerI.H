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
					  const RealVect&                  a_origin){

  m_grids   = a_grids;
  m_domains = a_domains;
  m_ref_rat = a_ref_rat;
  m_origin  = a_origin;

  m_dx.resize(a_dx.size());
  for (int lvl = 0; lvl < a_dx.size(); lvl++){
    m_dx[lvl] = a_dx[lvl]*RealVect::Unit;
  }
}

template <class P>
particle_container<P>::~particle_container(){

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
  MayDay::Abort("particle_container::add_particles - not implemented");
}

template <class P>
void particle_container<P>::add_particles_destructive(List<P>& a_particles){
  MayDay::Abort("particle_container::add_particles_destructive - not implemented");
}
