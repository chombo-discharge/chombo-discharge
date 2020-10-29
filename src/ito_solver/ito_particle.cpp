/*!
  @file   ito_particle.cpp
  @brief  Implementation of ito_particle.H
  @author Robert Marskar
  @date   April 2020
*/

#include "ito_particle.H"

int ito_particle::s_num_runtime_scalars = 0;
int ito_particle::s_num_runtime_vectors = 0;

void ito_particle::set_num_runtime_scalars(const int a_num){
  s_num_runtime_scalars = a_num;
}

void ito_particle::set_num_runtime_vectors(const int a_num){
  s_num_runtime_vectors = a_num;
}

ito_particle::ito_particle() : BinItem(){

  allocateRuntimeBuffers();
}

ito_particle::ito_particle(const ito_particle& a_other) {
  m_position    = a_other.m_position;
  m_mass        = a_other.m_mass;
  m_velocity    = a_other.m_velocity;
  m_diffusion   = a_other.m_diffusion;
  m_oldPosition = a_other.m_oldPosition;
  m_mobility    = a_other.m_mobility;
  m_energy      = a_other.m_energy;

  allocateRuntimeBuffers();
}

ito_particle::ito_particle(const Real      a_mass,
			   const RealVect& a_position,
			   const RealVect& a_velocity,
			   const Real      a_diffusion,
			   const Real      a_mobility,
			   const Real      a_energy) : BinItem(a_position){
  m_mass        = a_mass;
  m_velocity    = a_velocity;
  m_diffusion   = a_diffusion;
  m_oldPosition = a_position;
  m_mobility    = a_mobility;
  m_energy      = a_energy;

  allocateRuntimeBuffers();
}

ito_particle::~ito_particle(){

  if(m_runtimeScalars != nullptr){
    delete [] m_runtimeScalars;
    m_runtimeScalars = nullptr;
  }

  if(m_runtimeVectors != nullptr){
    delete [] m_runtimeVectors;
    m_runtimeVectors = nullptr;
  }
}

void ito_particle::allocateRuntimeBuffers(){
  m_runtimeScalars = nullptr;
  m_runtimeVectors = nullptr;
  
  if(s_num_runtime_scalars >= 0){
    m_runtimeScalars = new Real [s_num_runtime_scalars];
  }

  if(s_num_runtime_vectors >= 0){
    m_runtimeVectors = new RealVect [s_num_runtime_vectors];
  }
}

void ito_particle::define(const Real      a_mass,
			  const RealVect& a_position,
			  const RealVect& a_velocity,
			  const Real      a_diffusion,
			  const Real      a_mobility,
			  const Real      a_energy){
  setMass(a_mass);
  setPosition(a_position);
  setVelocity(a_velocity);
  setDiffusion(a_diffusion);
  setMobility(a_mobility);
  setEnergy(a_energy);
}

void ito_particle::setOldPosition(const RealVect a_oldPosition){
  m_oldPosition = a_oldPosition;
}

RealVect& ito_particle::oldPosition(){
  return m_oldPosition;
}

const RealVect& ito_particle::oldPosition() const{
  return m_oldPosition;
}

void ito_particle::setMass(const Real a_mass){
  m_mass = a_mass;
}

Real& ito_particle::mass(){
  return m_mass;
}

const Real& ito_particle::mass() const{
  return m_mass;
}

void ito_particle::setDiffusion(const Real a_diffusion){
  m_diffusion = a_diffusion;
}

Real& ito_particle::diffusion(){
  return m_diffusion;
}

const Real& ito_particle::diffusion() const{
  return m_diffusion;
}

void ito_particle::setMobility(const Real a_mobility){
  m_mobility = a_mobility;
}


Real& ito_particle::mobility(){
  return m_mobility;
}

const Real& ito_particle::mobility() const{
  return m_mobility;
}

void ito_particle::setEnergy(const Real a_energy){
  m_energy = a_energy;
}

Real& ito_particle::energy(){
  return m_energy;
}

const Real& ito_particle::energy() const{
  return m_energy;
}

void ito_particle::setVelocity(const RealVect& a_velocity){
  m_velocity = a_velocity;
}

void ito_particle::setVelocity(const Real& a_velocity, const int a_dir){
  m_velocity[a_dir] = a_velocity;
}

RealVect& ito_particle::velocity(){
  return m_velocity;
}

const RealVect& ito_particle::velocity() const{
  return m_velocity;
}

Real ito_particle::velocity(const int a_dir) const{
  return m_velocity[a_dir];
}

Real& ito_particle::tmp(){
  return m_tmp;
}

const Real& ito_particle::tmp() const{
  return m_tmp;
}

Real& ito_particle::runtime_scalar(const int a_num){
  return m_runtimeScalars[a_num];
}

const Real& ito_particle::runtime_scalar(const int a_num) const{
  return m_runtimeScalars[a_num];
}

RealVect& ito_particle::runtime_vector(const int a_num){
  return m_runtimeVectors[a_num];
}

const RealVect& ito_particle::runtime_vector(const int a_num) const {
  return m_runtimeVectors[a_num];
}

bool ito_particle::operator==(const ito_particle& a_p) const{
  return ( m_mass      == a_p.m_mass     &&
           m_position  == a_p.m_position &&
           m_velocity  == a_p.m_velocity &&
	   m_diffusion == a_p.m_diffusion);
}

bool ito_particle::operator==(const ito_particle* a_p) const{
  return (*this == *a_p);
}

bool ito_particle::operator!=(const ito_particle& a_p) const{
  return !(*this == a_p);
}

int ito_particle::size() const{
  const int compileSize = BinItem::size() + sizeof(m_mass) + sizeof(m_velocity) + sizeof(m_diffusion) + sizeof(m_oldPosition) + sizeof(m_mobility) + sizeof(m_energy);
  //  const int runtimeSize = sizeof(*m_runtimeScalars) + sizeof(*m_runtimeVectors);
  const int runtimeSize = s_num_runtime_scalars*sizeof(Real) + s_num_runtime_vectors*sizeof(RealVect);

  //  std::cout << s_num_runtime_vectors*sizeof(RealVect) << std::endl;

  return compileSize + runtimeSize;
}

void ito_particle::linearOut(void* buf) const{
  Real* buffer = (Real*)buf;
  D_TERM6( *buffer++ = m_position[0];,
	   *buffer++ = m_position[1];,
	   *buffer++ = m_position[2];,
	   *buffer++ = m_position[3];,
	   *buffer++ = m_position[4];,
	   *buffer++ = m_position[5];);

  D_TERM6( *buffer++ = m_velocity[0];,
	   *buffer++ = m_velocity[1];,
	   *buffer++ = m_velocity[2];,
	   *buffer++ = m_velocity[3];,
	   *buffer++ = m_velocity[4];,
	   *buffer++ = m_velocity[5];);

  D_TERM6( *buffer++ = m_oldPosition[0];,
	   *buffer++ = m_oldPosition[1];,
	   *buffer++ = m_oldPosition[2];,
	   *buffer++ = m_oldPosition[3];,
	   *buffer++ = m_oldPosition[4];,
	   *buffer++ = m_oldPosition[5];);

  *buffer++ = m_diffusion;
  *buffer++ = m_mass;
  *buffer++ = m_mobility;


  // Go through run-time memory
  for (int i = 0; i < s_num_runtime_scalars; i++){
    *buffer++ = m_runtimeScalars[i];
  }

  for (int i = 0; i < s_num_runtime_vectors; i++){
  D_TERM6( *buffer++ = m_runtimeVectors[i][0];,
	   *buffer++ = m_runtimeVectors[i][1];,
	   *buffer++ = m_runtimeVectors[i][2];,
	   *buffer++ = m_runtimeVectors[i][3];,
	   *buffer++ = m_runtimeVectors[i][4];,
	   *buffer++ = m_runtimeVectors[i][5];);
  }

  *buffer = m_energy;
}

void ito_particle::linearIn(void* buf){
  Real* buffer = (Real*)buf;
  D_TERM6( m_position[0] = *buffer++;,
	   m_position[1] = *buffer++;,
	   m_position[2] = *buffer++;,
	   m_position[3] = *buffer++;,
	   m_position[4] = *buffer++;,
	   m_position[5] = *buffer++;);

  D_TERM6( m_velocity[0] = *buffer++;,
	   m_velocity[1] = *buffer++;,
	   m_velocity[2] = *buffer++;,
	   m_velocity[3] = *buffer++;,
	   m_velocity[4] = *buffer++;,
	   m_velocity[5] = *buffer++;);

  D_TERM6( m_oldPosition[0] = *buffer++;,
	   m_oldPosition[1] = *buffer++;,
	   m_oldPosition[2] = *buffer++;,
	   m_oldPosition[3] = *buffer++;,
	   m_oldPosition[4] = *buffer++;,
	   m_oldPosition[5] = *buffer++;);

  m_diffusion = *buffer++;
  m_mass      = *buffer++;
  m_mobility  = *buffer++;

  // Go through run-time memory
  for (int i = 0; i < s_num_runtime_scalars; i++){
    m_runtimeScalars[i] = *buffer++;
  }
  for (int i = 0; i < s_num_runtime_vectors; i++){
    D_TERM6( m_runtimeVectors[i][0] = *buffer++;,
	     m_runtimeVectors[i][1] = *buffer++;,
	     m_runtimeVectors[i][2] = *buffer++;,
	     m_runtimeVectors[i][3] = *buffer++;,
	     m_runtimeVectors[i][4] = *buffer++;,
	     m_runtimeVectors[i][5] = *buffer++;);
  }

  m_energy    = *buffer;
}

std::ostream & operator<<(std::ostream& ostr, const ito_particle& p){
  ostr << " ito_particle : " << std::endl;
  ostr << " mass " << p.mass() << std::endl;
  ostr << " position ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.position(i); }
  ostr << " ) ";
  ostr << std::endl << " velocity ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.velocity(i); }
  ostr << " ) ";
  ostr << std::endl << " diffusion " << p.diffusion() << std::endl;
  ostr << std::endl << " mobility " << p.mobility() << std::endl;
  ostr << std::endl << " energy " << p.energy() << std::endl;
  return ostr;
}
