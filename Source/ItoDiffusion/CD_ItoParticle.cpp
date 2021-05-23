/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoParticle.cpp
  @brief  Implementation of CD_ItoParticle.H
  @author Robert Marskar
*/

// Our includes
#include <CD_ItoParticle.H>
#include <CD_NamespaceHeader.H>
  
int ItoParticle::s_num_runtimeScalars = 0;
int ItoParticle::s_num_runtimeVectors = 0;

void ItoParticle::setNumRunTimeScalars(const int a_num){
  s_num_runtimeScalars = a_num;
}

void ItoParticle::setNumRunTimeVectors(const int a_num){
  s_num_runtimeVectors = a_num;
}

ItoParticle::ItoParticle() : BinItem(){

  allocateRuntimeBuffers();
}

ItoParticle::ItoParticle(const ItoParticle& a_other) {
  m_position    = a_other.m_position;
  m_mass        = a_other.m_mass;
  m_velocity    = a_other.m_velocity;
  m_diffusion   = a_other.m_diffusion;
  m_oldPosition = a_other.m_oldPosition;
  m_mobility    = a_other.m_mobility;
  m_energy      = a_other.m_energy;

  allocateRuntimeBuffers();
}

ItoParticle::ItoParticle(const Real      a_mass,
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

ItoParticle::~ItoParticle(){

  if(m_runtimeScalars != nullptr){
    delete [] m_runtimeScalars;
    m_runtimeScalars = nullptr;
  }

  if(m_runtimeVectors != nullptr){
    delete [] m_runtimeVectors;
    m_runtimeVectors = nullptr;
  }
}

void ItoParticle::allocateRuntimeBuffers(){
  m_runtimeScalars = nullptr;
  m_runtimeVectors = nullptr;
  
  if(s_num_runtimeScalars >= 0){
    m_runtimeScalars = new Real [s_num_runtimeScalars];
  }

  if(s_num_runtimeVectors >= 0){
    m_runtimeVectors = new RealVect [s_num_runtimeVectors];
  }
}

void ItoParticle::define(const Real      a_mass,
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

void ItoParticle::setOldPosition(const RealVect a_oldPosition){
  m_oldPosition = a_oldPosition;
}

RealVect& ItoParticle::oldPosition(){
  return m_oldPosition;
}

const RealVect& ItoParticle::oldPosition() const{
  return m_oldPosition;
}

void ItoParticle::setMass(const Real a_mass){
  m_mass = a_mass;
}

Real& ItoParticle::mass(){
  return m_mass;
}

const Real& ItoParticle::mass() const{
  return m_mass;
}

void ItoParticle::setDiffusion(const Real a_diffusion){
  m_diffusion = a_diffusion;
}

Real& ItoParticle::diffusion(){
  return m_diffusion;
}

const Real& ItoParticle::diffusion() const{
  return m_diffusion;
}

void ItoParticle::setMobility(const Real a_mobility){
  m_mobility = a_mobility;
}


Real& ItoParticle::mobility(){
  return m_mobility;
}

const Real& ItoParticle::mobility() const{
  return m_mobility;
}

void ItoParticle::setEnergy(const Real a_energy){
  m_energy = a_energy;
}

Real& ItoParticle::energy(){
  return m_energy;
}

const Real& ItoParticle::energy() const{
  return m_energy;
}

void ItoParticle::setVelocity(const RealVect& a_velocity){
  m_velocity = a_velocity;
}

void ItoParticle::setVelocity(const Real& a_velocity, const int a_dir){
  m_velocity[a_dir] = a_velocity;
}

RealVect& ItoParticle::velocity(){
  return m_velocity;
}

const RealVect& ItoParticle::velocity() const{
  return m_velocity;
}

Real ItoParticle::velocity(const int a_dir) const{
  return m_velocity[a_dir];
}

Real& ItoParticle::tmp(){
  return m_tmp;
}

const Real& ItoParticle::tmp() const{
  return m_tmp;
}

Real& ItoParticle::runtimeScalar(const int a_num){
  return m_runtimeScalars[a_num];
}

const Real& ItoParticle::runtimeScalar(const int a_num) const{
  return m_runtimeScalars[a_num];
}

RealVect& ItoParticle::runtimeVector(const int a_num){
  return m_runtimeVectors[a_num];
}

const RealVect& ItoParticle::runtimeVector(const int a_num) const {
  return m_runtimeVectors[a_num];
}

bool ItoParticle::operator==(const ItoParticle& a_p) const{
  return ( m_mass      == a_p.m_mass     &&
	   m_position  == a_p.m_position &&
	   m_velocity  == a_p.m_velocity &&
	   m_diffusion == a_p.m_diffusion);
}

bool ItoParticle::operator==(const ItoParticle* a_p) const{
  return (*this == *a_p);
}

bool ItoParticle::operator!=(const ItoParticle& a_p) const{
  return !(*this == a_p);
}

int ItoParticle::size() const{
  const int compileSize = BinItem::size() + sizeof(m_mass) + sizeof(m_velocity) + sizeof(m_diffusion) + sizeof(m_oldPosition) + sizeof(m_mobility) + sizeof(m_energy);
  //  const int runtimeSize = sizeof(*m_runtimeScalars) + sizeof(*m_runtimeVectors);
  const int runtimeSize = s_num_runtimeScalars*sizeof(Real) + s_num_runtimeVectors*sizeof(RealVect);

  //  std::cout << s_num_runtimeVectors*sizeof(RealVect) << std::endl;

  return compileSize + runtimeSize;
}

void ItoParticle::linearOut(void* buf) const{
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
  for (int i = 0; i < s_num_runtimeScalars; i++){
    *buffer++ = m_runtimeScalars[i];
  }

  for (int i = 0; i < s_num_runtimeVectors; i++){
    D_TERM6( *buffer++ = m_runtimeVectors[i][0];,
	     *buffer++ = m_runtimeVectors[i][1];,
	     *buffer++ = m_runtimeVectors[i][2];,
	     *buffer++ = m_runtimeVectors[i][3];,
	     *buffer++ = m_runtimeVectors[i][4];,
	     *buffer++ = m_runtimeVectors[i][5];);
  }

  *buffer = m_energy;
}

void ItoParticle::linearIn(void* buf){
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
  for (int i = 0; i < s_num_runtimeScalars; i++){
    m_runtimeScalars[i] = *buffer++;
  }
  for (int i = 0; i < s_num_runtimeVectors; i++){
    D_TERM6( m_runtimeVectors[i][0] = *buffer++;,
	     m_runtimeVectors[i][1] = *buffer++;,
	     m_runtimeVectors[i][2] = *buffer++;,
	     m_runtimeVectors[i][3] = *buffer++;,
	     m_runtimeVectors[i][4] = *buffer++;,
	     m_runtimeVectors[i][5] = *buffer++;);
  }

  m_energy    = *buffer;
}

std::ostream & operator<<(std::ostream& ostr, const ItoParticle& p){
  ostr << " ItoParticle : " << std::endl;
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

#include <CD_NamespaceFooter.H>
