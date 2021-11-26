/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoParticle.cpp
  @brief  Implementation of CD_ItoParticle.H
  @author Robert Marskar
*/

#ifndef CD_ItoParticleImplem_H
#define CD_ItoParticleImplem_H

// Our includes
#include <CD_ItoParticle.H>
#include <CD_NamespaceHeader.H>

void ItoParticle::define(const Real      a_mass,
			 const RealVect& a_position,
			 const RealVect& a_velocity,
			 const Real      a_diffusion,
			 const Real      a_mobility,
			 const Real      a_averageEnergy){
  this->setMass     (a_mass         );
  this->setPosition (a_position     );
  this->setVelocity (a_velocity     );
  this->setDiffusion(a_diffusion    );
  this->setMobility (a_mobility     );
  this->setEnergy   (a_averageEnergy);
}

Real& ItoParticle::mass(){
  return m_mass;
}

const Real& ItoParticle::mass() const{
  return m_mass;
}

void ItoParticle::setMass(const Real& a_mass){
  m_mass = a_mass;
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

Real ItoParticle::conductivity() const {
  return m_mass*m_mobility;
}

Real ItoParticle::diffusivity() const {
  return m_mass*m_diffusion;
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

void ItoParticle::setEnergy(const Real a_averageEnergy){
  m_averageEnergy = a_averageEnergy;
}

Real& ItoParticle::energy(){
  return m_averageEnergy;
}

const Real& ItoParticle::energy() const{
  return m_averageEnergy;
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

#include <CD_NamespaceFooter.H>

#endif
