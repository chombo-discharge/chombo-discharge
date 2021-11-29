/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoSolverImplem.H
  @brief  Implementation of CD_ItoSolver.H
  @author Robert Marskar
  @date   April 2020
*/

#ifndef CD_ItoSolverImplem_H
#define CD_ItoSolverImplem_H

// Std includes
#include <chrono>

// Chombo includes
#include <EBAlias.H>
#include <PolyGeom.H>

// Our includes
#include <CD_ItoSolver.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

inline
Real ItoSolver::sign(const Real& a) const{
  return (a > 0) - (a < 0);
}

inline
RealVect ItoSolver::randomGaussian() {

  RealVect r = RealVect::Zero;
  for (int i = 0; i < SpaceDim; i++){
    r[i] = m_normalDistribution01(m_rng);
    r[i] = sign(r[i])*Min(Abs(r[i]), m_normalDistributionTruncation);
  }

  return r;
}

inline
RealVect ItoSolver::randomDirection(){
  constexpr Real EPS = 1.E-8;
#if CH_SPACEDIM==2
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = m_uniformDistribution11(m_rng);
    x2 = m_uniformDistribution11(m_rng);
    r  = x1*x1 + x2*x2;
  }

  return RealVect(x1,x2)/sqrt(r);
#elif CH_SPACEDIM==3
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = m_uniformDistribution11(m_rng);
    x2 = m_uniformDistribution11(m_rng);
    r  = x1*x1 + x2*x2;
  }

  const Real x = 2*x1*sqrt(1-r);
  const Real y = 2*x2*sqrt(1-r);
  const Real z = 1 - 2*r;

  return RealVect(x,y,z);
#endif
}

inline
RealVect ItoSolver::randomPosition(const RealVect a_pos,
				   const RealVect a_lo,
				   const RealVect a_hi,
				   const RealVect a_bndryCentroid,
				   const RealVect a_bndryNormal,
				   const Real     a_dx,
				   const Real     a_kappa) const {

  RealVect pos;
  if(a_kappa < 1.0){ // Rejection sampling. 
    pos = this->randomPosition(a_lo, a_hi, a_bndryCentroid, a_bndryNormal);
  }
  else{ // Regular cell. Get a position. 
    pos = this->randomPosition(a_lo, a_hi);
  }

  pos = a_pos + pos*a_dx;

  return pos;
}

inline
RealVect ItoSolver::randomPosition(const RealVect a_lo,
				   const RealVect a_hi,
				   const RealVect a_bndryCentroid,
				   const RealVect a_bndryNormal) const {
  RealVect pos = this->randomPosition(a_lo, a_hi);
  bool valid   = PolyGeom::dot(pos-a_bndryCentroid, a_bndryNormal) >= 0.0;

  while(!valid){
    pos    = this->randomPosition(a_lo, a_hi);
    valid = PolyGeom::dot(pos-a_bndryCentroid, a_bndryNormal) >= 0.0;
  }

  return pos;
}

inline
RealVect ItoSolver::randomPosition(const RealVect a_lo, const RealVect a_hi) const {

  RealVect pos = RealVect::Unit;

  for (int dir = 0; dir < SpaceDim; dir++){
    pos[dir] = a_lo[dir] + 0.5*(1.0 + m_uniformDistribution11(m_rng))*(a_hi[dir] - a_lo[dir]);
  }

  return pos;
}

template <class T>
void ItoSolver::depositParticles(EBAMRCellData& a_phi, const ParticleContainer<T>& a_particles){
  CH_TIME("ItoSolver::depositParticles");
  if(m_verbosity > 5){
    pout() << m_name + "::depositParticles" << endl;
  }

  this->depositParticles(a_phi, a_particles, m_deposition);
}

template <class T>
void ItoSolver::depositParticles(EBAMRCellData&               a_phi,
				 const ParticleContainer<T>& a_particles,
				 const DepositionType  a_deposition){
  CH_TIME("ItoSolver::depositParticles");
  if(m_verbosity > 5){
    pout() << m_name + "::depositParticles" << endl;
  }
           
  this->depositKappaConservative(a_phi, a_particles, a_deposition); // a_phi contains only weights, i.e. not divided by kappa
  if(m_useRedistribution){
    this->depositNonConservative(m_depositionNC, a_phi);              // Compute m_depositionNC = sum(kappa*Wc)/sum(kappa)
    this->depositHybrid(a_phi, m_massDiff, m_depositionNC);           // Compute hybrid deposition, including mass differnce
    this->incrementRedist(m_massDiff);                                  // Increment level redistribution register

    // Do the redistribution magic
    const bool ebcf = m_amr->getEbCf();
    if(ebcf){ // Mucho stuff to do here...
      this->coarseFineIncrement(m_massDiff);       // Compute C2F, F2C, and C2C mass transfers
      this->levelRedist(a_phi);           // Level redistribution. Weights is a dummy parameter
      this->coarseFineRedistribution(a_phi);     // Do the coarse-fine redistribution
    }
    else{ // Very simple, redistribute this level.
      this->levelRedist(a_phi);
    }
  }

  // Average down and interpolate
  m_amr->averageDown(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);
}

template <class P>
void ItoSolver::depositKappaConservative(EBAMRCellData&               a_phi,
					   const ParticleContainer<P>& a_particles,
					   const DepositionType  a_deposition){
  CH_TIME("ItoSolver::depositKappaConservative");
  if(m_verbosity > 5){
    pout() << m_name + "::depositKappaConservative" << endl;
  }

  CH_assert(m_haloBuffer == 0 || m_pvrBuffer == 0);

  if(m_haloBuffer > 0 && m_pvrBuffer > 0){
    MayDay::Error("ItoSolver::depositKappaConservative - conflicting inputs, halo/pvr buffers both > 0 but one of them must be 0");
  }  



  // Figure out how to handle the coarse-fine deposition.
  CoarseFineDeposition coarseFineDeposition;  
  if(m_pvrBuffer > 0){
    coarseFineDeposition = CoarseFineDeposition::PVR;
  }
  else if(m_haloBuffer > 0){
    if(m_forceHaloNGP){
      coarseFineDeposition = CoarseFineDeposition::HaloNGP;
    }
    else{
      coarseFineDeposition = CoarseFineDeposition::Halo;
    }
  }
  else{
    MayDay::Error("ItoSolver::depositKappaConservative -- logic bust");
  }


  // Now do the deposition. Recall that when we deposit with "halos", we need to fetch the subset of coarse-level particles that surround the
  // refinement boundary. 
  switch(coarseFineDeposition){
  case CoarseFineDeposition::PVR:
    m_amr->depositParticles<P, &P::mass>(a_phi, m_realm, m_phase, a_particles, a_deposition, coarseFineDeposition, m_forceIrregDepositionNGP);
    break;
  case CoarseFineDeposition::Halo:
    {
      const AMRMask& mask = m_amr->getMask(s_particle_halo, m_haloBuffer, m_realm);
      a_particles.copyMaskParticles(mask);
  
      m_amr->depositParticles<P, &P::mass>(a_phi, m_realm, m_phase, a_particles, a_deposition, CoarseFineDeposition::Halo, m_forceIrregDepositionNGP);

      a_particles.clearMaskParticles();
      
      break;
    }
  case CoarseFineDeposition::HaloNGP:
    {
      const AMRMask& mask = m_amr->getMask(s_particle_halo, m_haloBuffer, m_realm);

      a_particles.copyMaskParticles(mask);
      a_particles.copyNonMaskParticles(mask);

      m_amr->depositParticles<P, &P::mass>(a_phi, m_realm, m_phase, a_particles, a_deposition, CoarseFineDeposition::HaloNGP, m_forceIrregDepositionNGP);  
  
      a_particles.clearMaskParticles();
      a_particles.clearNonMaskParticles();

      break;
    }
  default:
    MayDay::Error("ItoSolverImplem.H in function ItoSolver::depositKappaConservative -- logic bust!");
    break;
  }
}

#include <CD_NamespaceFooter.H>

#endif
