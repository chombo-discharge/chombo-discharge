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

// Our includes
#include <CD_ItoSolver.H>
#include <CD_DataOps.H>
#include <CD_EBParticleMesh.H>
#include <CD_EBCoarseFineParticleMesh.H>
#include <CD_NamespaceHeader.H>

inline
Real ItoSolver::sign(const Real& a) const{
  return (a > 0) - (a < 0);
}

inline
RealVect ItoSolver::randomGaussian() {

  RealVect r = RealVect::Zero;
  for (int i = 0; i < SpaceDim; i++){
    r[i] = m_gauss01(m_rng);
    r[i] = sign(r[i])*Min(Abs(r[i]), m_normal_max);
  }

  return r;
}

inline
RealVect ItoSolver::randomDirection(){
  const Real EPS = 1.E-8;
#if CH_SPACEDIM==2
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = m_udist11(m_rng);
    x2 = m_udist11(m_rng);
    r  = x1*x1 + x2*x2;
  }

  return RealVect(x1,x2)/sqrt(r);
#elif CH_SPACEDIM==3
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = m_udist11(m_rng);
    x2 = m_udist11(m_rng);
    r  = x1*x1 + x2*x2;
  }

  const Real x = 2*x1*sqrt(1-r);
  const Real y = 2*x2*sqrt(1-r);
  const Real z = 1 - 2*r;

  return RealVect(x,y,z);
#endif
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
  if(m_redistribute){
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

template <class T>
void ItoSolver::depositKappaConservative(EBAMRCellData&               a_phi,
					   const ParticleContainer<T>& a_particles,
					   const DepositionType  a_deposition){
  CH_TIME("ItoSolver::depositKappaConservative");
  if(m_verbosity > 5){
    pout() << m_name + "::depositKappaConservative" << endl;
  }

  if(m_halo_buffer > 0 && m_pvr_buffer > 0){
    MayDay::Abort("ItoSolver::depositKappaConservative - conflicting inputs, halo/pvr buffers both > 0");
  }

  if(m_pvr_buffer > 0){
    this->depositKappaConservativeWithPVR(a_phi, a_particles, a_deposition);
  }
  else{
    if(m_ngp_halo || a_deposition == DepositionType::NGP){
      this->depositKappaConservativeNoPVR_NGP(a_phi, a_particles, a_deposition);
    }
    else{
      this->depositKappaConservativeNoPVR_native(a_phi, a_particles, a_deposition);
    }
  }
}

template <class T>
void ItoSolver::depositKappaConservativeWithPVR(EBAMRCellData&               a_phi,
						  const ParticleContainer<T>& a_particles,
						  const DepositionType  a_deposition){
  CH_TIME("ItoSolver::depositKappaConservativeWithPVR");
  if(m_verbosity > 5){
    pout() << m_name + "::depositKappaConservativeWithPVR" << endl;
  }

  m_amr->depositParticles<T, &T::mass>(a_phi, m_realm, m_phase, a_particles, a_deposition, CoarseFineDeposition::PVR, m_irreg_ngp_deposition);  
}

template <class T>
void ItoSolver::depositKappaConservativeNoPVR_native(EBAMRCellData&               a_phi,
						       const ParticleContainer<T>& a_particles,
						       const DepositionType  a_deposition){
  CH_TIME("ItoSolver::depositKappaConservativeNoPVR_native");
  if(m_verbosity > 5){
    pout() << m_name + "::depositKappaConservativeNoPVR_native" << endl;
  }

  const AMRMask& mask = m_amr->getMask(s_particle_halo, m_halo_buffer, m_realm);
  a_particles.copyMaskParticles(mask);
  
  m_amr->depositParticles<T, &T::mass>(a_phi, m_realm, m_phase, a_particles, a_deposition, CoarseFineDeposition::Halo, m_irreg_ngp_deposition);

  a_particles.clearMaskParticles();  
}

template <class T>
void ItoSolver::depositKappaConservativeNoPVR_NGP(EBAMRCellData&               a_phi,
						    const ParticleContainer<T>& a_particles,
						    const DepositionType  a_deposition){
  CH_TIME("ItoSolver::depositKappaConservativeNoPVR_NGP");
  if(m_verbosity > 5){
    pout() << m_name + "::depositKappaConservativeNoPVR_NGP" << endl;
  }

  const AMRMask& mask = m_amr->getMask(s_particle_halo, m_halo_buffer, m_realm);

  a_particles.copyMaskParticles(mask);
  a_particles.copyNonMaskParticles(mask);

  m_amr->depositParticles<T, &T::mass>(a_phi, m_realm, m_phase, a_particles, a_deposition, CoarseFineDeposition::HaloNGP, m_irreg_ngp_deposition);  
  
  a_particles.clearMaskParticles();
  a_particles.clearNonMaskParticles();  
}

#include <CD_NamespaceFooter.H>

#endif
