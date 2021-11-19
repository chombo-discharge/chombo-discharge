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
#include <CD_EbParticleInterp.H>
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
  
  const int comp = 0;
  const Interval interv(comp, comp);

  const RealVect origin  = m_amr->getProbLo();
  const int finest_level = m_amr->getFinestLevel();

  DataOps::setValue(a_phi,    0.0);
  DataOps::setValue(m_scratch,  0.0);

  const AMRParticles<T>& particles      = a_particles.getParticles();
  const AMRParticles<T>& halo_particles = a_particles.getMaskParticles();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->getDx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& dom     = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const RefCountedPtr<EBLevelGrid>& eblg = m_amr->getEBLevelGrid(m_realm, m_phase)[lvl];

    const bool has_coar = (lvl > 0);
    const bool has_fine = (lvl < finest_level);

    // 1. If we have a coarser level whose cloud extends beneath this level, interpolate that result here first. 
    if(has_coar){
      RefCountedPtr<EBMGInterp>& interp = m_amr->getEBMGInterp(m_realm, m_phase)[lvl];
      interp->pwcInterp(*a_phi[lvl], *m_scratch[lvl-1], interv);
    }

    // 2. Deposit this levels particles. Note that this will deposit into ghost cells, which must later
    //    be added to neighboring patches
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      EbParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin);
      interp.deposit<ItoParticle, &ItoParticle::mass>((*particles[lvl])[dit()].listItems(), (*a_phi[lvl])[dit()], m_deposition, m_irreg_ngp_deposition);      
    }

    // This code adds contributions from ghost cells into the valid region
    const RefCountedPtr<Copier>& reversecopier = m_amr->getReverseCopier(m_realm, m_phase)[lvl];
    LDaddOp<FArrayBox> addOp;
    LevelData<FArrayBox> aliasFAB;
    aliasEB(aliasFAB, *a_phi[lvl]);
    aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);

    // 3. If we have a finer level, copy contributions from this level to the temporary holder that is used for
    //    interpolation of "hanging clouds"
    if(has_fine){
      a_phi[lvl]->localCopyTo(*m_scratch[lvl]);
    }
  }
}

template <class T>
void ItoSolver::depositKappaConservativeNoPVR_native(EBAMRCellData&               a_phi,
						       const ParticleContainer<T>& a_particles,
						       const DepositionType  a_deposition){
  CH_TIME("ItoSolver::depositKappaConservativeNoPVR_native");
  if(m_verbosity > 5){
    pout() << m_name + "::depositKappaConservativeNoPVR_native" << endl;
  }

  const int comp = 0;
  const Interval interv(comp, comp);

  const RealVect origin  = m_amr->getProbLo();
  const int finest_level = m_amr->getFinestLevel();

  DataOps::setValue(a_phi,    0.0);

  //  const EBAMRBool& mask = m_amr->getMask(s_particle_halo, m_halo_buffer, m_realm);
  const AMRMask& mask = m_amr->getMask(s_particle_halo, m_halo_buffer, m_realm);
  a_particles.copyMaskParticles(mask);
    
  const AMRParticles<T>& particles          = a_particles.getParticles();
  const AMRParticles<T>& halo_particles     = a_particles.getMaskParticles();
  const AMRParticles<T>& non_halo_particles = a_particles.getNonMaskParticles();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->getDx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& dom     = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    const bool has_coar = (lvl > 0);
    const bool has_fine = (lvl < finest_level);

    // 1. Deposit this levels particles. Note that this will deposit into ghost cells. 
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      EbParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin);
      interp.deposit<ItoParticle, &ItoParticle::mass>((*particles[lvl])[dit()].listItems(), (*a_phi[lvl])[dit()], m_deposition, m_irreg_ngp_deposition);      
    }

    // 2. Exchange ghost cells on this level. 
    const RefCountedPtr<Copier>& reversecopier = m_amr->getReverseCopier(m_realm, m_phase)[lvl];
    LDaddOp<FArrayBox> addOp;
    LevelData<FArrayBox> aliasFAB;
    aliasEB(aliasFAB, *a_phi[lvl]);
    aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);

    // 3. Deposition into ghost cells across the fine boundary should end up on the coarse level. Add that
    //    mass to the coarse level right now. 
    if(has_coar) {
      EBCoarseFineParticleMesh& ghostcloud = *(m_amr->getEBCoarseFineParticleMesh(m_realm, m_phase)[lvl]);
      ghostcloud.addFineGhostsToCoarse(*a_phi[lvl-1], *a_phi[lvl]);
    }

    // 4. The particles on the coarse side of the refinement boundary should deposit onto this level.
    //    The halo particles have been copied, and we run a special deposition procedure here to ensure
    //    that the mass ends up on the fine level. 
    if(has_coar){
      const int refRat = m_amr->getRefinementRatios()[lvl-1];
      EBCoarseFineParticleMesh& ghostcloud         = *(m_amr->getEBCoarseFineParticleMesh(m_realm, m_phase)[lvl]);
      LevelData<EBCellFAB>& buffer = ghostcloud.getFiCoBuffer();
      const EBLevelGrid& eblgFiCo      = ghostcloud.getEblgFiCo();
      
      for (DataIterator dit = buffer.dataIterator(); dit.ok(); ++dit){
	const Box box = buffer.box(dit());
	const EBISBox& ebisbox = eblgFiCo.getEBISL()[dit()];

	buffer[dit()].setVal(0.0);
	EbParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin);

	if(refRat == 2){
	  interp.deposit2<ItoParticle, &ItoParticle::mass>((*halo_particles[lvl-1])[dit()].listItems(), buffer[dit()], m_deposition, m_irreg_ngp_deposition);
	}
	else if (refRat == 4){
	  interp.deposit4<ItoParticle, &ItoParticle::mass>((*halo_particles[lvl-1])[dit()].listItems(), buffer[dit()], m_deposition, m_irreg_ngp_deposition);	  
	}
	else{
	  MayDay::Abort("ItoSolverI::depositKappaConservativeNoPVR - logic bust");
	}
      }

      // Add the result of halo particle depositions to the fine level. 
      ghostcloud.addFiCoDataToFine(*a_phi[lvl], buffer);
    }
  }

  // Clear halo and non-halo particles
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

  const int comp = 0;
  const Interval interv(comp, comp);

  const RealVect origin  = m_amr->getProbLo();
  const int finest_level = m_amr->getFinestLevel();

  DataOps::setValue(a_phi,    0.0);

  const AMRMask& mask = m_amr->getMask(s_particle_halo, m_halo_buffer, m_realm);

  a_particles.copyMaskParticles(mask);
  a_particles.copyNonMaskParticles(mask);
    
  const AMRParticles<T>& particles          = a_particles.getParticles();
  const AMRParticles<T>& halo_particles     = a_particles.getMaskParticles();
  const AMRParticles<T>& non_halo_particles = a_particles.getNonMaskParticles(); 

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->getDx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& dom     = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    const bool has_coar = (lvl > 0);
    const bool has_fine = (lvl < finest_level);

    // 1. Deposit "internal particles" with the desired deposition scheme. 
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      EbParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin);
      interp.deposit<ItoParticle, &ItoParticle::mass>((*non_halo_particles[lvl])[dit()].listItems(), (*a_phi[lvl])[dit()], m_deposition, m_irreg_ngp_deposition);    }

    // 2. Deposition into ghost cells across the fine boundary should end up on the coarse level. 
    if(has_coar) {
      EBCoarseFineParticleMesh& ghostcloud = *(m_amr->getEBCoarseFineParticleMesh(m_realm, m_phase)[lvl]);
      ghostcloud.addFineGhostsToCoarse(*a_phi[lvl-1], *a_phi[lvl]);
    }

    // 3. Exchange ghost cells on this level. 
    const RefCountedPtr<Copier>& reversecopier = m_amr->getReverseCopier(m_realm, m_phase)[lvl];
    LDaddOp<FArrayBox> addOp;
    LevelData<FArrayBox> aliasFAB;
    aliasEB(aliasFAB, *a_phi[lvl]);
    aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);

    // 4. Deposit particles in the halo with an NGP method
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      EbParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin);
      interp.deposit<ItoParticle, &ItoParticle::mass>((*halo_particles[lvl])[dit()].listItems(), (*a_phi[lvl])[dit()], DepositionType::NGP, m_irreg_ngp_deposition);      
    }    
  }

  // Clear halo and non-halo particles
  a_particles.clearMaskParticles();
  a_particles.clearNonMaskParticles();
}

#include <CD_NamespaceFooter.H>

#endif
