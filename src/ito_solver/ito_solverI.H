/*!
  @file   ito_solverI.H
  @brief  Templated functions for ito_solver
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _ITO_SOLVERI_H_
#define _ITO_SOLVERI_H_

#include "ito_solver.H"
#include "simple_ito_particle.H"
#include "ito_solver.H"
#include "data_ops.H"
#include "EBParticleInterp.H"
#include "units.H"
#include "EBGhostCloud.H"
#include "ito_layout.H"


#include "particle_ops.H"

#include <EBArith.H>
#include <ParmParse.H>
#include <EBAlias.H>
#include <BaseEBCellFactory.H>
#include <ParticleIO.H>

#include <chrono>

template <class T>
void ito_solver::deposit_particles(EBAMRCellData& a_state, const particle_container<T>& a_particles){
  CH_TIME("ito_solver::deposit_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_particles" << endl;
  }

  this->deposit_particles(a_state, a_particles, m_deposition);
}

template <class T>
void ito_solver::deposit_particles(EBAMRCellData&               a_state,
				   const particle_container<T>& a_particles,
				   const DepositionType::Which  a_deposition){
  CH_TIME("ito_solver::deposit_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_particles" << endl;
  }
           
  this->deposit_kappaConservative(a_state, a_particles, a_deposition); // a_state contains only weights, i.e. not divided by kappa
  if(m_redistribute){
    this->deposit_nonConservative(m_depositionNC, a_state);              // Compute m_depositionNC = sum(kappa*Wc)/sum(kappa)
    this->deposit_hybrid(a_state, m_massDiff, m_depositionNC);           // Compute hybrid deposition, including mass differnce
    this->increment_redist(m_massDiff);                                  // Increment level redistribution register

    // Do the redistribution magic
    const bool ebcf = m_amr->get_ebcf();
    if(ebcf){ // Mucho stuff to do here...
      this->coarse_fine_increment(m_massDiff);       // Compute C2F, F2C, and C2C mass transfers
      this->level_redistribution(a_state);           // Level redistribution. Weights is a dummy parameter
      this->coarse_fine_redistribution(a_state);     // Do the coarse-fine redistribution
    }
    else{ // Very simple, redistribute this level.
      this->level_redistribution(a_state);
    }
  }

  // Average down and interpolate
  m_amr->average_down(a_state, m_realm, m_phase);
  m_amr->interp_ghost(a_state, m_realm, m_phase);
}

template <class T>
void ito_solver::deposit_kappaConservative(EBAMRCellData&               a_state,
					   const particle_container<T>& a_particles,
					   const DepositionType::Which  a_deposition){
  CH_TIME("ito_solver::deposit_kappaConservative");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_kappaConservative" << endl;
  }

  if(m_pvr_buffer > 0){
    this->deposit_kappaConservativeWithPVR(a_state, a_particles, a_deposition);
  }
  else{
    if(m_ngp_halo || a_deposition == DepositionType::NGP){
      this->deposit_kappaConservativeNoPVR_NGP(a_state, a_particles, a_deposition);
    }
    else{
      this->deposit_kappaConservativeNoPVR_native(a_state, a_particles, a_deposition);
    }

  }
}

template <class T>
void ito_solver::deposit_kappaConservativeWithPVR(EBAMRCellData&               a_state,
						  const particle_container<T>& a_particles,
						  const DepositionType::Which  a_deposition){
  CH_TIME("ito_solver::deposit_kappaConservative");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_kappaConservative" << endl;
  }

  const int comp = 0;
  const Interval interv(comp, comp);

  const RealVect origin  = m_amr->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  data_ops::set_value(a_state,    0.0);
  data_ops::set_value(m_scratch,  0.0);

  const AMRParticles<T>& particles = a_particles.get_particles();
  const AMRParticles<T>& halo_particles = a_particles.get_halo_particles();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const ProblemDomain& dom     = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[lvl];
    const RefCountedPtr<EBLevelGrid>& eblg = m_amr->get_eblg(m_realm, m_phase)[lvl];

    const bool has_coar = (lvl > 0);
    const bool has_fine = (lvl < finest_level);

    // 1. If we have a coarser level whose cloud extends beneath this level, interpolate that result here first. 
    if(has_coar){
      RefCountedPtr<EBMGInterp>& interp = m_amr->get_eb_mg_interp(m_realm, m_phase)[lvl];
      interp->pwcInterp(*a_state[lvl], *m_scratch[lvl-1], interv);
    }
    
    // 2. Deposit this levels particles. Note that this will deposit into ghost cells, which must later
    //    be added to neighboring patches
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      EBParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin, m_redistribute);
      interp.deposit((*particles[lvl])[dit()].listItems(), (*a_state[lvl])[dit()].getFArrayBox(), m_deposition);
    }

    // This code adds contributions from ghost cells into the valid region
    const RefCountedPtr<Copier>& reversecopier = m_amr->get_reverse_copier(m_realm, m_phase)[lvl];
    LDaddOp<FArrayBox> addOp;
    LevelData<FArrayBox> aliasFAB;
    aliasEB(aliasFAB, *a_state[lvl]);
    aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);

    // 3. If we have a finer level, copy contributions from this level to the temporary holder that is used for
    //    interpolation of "hanging clouds"
    if(has_fine){
      a_state[lvl]->localCopyTo(*m_scratch[lvl]);
    }
  }
}

template <class T>
void ito_solver::deposit_kappaConservativeNoPVR_native(EBAMRCellData&               a_state,
						  const particle_container<T>& a_particles,
						  const DepositionType::Which  a_deposition){
  CH_TIME("ito_solver::deposit_kappaConservativeNoPVR_native");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_kappaConservativeNoPVR_native" << endl;
  }

  const int comp = 0;
  const Interval interv(comp, comp);

  const RealVect origin  = m_amr->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  data_ops::set_value(a_state,    0.0);

  a_particles.copy_halo_particles();
    
  const AMRParticles<T>& particles          = a_particles.get_particles();
  const AMRParticles<T>& halo_particles     = a_particles.get_halo_particles();
  const AMRParticles<T>& non_halo_particles = a_particles.get_non_halo_particles();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const ProblemDomain& dom     = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[lvl];

    const bool has_coar = (lvl > 0);
    const bool has_fine = (lvl < finest_level);

    // 1. Deposit this levels particles. Note that this will deposit into ghost cells. 
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      EBParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin, m_redistribute);
      interp.deposit((*particles[lvl])[dit()].listItems(), (*a_state[lvl])[dit()].getFArrayBox(), m_deposition);
    }

    // 2. Exchange ghost cells on this level. 
    const RefCountedPtr<Copier>& reversecopier = m_amr->get_reverse_copier(m_realm, m_phase)[lvl];
    LDaddOp<FArrayBox> addOp;
    LevelData<FArrayBox> aliasFAB;
    aliasEB(aliasFAB, *a_state[lvl]);
    aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);

    // 3. Deposition into ghost cells across the fine boundary should end up on the coarse level. Add that
    //    mass to the coarse level right now. 
    if(has_coar) {
      EBGhostCloud& ghostcloud = *(m_amr->get_ghostcloud(m_realm, m_phase)[lvl]);
      ghostcloud.addFineGhostsToCoarse(*a_state[lvl-1], *a_state[lvl]);
    }

    // 4. The particles on the coarse side of the refinement boundary should deposit onto this level.
    //    The halo particles have been copied, and we run a special deposition procedure here to ensure
    //    that the mass ends up on the fine level. 
    if(has_coar){
      const int refRat = m_amr->get_ref_rat()[lvl-1];
      EBGhostCloud& ghostcloud         = *(m_amr->get_ghostcloud(m_realm, m_phase)[lvl]);
      BoxLayoutData<FArrayBox>& buffer = ghostcloud.getFiCoBuffer();
      const EBLevelGrid& eblgFiCo      = ghostcloud.getEblgFiCo();
      
      for (DataIterator dit = buffer.dataIterator(); dit.ok(); ++dit){
	const Box box = buffer.box(dit());
	const EBISBox& ebisbox = eblgFiCo.getEBISL()[dit()];

	buffer[dit()].setVal(0.0);
	EBParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin, m_redistribute);

	if(refRat == 2){
	  interp.deposit2((*halo_particles[lvl-1])[dit()].listItems(), buffer[dit()], m_deposition);
	}
	else if (refRat == 4){
	  interp.deposit4((*halo_particles[lvl-1])[dit()].listItems(), buffer[dit()], m_deposition);
	}
	else{
	  MayDay::Abort("ito_solverI::deposit_kappaConservativeNoPVR - logic bust");
	}
      }

      // Add the result of halo particle depositions to the fine level. 
      ghostcloud.addFiCoDataToFine(*a_state[lvl], buffer);
    }
  }

  // Clear halo and non-halo particles
  a_particles.clear_halo_particles();
}

template <class T>
void ito_solver::deposit_kappaConservativeNoPVR_NGP(EBAMRCellData&               a_state,
						  const particle_container<T>& a_particles,
						  const DepositionType::Which  a_deposition){
  CH_TIME("ito_solver::deposit_kappaConservativeNoPVR_NGP");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_kappaConservativeNoPVR_NGP" << endl;
  }

  const int comp = 0;
  const Interval interv(comp, comp);

  const RealVect origin  = m_amr->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  data_ops::set_value(a_state,    0.0);

  a_particles.copy_halo_particles();
  a_particles.copy_non_halo_particles();
    
  const AMRParticles<T>& particles          = a_particles.get_particles();
  const AMRParticles<T>& halo_particles     = a_particles.get_halo_particles();
  const AMRParticles<T>& non_halo_particles = a_particles.get_non_halo_particles(); 

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const ProblemDomain& dom     = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[lvl];

    const bool has_coar = (lvl > 0);
    const bool has_fine = (lvl < finest_level);

    // 1. Deposit "internal particles" with the desired deposition scheme. 
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      EBParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin, m_redistribute);
      interp.deposit((*non_halo_particles[lvl])[dit()].listItems(), (*a_state[lvl])[dit()].getFArrayBox(), m_deposition);
    }

    // 2. Deposition into ghost cells across the fine boundary should end up on the coarse level. 
    if(has_coar) {
      EBGhostCloud& ghostcloud = *(m_amr->get_ghostcloud(m_realm, m_phase)[lvl]);
      ghostcloud.addFineGhostsToCoarse(*a_state[lvl-1], *a_state[lvl]);
    }

    // 3. Exchange ghost cells on this level. 
    const RefCountedPtr<Copier>& reversecopier = m_amr->get_reverse_copier(m_realm, m_phase)[lvl];
    LDaddOp<FArrayBox> addOp;
    LevelData<FArrayBox> aliasFAB;
    aliasEB(aliasFAB, *a_state[lvl]);
    aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);

    // 4. Deposit particles in the halo with an NGP method
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      EBParticleInterp interp(box, ebisbox, dx*RealVect::Unit, origin, m_redistribute);
      interp.deposit((*halo_particles[lvl])[dit()].listItems(), (*a_state[lvl])[dit()].getFArrayBox(), DepositionType::NGP);
    }    
  }

  // Clear halo and non-halo particles
  a_particles.clear_halo_particles();
  a_particles.clear_non_halo_particles();
}


#endif
