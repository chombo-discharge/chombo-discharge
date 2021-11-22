/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRParticleMeshImplem.H
  @brief  Implementation of CD_EBAMRParticleMesh.H
  @author Robert Marskar
*/

#ifndef CD_EBAMRParticleMeshImplem_H
#define CD_EBAMRParticleMeshImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBAMRParticleMesh.H>
#include <CD_EBAddOp.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

template <class T, const Real&(T::*particleScalarField)() const>
void EBAMRParticleMesh::deposit(EBAMRCellData&              a_meshData,
				const ParticleContainer<T>& a_particles,
				const DepositionType        a_depositionType,
				const CoarseFineDeposition  a_coarseFineDeposition,
				const bool                  a_forceIrregNGP) {
  CH_TIME("EBAMRParticleMesh::deposit");

  switch(a_coarseFineDeposition){
  case CoarseFineDeposition::PVR:
    EBAMRParticleMesh::depositPVR<T, particleScalarField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);
    break;
  case CoarseFineDeposition::Halo:
    EBAMRParticleMesh::depositHalo<T, particleScalarField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);    
    break;    
  default:
    MayDay::Error("EBAMRParticleMesh::deposit - logic bust");
    break;
  }
}

template <class T, const Real&(T::*particleScalarField)() const>
void EBAMRParticleMesh::depositPVR(EBAMRCellData&              a_meshData,
				   const ParticleContainer<T>& a_particles,
				   const DepositionType        a_depositionType,
				   const bool                  a_forceIrregNGP){
  CH_TIME("EBAMRParticleMesh::depositPVR");
  
  CH_assert(a_meshData[0]->nComp() == 1);

  // TLDR: On each level the particles are deposited as usual. This routine assumes that none of the fine-level particles
  //       have deposition clouds that hang over the refinement boundary, so there is no fine-to-coarse mass transfer from
  //       fine-level particles. On the coarse level the particles might have deposition clouds that extend beneath the
  //       fine grid. That mass needs to find it's way into the fine level. The way we handle this is that we do piecewise
  //       constant interpolation of the coarse-grid data to the fine grid. 

  DataOps::setValue(a_meshData, 0.0);

  const Interval interv(0,0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    
    const ParticleData<T>&   particles = *a_particles[lvl];
    const DisjointBoxLayout& dbl       = m_eblgs[lvl]->getDBL   ();
    const EBISLayout&        ebisl     = m_eblgs[lvl]->getEBISL ();
    const ProblemDomain&     domain    = m_eblgs[lvl]->getDomain();

    const bool hasCoar = lvl > 0;
    
    // 1. Deposit particles on this level.
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box      cellBox = dbl  [dit()];
      const EBISBox& ebisBox = ebisl[dit()];

      EBCellFAB&        meshData     = (*a_meshData [lvl])[dit()];
      const ListBox<T>& boxParticles = (*a_particles[lvl])[dit()];

      // Note: EbParticleInterp construction is cheap.
      EbParticleInterp interp(cellBox, ebisBox, m_dx[lvl]*RealVect::Unit, m_probLo);

      interp.deposit<T, particleScalarField> (boxParticles.listItems(), meshData, a_depositionType, a_forceIrregNGP);
    }


    // 2. Exchange ghost data on this level. After this, all the mass should be on the current level. 
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. If the particles deposited over the coarse-fine boundary, add the mass to the coarse level.
    if(hasCoar){
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl-1], *a_meshData[lvl]);
    }

    // 4. The mass on the part of the coarse level underneath this grid should be added to this level. Fortunately,
    //    EBCoarseFineParticleMesh has buffers and functions that handle that type of interpolation. 
    if(hasCoar){
      m_coarseFinePM[lvl]->addInvalidCoarseToFine(*a_meshData[lvl-1], *a_meshData[lvl]);
    }
  }
}

template <class T, const Real&(T::*particleScalarField)() const>
void EBAMRParticleMesh::depositHalo(EBAMRCellData&              a_meshData,
				    const ParticleContainer<T>& a_particles,
				    const DepositionType        a_depositionType,
				    const bool                  a_forceIrregNGP){
  CH_TIME("EBAMRParticleMesh::depositHalo");

  CH_assert(a_meshData[0]->nComp() == 1);

  // TLDR: On each level the particles are deposited as usual. This routine permits particles to live right next to the refinement boundary
  //       so both the coarse-level and fine-level particles have clouds that hang into the fine/coarse level. The way we handle this is as follows:
  //       The fine-level particles will have deposited into ghost cells that lie around the refinement boundary. This mass is put on the coarse grid.
  //       On the coarse level the particles have clouds that extend beneath the fine grid. Instead of interpolating this mass to the fine grid we
  //       run a special deposition procedure which allows us to deposit directly on the fine grid using 2x or 4x the particle width (depending on the
  //       refinement ratio).

  DataOps::setValue(a_meshData, 0.0);
  
  const Interval interv(0,0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const ParticleData<T>&   particles = *a_particles[lvl];
    const DisjointBoxLayout& dbl       = m_eblgs[lvl]->getDBL   ();
    const EBISLayout&        ebisl     = m_eblgs[lvl]->getEBISL ();
    const ProblemDomain&     domain    = m_eblgs[lvl]->getDomain();

    const bool hasCoar = (lvl > 0);

    // 1. Deposit particles on this level.
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box      cellBox = dbl  [dit()];
      const EBISBox& ebisBox = ebisl[dit()];

      EBCellFAB&        meshData     = (*a_meshData [lvl])[dit()];
      const ListBox<T>& boxParticles = (*a_particles[lvl])[dit()];

      // Note: EbParticleInterp construction is cheap.
      EbParticleInterp interp(cellBox, ebisBox, m_dx[lvl]*RealVect::Unit, m_probLo);

      interp.deposit<T, particleScalarField> (boxParticles.listItems(), meshData, a_depositionType, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level. 
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary should end up to the coarse level. Add that mass right now.
    if(hasCoar){
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl-1], *a_meshData[lvl]);
    }

    // 4. The particles on the other side of the refinement bounadry should deposit into this level.
    if(hasCoar){
      const ParticleData<T>&   coarHaloParticles = *a_particles.getMaskParticles()[lvl-1];  // These are the particles that live "around" the refinement boundary.
      const int                refRat            =  m_refRat[lvl-1];                        // Refinement between coarse level and this level.
      const EBLevelGrid&       eblgFiCo          =  m_coarseFinePM[lvl]->getEblgFiCo();     // This is the refinement of the coarse grid.
      const DisjointBoxLayout& dblFiCo           =  eblgFiCo.getDBL();                      // Refinement of coarse grid
      const ProblemDomain&     domainFiCo        =  eblgFiCo.getDomain();                   // Refinement of coarse domain
      const EBISLayout&        ebislFiCo         =  eblgFiCo.getEBISL();                    // Refinement of coarse EBISL
      
      LevelData<EBCellFAB>&  bufferFiCo          = m_coarseFinePM[lvl]->getFiCoBuffer();   // Data that is a refinement of the coarse level.

      for (DataIterator dit(dbl); dit.ok(); ++dit){
	EBCellFAB&      dataFiCo       = bufferFiCo       [dit()];
	const ListBox<T> haloParticles = coarHaloParticles[dit()];
	const Box        cellBox       = dblFiCo          [dit()];
	const EBISBox&   ebisBox       = ebislFiCo        [dit()];
	
	dataFiCo.setVal(0.0);

	EbParticleInterp interp(cellBox, ebisBox, m_dx[lvl-1]*RealVect::Unit, m_probLo);

	switch(refRat){
	case 2:
	  interp.deposit2<T, particleScalarField> (haloParticles.listItems(), dataFiCo, a_depositionType, a_forceIrregNGP);
	  break;
	case 4:
	  interp.deposit4<T, particleScalarField> (haloParticles.listItems(), dataFiCo, a_depositionType, a_forceIrregNGP);
	  break;
	default:
	  MayDay::Error("CD_EBAMRParticleMeshImplem.H - logic bust in EBAMRParticleMesh::depositHalo");
	  break;
	}
      }

      // Add the result of the buffer deposition to this level.
      m_coarseFinePM[lvl]->addFiCoDataToFine(*a_meshData[lvl-1], bufferFiCo);
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
