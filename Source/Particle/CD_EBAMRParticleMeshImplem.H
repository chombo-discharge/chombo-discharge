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

  MayDay::Warning("EBAMRParticleMesh::depositPVR -- not implemented");

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

    // 3. The mass on the part of the coarse level underneath this grid should be added to this level. Fortunately,
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

  MayDay::Warning("EBAMRParticleMesh::depositHalo -- not implemented");
}

#include <CD_NamespaceFooter.H>

#endif
