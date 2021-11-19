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

  MayDay::Warning("EBAMRParticleMesh::depositPVR -- not implemented");  
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
