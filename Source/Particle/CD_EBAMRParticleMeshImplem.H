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
#include <EBAlias.H>
#include <LevelData.H>
#include <EBCellFactory.H>

// Our includes
#include <CD_EBAMRParticleMesh.H>
#include <CD_EBAddOp.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

template <class P, const Real& (P::*particleScalarField)() const>
void
EBAMRParticleMesh::deposit(EBAMRCellData&              a_meshData,
                           const ParticleContainer<P>& a_particles,
                           const DepositionType        a_depositionType,
                           const CoarseFineDeposition  a_coarseFineDeposition,
                           const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::deposit(scalar)");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::deposit(scalar)" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == 1);

  switch (a_coarseFineDeposition) {
  case CoarseFineDeposition::Interp: {
    this->depositInterp<P, particleScalarField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::Halo: {
    this->depositHalo<P, particleScalarField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::HaloNGP: {
    this->depositHaloNGP<P, particleScalarField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::Transition: {
    this->depositTransition<P, particleScalarField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  default: {
    MayDay::Error("EBAMRParticleMesh::deposit(scalar) - logic bust");

    break;
  }
  }
}

template <class P, Real (P::*particleScalarField)() const>
void
EBAMRParticleMesh::deposit(EBAMRCellData&              a_meshData,
                           const ParticleContainer<P>& a_particles,
                           const DepositionType        a_depositionType,
                           const CoarseFineDeposition  a_coarseFineDeposition,
                           const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::deposit(scalar)");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::deposit(scalar)" << endl;
  }

  // Note: This is just like the other version except that the function signature is Real P::particleScalarField() const
  CH_assert(a_meshData[0]->nComp() == 1);

  switch (a_coarseFineDeposition) {
  case CoarseFineDeposition::Interp: {
    this->depositInterp<P, particleScalarField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::Halo: {
    this->depositHalo<P, particleScalarField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::HaloNGP: {
    this->depositHaloNGP<P, particleScalarField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::Transition: {
    this->depositTransition<P, particleScalarField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  default: {
    MayDay::Error("EBAMRParticleMesh::deposit(scalar) - logic bust");

    break;
  }
  }
}

template <class P, const RealVect& (P::*particleVectorField)() const>
void
EBAMRParticleMesh::deposit(EBAMRCellData&              a_meshData,
                           const ParticleContainer<P>& a_particles,
                           const DepositionType        a_depositionType,
                           const CoarseFineDeposition  a_coarseFineDeposition,
                           const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::deposit(vector)");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::deposit(vector)" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == SpaceDim);

  switch (a_coarseFineDeposition) {
  case CoarseFineDeposition::Interp: {
    this->depositInterp<P, particleVectorField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::Halo: {
    this->depositHalo<P, particleVectorField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::HaloNGP: {
    this->depositHaloNGP<P, particleVectorField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::Transition: {
    this->depositTransition<P, particleVectorField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  default: {
    MayDay::Error("EBAMRParticleMesh::deposit(RealVect) - logic bust");

    break;
  }
  }
}

template <class P, RealVect (P::*particleVectorField)() const>
void
EBAMRParticleMesh::deposit(EBAMRCellData&              a_meshData,
                           const ParticleContainer<P>& a_particles,
                           const DepositionType        a_depositionType,
                           const CoarseFineDeposition  a_coarseFineDeposition,
                           const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::deposit(vector)");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::deposit(vector)" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == SpaceDim);

  switch (a_coarseFineDeposition) {
  case CoarseFineDeposition::Interp: {
    this->depositInterp<P, particleVectorField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::Halo: {
    this->depositHalo<P, particleVectorField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::HaloNGP: {
    this->depositHaloNGP<P, particleVectorField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::Transition: {
    this->depositTransition<P, particleVectorField>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  default: {
    MayDay::Error("EBAMRParticleMesh::deposit(RealVect) - logic bust");

    break;
  }
  }
}

template <class P, const Real& (P::*particleScalarField)() const>
void
EBAMRParticleMesh::depositInterp(EBAMRCellData&              a_meshData,
                                 const ParticleContainer<P>& a_particles,
                                 const DepositionType        a_depositionType,
                                 const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositInterp");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositInterp" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == 1);

  // TLDR: On each level the particles are deposited as usual. This routine assumes that none of the fine-level particles
  //       have deposition clouds that hang over the refinement boundary, so there is no fine-to-coarse mass transfer from
  //       fine-level particles. On the coarse level the particles might have deposition clouds that extend beneath the
  //       fine grid. That mass needs to find it's way into the fine level. The way we handle this is that we do piecewise
  //       constant interpolation of the coarse-grid data to the fine grid.

  DataOps::setValue(a_meshData, 0.0);

  const Interval interv(0, 0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl      = m_eblgs[lvl]->getDBL();
    const bool               hasCoar  = lvl > 0;
    const DataIterator&      dit      = dbl.dataIterator();
    const int                numBoxes = dit.size();

    // 1. Deposit particles on this level.
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&        meshData     = (*a_meshData[lvl])[din];
      const ListBox<P>& boxParticles = a_particles[lvl][din];

      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleScalarField>(boxParticles.listItems(), meshData, a_depositionType, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level. After this, all the mass should be on the current level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. If the particles deposited over the coarse-fine boundary, add the mass to the coarse level.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
    }

    // 4. The mass on the part of the coarse level underneath this grid should be added to this level. Fortunately,
    //    EBCoarseFineParticleMesh has buffers and functions that handle that type of interpolation.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addInvalidCoarseToFine(*a_meshData[lvl], *a_meshData[lvl - 1]);
    }
  }
}

template <class P, Real (P::*particleScalarField)() const>
void
EBAMRParticleMesh::depositInterp(EBAMRCellData&              a_meshData,
                                 const ParticleContainer<P>& a_particles,
                                 const DepositionType        a_depositionType,
                                 const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositInterp");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositInterp" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == 1);

  // TLDR: On each level the particles are deposited as usual. This routine assumes that none of the fine-level particles
  //       have deposition clouds that hang over the refinement boundary, so there is no fine-to-coarse mass transfer from
  //       fine-level particles. On the coarse level the particles might have deposition clouds that extend beneath the
  //       fine grid. That mass needs to find it's way into the fine level. The way we handle this is that we do piecewise
  //       constant interpolation of the coarse-grid data to the fine grid.

  DataOps::setValue(a_meshData, 0.0);

  const Interval interv(0, 0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl   = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit     = dbl.dataIterator();
    const bool               hasCoar = lvl > 0;

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&        meshData     = (*a_meshData[lvl])[din];
      const ListBox<P>& boxParticles = a_particles[lvl][din];

      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleScalarField>(boxParticles.listItems(), meshData, a_depositionType, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level. After this, all the mass should be on the current level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. If the particles deposited over the coarse-fine boundary, add the mass to the coarse level.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
    }

    // 4. The mass on the part of the coarse level underneath this grid should be added to this level. Fortunately,
    //    EBCoarseFineParticleMesh has buffers and functions that handle that type of interpolation.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addInvalidCoarseToFine(*a_meshData[lvl], *a_meshData[lvl - 1]);
    }
  }
}

template <class P, const Real& (P::*particleScalarField)() const>
void
EBAMRParticleMesh::depositHalo(EBAMRCellData&              a_meshData,
                               const ParticleContainer<P>& a_particles,
                               const DepositionType        a_depositionType,
                               const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositHalo");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositHalo" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == 1);

  // TLDR: On each level the particles are deposited as usual. This routine permits particles to live right next to the refinement boundary
  //       so both the coarse-level and fine-level particles have clouds that hang into the fine/coarse level. The way we handle this is as follows:
  //       The fine-level particles will have deposited into ghost cells that lie around the refinement boundary. This mass is put on the coarse grid.
  //       On the coarse level the particles have clouds that extend beneath the fine grid. Instead of interpolating this mass to the fine grid we
  //       run a special deposition procedure which allows us to deposit directly on the fine grid using 2x or 4x the particle width (depending on the
  //       refinement ratio).

  DataOps::setValue(a_meshData, 0.0);

  // Copy the required particles to a masked particle data holder -- this takes the first strip of cells on the coarse side of the refinemnt boundary.
  constexpr int coarseMaskWidth = 1;

  a_particles.copyMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

  const Interval interv(0, 0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl   = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit     = dbl.dataIterator();
    const bool               hasCoar = (lvl > 0);

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&        meshData     = (*a_meshData[lvl])[din];
      const ListBox<P>& boxParticles = a_particles[lvl][din];

      // Note: EBParticleMesh construction is cheap.
      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleScalarField>(boxParticles.listItems(), meshData, a_depositionType, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary should end up to the coarse level. Add that mass right now.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
    }

    // 4. The particles on the other side of the refinement boundary should deposit into this level.
    if (hasCoar) {
      // Particles that live "around" the refinement boundary. Must have
      // been filled by the user.
      const ParticleData<P>& coarHaloParticles = *a_particles.getMaskParticles()[lvl - 1];

      // Get the refined coarse grid stuff.
      const int                refRat    = m_refRat[lvl - 1];
      const EBLevelGrid&       eblgFiCo  = m_coarseFinePM[lvl]->getEblgFiCo();
      const DisjointBoxLayout& dblFiCo   = eblgFiCo.getDBL();
      const EBISLayout&        ebislFiCo = eblgFiCo.getEBISL();
      const DataIterator&      ditFiCo   = dblFiCo.dataIterator();

      // Make a buffer we can deposit into.
      LevelData<EBCellFAB> bufferFiCo(dblFiCo, 1, m_ghost * IntVect::Unit, EBCellFactory(ebislFiCo));

      const int numBoxesFiCo = ditFiCo.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxesFiCo; mybox++) {
        const DataIndex& din = ditFiCo[mybox];

        EBCellFAB&        dataFiCo      = bufferFiCo[din];
        const ListBox<P>& haloParticles = coarHaloParticles[din];

        dataFiCo.setVal(0.0);

        const EBParticleMesh& interp = (*m_ebParticleMeshFiCo[lvl])[din];

        interp.deposit<P, const Real&, particleScalarField>(dataFiCo,
                                                            haloParticles.listItems(),
                                                            a_depositionType,
                                                            refRat,
                                                            a_forceIrregNGP);
      }

      // Add the result of the buffer deposition to this level.
      m_coarseFinePM[lvl]->addFiCoDataToFine(*a_meshData[lvl], bufferFiCo);
    }
  }

  a_particles.clearMaskParticles();
}

template <class P, Real (P::*particleScalarField)() const>
void
EBAMRParticleMesh::depositHalo(EBAMRCellData&              a_meshData,
                               const ParticleContainer<P>& a_particles,
                               const DepositionType        a_depositionType,
                               const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositHalo");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositHalo" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == 1);

  // TLDR: On each level the particles are deposited as usual. This routine permits particles to live right next to the refinement boundary
  //       so both the coarse-level and fine-level particles have clouds that hang into the fine/coarse level. The way we handle this is as follows:
  //       The fine-level particles will have deposited into ghost cells that lie around the refinement boundary. This mass is put on the coarse grid.
  //       On the coarse level the particles have clouds that extend beneath the fine grid. Instead of interpolating this mass to the fine grid we
  //       run a special deposition procedure which allows us to deposit directly on the fine grid using 2x or 4x the particle width (depending on the
  //       refinement ratio).

  DataOps::setValue(a_meshData, 0.0);

  const Interval interv(0, 0);

  // Copy the required particles to a masked particle data holder -- this takes the first strip of cells on the coarse side of the refinemnt boundary.
  constexpr int coarseMaskWidth = 1;

  a_particles.copyMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl   = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit     = dbl.dataIterator();
    const bool               hasCoar = (lvl > 0);

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&        meshData     = (*a_meshData[lvl])[din];
      const ListBox<P>& boxParticles = a_particles[lvl][din];

      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleScalarField>(boxParticles.listItems(), meshData, a_depositionType, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary should end up to the coarse level. Add that mass right now.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
    }

    // 4. The particles on the other side of the refinement boundary should deposit into this level.
    if (hasCoar) {
      // Particles that live "around" the refinement boundary. Must have
      // been filled by the user.
      const ParticleData<P>& coarHaloParticles = *a_particles.getMaskParticles()[lvl - 1];

      // Get the refined coarse grid stuff.
      const int                refRat    = m_refRat[lvl - 1];
      const EBLevelGrid&       eblgFiCo  = m_coarseFinePM[lvl]->getEblgFiCo();
      const DisjointBoxLayout& dblFiCo   = eblgFiCo.getDBL();
      const EBISLayout&        ebislFiCo = eblgFiCo.getEBISL();
      const DataIterator&      ditFiCo   = dblFiCo.dataIterator();

      // Make a buffer we can deposit into.
      LevelData<EBCellFAB> bufferFiCo(dblFiCo, 1, m_ghost * IntVect::Unit, EBCellFactory(ebislFiCo));

      const int numBoxesFiCo = ditFiCo.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxesFiCo; mybox++) {
        const DataIndex& din = ditFiCo[mybox];

        EBCellFAB&        dataFiCo      = bufferFiCo[din];
        const ListBox<P>& haloParticles = coarHaloParticles[din];

        dataFiCo.setVal(0.0);

        const EBParticleMesh& interp = (*m_ebParticleMeshFiCo[lvl])[din];

        interp.deposit<P, Real, particleScalarField>(dataFiCo,
                                                     haloParticles.listItems(),
                                                     a_depositionType,
                                                     1.0 * refRat,
                                                     a_forceIrregNGP);
      }

      // Add the result of the buffer deposition to this level.
      m_coarseFinePM[lvl]->addFiCoDataToFine(*a_meshData[lvl], bufferFiCo);
    }
  }

  a_particles.clearMaskParticles();
}

template <class P, const Real& (P::*particleScalarField)() const>
void
EBAMRParticleMesh::depositHaloNGP(EBAMRCellData&              a_meshData,
                                  const ParticleContainer<P>& a_particles,
                                  const DepositionType        a_depositionType,
                                  const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositHaloNGP");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositHaloNGP" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == 1);

  DataOps::setValue(a_meshData, 0.0);

  const Interval interv(0, 0);

  // Copy the required particles to a masked particle data holder -- this takes the first strip of cells on the coarse side of the refinemnt boundary.
  constexpr int coarseMaskWidth = 1;

  // Doing the nasty here...
  ParticleContainer<P>& particles = const_cast<ParticleContainer<P>&>(a_particles);

  particles.transferMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

  // nonHaloParticles will get deposited with the 'a_depositionType' scheme and haloParticles with an NGP scheme
  const AMRParticles<P>& nonHaloParticles = particles.getParticles();
  const AMRParticles<P>& haloParticles    = particles.getMaskParticles();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit   = dbl.dataIterator();

    const bool hasCoar = (lvl > 0);

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&        meshData       = (*a_meshData[lvl])[din];
      const ListBox<P>& particlesCloud = (*nonHaloParticles[lvl])[din];
      const ListBox<P>& particlesNGP   = (*haloParticles[lvl])[din];

      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleScalarField>(particlesCloud.listItems(), meshData, a_depositionType, a_forceIrregNGP);
      interp.deposit<P, particleScalarField>(particlesNGP.listItems(), meshData, DepositionType::NGP, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary should end up to the coarse level. Add that mass right now.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
    }
  }

  particles.transferParticles(particles.getMaskParticles());
}

template <class P, Real (P::*particleScalarField)() const>
void
EBAMRParticleMesh::depositHaloNGP(EBAMRCellData&              a_meshData,
                                  const ParticleContainer<P>& a_particles,
                                  const DepositionType        a_depositionType,
                                  const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositHaloNGP");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositHaloNGP" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == 1);

  DataOps::setValue(a_meshData, 0.0);

  const Interval interv(0, 0);

  // Doing the nasty here...
  ParticleContainer<P>& particles = const_cast<ParticleContainer<P>&>(a_particles);

  // Copy the required particles to a masked particle data holder -- this takes the first strip of cells on the coarse side of the refinemnt boundary.
  constexpr int coarseMaskWidth = 1;

  particles.transferMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

  // nonHaloParticles will get deposited with the 'a_depositionType' scheme and haloParticles with an NGP scheme
  const AMRParticles<P>& nonHaloParticles = particles.getParticles();
  const AMRParticles<P>& haloParticles    = particles.getMaskParticles();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit   = dbl.dataIterator();

    const bool hasCoar = (lvl > 0);

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&        meshData       = (*a_meshData[lvl])[din];
      const ListBox<P>& particlesCloud = (*nonHaloParticles[lvl])[din];
      const ListBox<P>& particlesNGP   = (*haloParticles[lvl])[din];

      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleScalarField>(particlesCloud.listItems(), meshData, a_depositionType, a_forceIrregNGP);
      interp.deposit<P, particleScalarField>(particlesNGP.listItems(), meshData, DepositionType::NGP, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary should end up to the coarse level. Add that mass right now.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
    }
  }

  particles.transferParticles(particles.getMaskParticles());
}

template <class P, const Real& (P::*particleScalarField)() const>
void
EBAMRParticleMesh::depositTransition(EBAMRCellData&              a_meshData,
                                     const ParticleContainer<P>& a_particles,
                                     const DepositionType        a_depositionType,
                                     const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositTransition");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositTransition" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == 1);

  // If we're calling with NGP, call a simpler function.
  if (a_depositionType == DepositionType::NGP) {
    depositInterp<P, particleScalarField>(a_meshData, a_particles, DepositionType::NGP, a_forceIrregNGP);

    return;
  }

  const int      numComp = 1;
  const int      comp    = 0;
  const Interval interv  = Interval(comp, comp);

  DataOps::setValue(a_meshData, 0.0);

  // This piece of code transfer the particles that lie on the coarse-side interface to a different particle container than
  // a_particles. We can not use ParticleContainer::transferMaskParticles because the mask is defined on the refined coarse
  // level. So we do this transfer directly.
  ParticleContainer<P>& particles = const_cast<ParticleContainer<P>&>(a_particles);

  this->transferMaskParticlesTransition(particles, a_depositionType);

  // Main deposition loop.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl      = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl    = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit      = dbl.dataIterator();
    const int                numBoxes = dit.size();

    const bool hasCoar = (lvl > 0);
    const bool hasFine = (lvl < m_finestLevel);

    // 1. Deposit particles on this level.
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&            meshData     = (*a_meshData[lvl])[din];
      const List<P>&        amrParticles = (*particles.getParticles()[lvl])[din].listItems();
      const EBParticleMesh& particleMesh = (*m_ebParticleMesh[lvl])[din];

      particleMesh.deposit<P, particleScalarField>(amrParticles, meshData, a_depositionType, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level so that all the mass is on this level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary must end up on the coarse level.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
    }

    // 4. If there is a finer level, not all of the particles on this level have yet been deposited because some of them were transferred into another
    //    particle container. We now run the following steps:
    //
    //    a) Deposit the particles around the refinement boundary on the refined coarse level, using the fine-grid particle width.
    //    b) Add mass from the refined coarse level to the fine level -- this puts some mass into the fine side of the refinement boundary
    //    c) Exchange data on the refined coarse level to update ghost cells, so that all special mass is contained in the refined coarse grid
    //    d) Coarsen mass from the refined level onto this level.
    if (hasFine) {
      // Note: eblgFiCo is for transferring data between lvl and lvl+1. It stores the refinement of grid level 'lvl' in eblgFiCo
      const EBLevelGrid&       eblgFiCo     = m_coarseFinePM[lvl + 1]->getEblgFiCo();
      const DisjointBoxLayout& dblFiCo      = eblgFiCo.getDBL();
      const EBISLayout&        ebislFiCo    = eblgFiCo.getEBISL();
      const DataIterator&      ditFiCo      = dblFiCo.dataIterator();
      const int                numBoxesFiCo = ditFiCo.size();

      // Make a buffer we can deposit into. This is a refined version of this level.
      LevelData<EBCellFAB> bufferFiCo(dblFiCo, 1, m_ghost * IntVect::Unit, EBCellFactory(ebislFiCo));

      const int maskWidth = this->getTransitionMaskWidth(a_depositionType, m_refRat[lvl]);

      // a) Deposit the particles on the refined coarse level.
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxesFiCo; mybox++) {
        const DataIndex& din = ditFiCo[mybox];

        EBCellFAB& dataFiCo = bufferFiCo[din];

        dataFiCo.setVal(0.0);

        const BaseFab<bool>&  mask          = (*m_transitionMasks.at(maskWidth)[lvl])[din];
        const List<P>&        maskParticles = (*particles.getMaskParticles()[lvl])[din].listItems();
        const EBParticleMesh& particleMesh  = (*m_ebParticleMeshFiCo[lvl + 1])[din];

        if (mask.isUsable()) {
          //          particleMesh.deposit<P, particleScalarField>(maskParticles, dataFiCo, DepositionType::NGP, a_forceIrregNGP);
          particleMesh.deposit<P, particleScalarField>(maskParticles, dataFiCo, a_depositionType, a_forceIrregNGP);
        }
      }

      // b) Add the data to the fine level. This moves from valid+ghost -> valid
      m_coarseFinePM[lvl + 1]->addFiCoDataToFine(*a_meshData[lvl + 1], bufferFiCo);

      // c) Exchange data on this level
      m_coarseFinePM[lvl + 1]->exchangeAndAddFiCoData(bufferFiCo);

      // d) Coarsen data from the refined grid to this grid.
      m_coarseFinePM[lvl + 1]->restrictAndAddFiCoDataToCoar(*a_meshData[lvl],
                                                            bufferFiCo,
                                                            EBCoarseFineParticleMesh::Average::Arithmetic);
    }
  }

  // Masked particles are but back in their correct mesh.
  particles.transferParticles(particles.getMaskParticles());
}

template <class P, Real (P::*particleScalarField)() const>
void
EBAMRParticleMesh::depositTransition(EBAMRCellData&              a_meshData,
                                     const ParticleContainer<P>& a_particles,
                                     const DepositionType        a_depositionType,
                                     const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositTransition");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositTransition" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == 1);

  // If we're calling with NGP, call a simpler function.
  if (a_depositionType == DepositionType::NGP) {
    depositInterp<P, particleScalarField>(a_meshData, a_particles, DepositionType::NGP, a_forceIrregNGP);

    return;
  }

  const int      numComp = 1;
  const int      comp    = 0;
  const Interval interv  = Interval(comp, comp);

  DataOps::setValue(a_meshData, 0.0);

  // This piece of code transfer the particles that lie on the coarse-side interface to a different particle container than
  // a_particles. We can not use ParticleContainer::transferMaskParticles because the mask is defined on the refined coarse
  // level. So we do this transfer directly.
  ParticleContainer<P>& particles = const_cast<ParticleContainer<P>&>(a_particles);

  this->transferMaskParticlesTransition(particles, a_depositionType);

  // Main deposition loop.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl      = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl    = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit      = dbl.dataIterator();
    const int                numBoxes = dit.size();

    const bool hasCoar = (lvl > 0);
    const bool hasFine = (lvl < m_finestLevel);

    // 1. Deposit particles on this level.
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&            meshData     = (*a_meshData[lvl])[din];
      const List<P>&        amrParticles = (*particles.getParticles()[lvl])[din].listItems();
      const EBParticleMesh& particleMesh = (*m_ebParticleMesh[lvl])[din];

      particleMesh.deposit<P, particleScalarField>(amrParticles, meshData, a_depositionType, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level so that all the mass is on this level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary must end up on the coarse level.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
    }

    // 4. If there is a finer level, not all of the particles on this level have yet been deposited because some of them were transferred into another
    //    particle container. We now run the following steps:
    //
    //    a) Deposit the particles around the refinement boundary on the refined coarse level, using the fine-grid particle width.
    //    b) Add mass from the refined coarse level to the fine level -- this puts some mass into the fine side of the refinement boundary
    //    c) Exchange data on the refined coarse level to update ghost cells, so that all special mass is contained in the refined coarse grid
    //    d) Coarsen mass from the refined level onto this level.
    if (hasFine) {
      // Note: eblgFiCo is for transferring data between lvl and lvl+1. It stores the refinement of grid level 'lvl' in eblgFiCo
      const EBLevelGrid&       eblgFiCo     = m_coarseFinePM[lvl + 1]->getEblgFiCo();
      const DisjointBoxLayout& dblFiCo      = eblgFiCo.getDBL();
      const EBISLayout&        ebislFiCo    = eblgFiCo.getEBISL();
      const DataIterator&      ditFiCo      = dblFiCo.dataIterator();
      const int                numBoxesFiCo = ditFiCo.size();

      // Make a buffer we can deposit into. This is a refined version of this level.
      LevelData<EBCellFAB> bufferFiCo(dblFiCo, 1, m_ghost * IntVect::Unit, EBCellFactory(ebislFiCo));

      const int maskWidth = this->getTransitionMaskWidth(a_depositionType, m_refRat[lvl]);

      // a) Deposit the particles on the refined coarse level.
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxesFiCo; mybox++) {
        const DataIndex& din = ditFiCo[mybox];

        EBCellFAB& dataFiCo = bufferFiCo[din];

        dataFiCo.setVal(0.0);

        const BaseFab<bool>&  mask          = (*m_transitionMasks.at(maskWidth)[lvl])[din];
        const List<P>&        maskParticles = (*particles.getMaskParticles()[lvl])[din].listItems();
        const EBParticleMesh& particleMesh  = (*m_ebParticleMeshFiCo[lvl + 1])[din];

        if (mask.isUsable()) {
          //          particleMesh.deposit<P, particleScalarField>(maskParticles, dataFiCo, DepositionType::NGP, a_forceIrregNGP);
          particleMesh.deposit<P, particleScalarField>(maskParticles, dataFiCo, a_depositionType, a_forceIrregNGP);
        }
      }

      // b) Add the data to the fine level. This moves from valid+ghost -> valid
      m_coarseFinePM[lvl + 1]->addFiCoDataToFine(*a_meshData[lvl + 1], bufferFiCo);

      // c) Exchange data on this level
      m_coarseFinePM[lvl + 1]->exchangeAndAddFiCoData(bufferFiCo);

      // d) Coarsen data from the refined grid to this grid.
      m_coarseFinePM[lvl + 1]->restrictAndAddFiCoDataToCoar(*a_meshData[lvl],
                                                            bufferFiCo,
                                                            EBCoarseFineParticleMesh::Average::Arithmetic);
    }
  }

  // Masked particles are but back in their correct mesh.
  particles.transferParticles(particles.getMaskParticles());
}

template <class P, Real& (P::*particleScalarField)()>
void
EBAMRParticleMesh::interpolate(ParticleContainer<P>& a_particles,
                               const EBAMRCellData&  a_meshScalarField,
                               const DepositionType  a_interpType,
                               const bool            a_forceIrregNGP) const
{
  CH_TIME("EBAMRParticleMesh::interpolate");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::interpolate" << endl;
  }

  CH_assert(m_isDefined);

  // TLDR: Run through each patch and interpolate to the particle positions.

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const EBLevelGrid&       eblg  = *m_eblgs[lvl];
    const DisjointBoxLayout& dbl   = eblg.getDBL();
    const EBISLayout&        ebisl = eblg.getEBISL();
    const DataIterator&      dit   = dbl.dataIterator();

    CH_assert(a_meshScalarField[lvl]->nComp() == 1);

    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& particles = a_particles[lvl][din].listItems();

      const EBCellFAB&      data   = (*a_meshScalarField[lvl])[din];
      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.interpolate<P, particleScalarField>(particles, data, a_interpType, a_forceIrregNGP);
    }
  }
}

template <class P, RealVect& (P::*particleVectorField)()>
void
EBAMRParticleMesh::interpolate(ParticleContainer<P>& a_particles,
                               const EBAMRCellData&  a_meshVectorField,
                               const DepositionType  a_interpType,
                               const bool            a_forceIrregNGP) const
{
  CH_TIME("EBAMRParticleMesh::interpolate");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::interpolate" << endl;
  }

  CH_assert(m_isDefined);

  // TLDR: Run through each patch and interpolate to the particle positions.

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const EBLevelGrid&       eblg   = *m_eblgs[lvl];
    const DisjointBoxLayout& dbl    = eblg.getDBL();
    const ProblemDomain&     domain = eblg.getDomain();
    const EBISLayout&        ebisl  = eblg.getEBISL();
    const Real&              dx     = m_dx[lvl];
    const DataIterator&      dit    = dbl.dataIterator();

    CH_assert(a_meshVectorField[lvl]->nComp() == SpaceDim);

    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& particles = a_particles[lvl][din].listItems();

      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];
      const EBCellFAB&      data   = (*a_meshVectorField[lvl])[din];

      interp.interpolate<P, particleVectorField>(particles, data, a_interpType, a_forceIrregNGP);
    }
  }
}

template <class P, const RealVect& (P::*particleVectorField)() const>
void
EBAMRParticleMesh::depositInterp(EBAMRCellData&              a_meshData,
                                 const ParticleContainer<P>& a_particles,
                                 const DepositionType        a_depositionType,
                                 const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositInterp");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositInterp" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == SpaceDim);

  DataOps::setValue(a_meshData, 0.0);

  const Interval interv(0, 0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
    const DataIterator&      dit     = dbl.dataIterator();
    const bool               hasCoar = lvl > 0;

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&        meshData     = (*a_meshData[lvl])[din];
      const ListBox<P>& boxParticles = a_particles[lvl][din];

      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleVectorField>(boxParticles.listItems(), meshData, a_depositionType, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level. After this, all the mass should be on the current level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. If the particles deposited over the coarse-fine boundary, add the mass to the coarse level. Since
    //    the coarse-fine choreography buffer only uses a single component in the memory, we need to alias the data.
    if (hasCoar) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        LevelData<EBCellFAB> coarAlias;
        LevelData<EBCellFAB> fineAlias;

        // This lets coarAlias be a 1-component data holder which references the dir-component
        // in a_meshData.
        aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(dir, dir));
        aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(dir, dir));

        // Add the mass that hangs from the fine level and over the refinement boundary onto the coarse level.
        m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);

        // Likewise, take the particles that deposited mass to underneath the current level
        // and put that mas on the fine level.
        m_coarseFinePM[lvl]->addInvalidCoarseToFine(fineAlias, coarAlias);
      }
    }
  }
}

template <class P, RealVect (P::*particleVectorField)() const>
void
EBAMRParticleMesh::depositInterp(EBAMRCellData&              a_meshData,
                                 const ParticleContainer<P>& a_particles,
                                 const DepositionType        a_depositionType,
                                 const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositInterp");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositInterp" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == SpaceDim);

  DataOps::setValue(a_meshData, 0.0);

  const Interval interv(0, 0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
    const DataIterator&      dit     = dbl.dataIterator();
    const bool               hasCoar = lvl > 0;

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&            meshData     = (*a_meshData[lvl])[din];
      const ListBox<P>&     boxParticles = a_particles[lvl][din];
      const EBParticleMesh& interp       = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleVectorField>(boxParticles.listItems(), meshData, a_depositionType, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level. After this, all the mass should be on the current level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. If the particles deposited over the coarse-fine boundary, add the mass to the coarse level. Since
    //    the coarse-fine choreography buffer only uses a single component in the memory, we need to alias the data.
    if (hasCoar) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        LevelData<EBCellFAB> coarAlias;
        LevelData<EBCellFAB> fineAlias;

        // This lets coarAlias be a 1-component data holder which references the dir-component
        // in a_meshData.
        aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(dir, dir));
        aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(dir, dir));

        // Add the mass that hangs from the fine level and over the refinement boundary onto the coarse level.
        m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);

        // Likewise, take the particles that deposited mass to underneath the current level
        // and put that mas on the fine level.
        m_coarseFinePM[lvl]->addInvalidCoarseToFine(fineAlias, coarAlias);
      }
    }
  }
}

template <class P, const RealVect& (P::*particleVectorField)() const>
void
EBAMRParticleMesh::depositHalo(EBAMRCellData&              a_meshData,
                               const ParticleContainer<P>& a_particles,
                               const DepositionType        a_depositionType,
                               const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositHalo");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositHalo" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == SpaceDim);

  DataOps::setValue(a_meshData, 0.0);

  // Copy the required particles to a masked particle data holder -- this takes the first strip of cells on the coarse side of the refinemnt boundary.
  constexpr int coarseMaskWidth = 1;

  a_particles.copyMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

  const Interval interv(0, 0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit   = dbl.dataIterator();

    const bool hasCoar = (lvl > 0);

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&        meshData     = (*a_meshData[lvl])[din];
      const ListBox<P>& boxParticles = a_particles[lvl][din];

      // Note: EBParticleMesh construction is cheap.
      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleVectorField>(boxParticles.listItems(), meshData, a_depositionType, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary should end up to the coarse level. Add that mass right now.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
    }

    // 4. The particles on the other side of the refinement boundary should deposit into this level.
    if (hasCoar) {
      // Particles that live "around" the refinement boundary. Must have
      // been filled by the user.
      const ParticleData<P>& coarHaloParticles = *a_particles.getMaskParticles()[lvl - 1];

      // Get the refined coarse grid stuff.
      const int                refRat    = m_refRat[lvl - 1];
      const EBLevelGrid&       eblgFiCo  = m_coarseFinePM[lvl]->getEblgFiCo();
      const DisjointBoxLayout& dblFiCo   = eblgFiCo.getDBL();
      const EBISLayout&        ebislFiCo = eblgFiCo.getEBISL();
      const DataIterator&      ditFiCo   = dblFiCo.dataIterator();

      // Make a buffer we can deposit into.
      LevelData<EBCellFAB> bufferFiCo(dblFiCo, SpaceDim, m_ghost * IntVect::Unit, EBCellFactory(ebislFiCo));

      const int numBoxesFiCo = ditFiCo.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxesFiCo; mybox++) {
        const DataIndex&  din           = ditFiCo[mybox];
        EBCellFAB&        dataFiCo      = bufferFiCo[din];
        const ListBox<P>& haloParticles = coarHaloParticles[din];

        dataFiCo.setVal(0.0);

        const EBParticleMesh& interp = (*m_ebParticleMeshFiCo[lvl])[din];

        switch (refRat) {
        case 2: {
          interp.deposit2<P, particleVectorField>(haloParticles.listItems(),
                                                  dataFiCo,
                                                  a_depositionType,
                                                  a_forceIrregNGP);

          break;
        }
        case 4: {
          interp.deposit4<P, particleVectorField>(haloParticles.listItems(),
                                                  dataFiCo,
                                                  a_depositionType,
                                                  a_forceIrregNGP);

          break;
        }
        default: {
          MayDay::Error("CD_EBAMRParticleMeshImplem.H - logic bust in EBAMRParticleMesh::depositHalo(RealVect)");

          break;
        }
        }
      }

      // Add the result of the buffer deposition to this level. The buffer above had SpaceDim// components but coarseFinePM works with only one component, so we just alias the data.
      for (int dir = 0; dir < SpaceDim; dir++) {
        LevelData<EBCellFAB> meshDataAlias;
        LevelData<EBCellFAB> bufferFiCoAlias;

        aliasLevelData<EBCellFAB>(meshDataAlias, &(*a_meshData[lvl]), Interval(dir, dir));
        aliasLevelData<EBCellFAB>(bufferFiCoAlias, &bufferFiCo, Interval(dir, dir));

        m_coarseFinePM[lvl]->addFiCoDataToFine(meshDataAlias, bufferFiCoAlias);
      }
    }
  }

  a_particles.clearMaskParticles();
}

template <class P, RealVect (P::*particleVectorField)() const>
void
EBAMRParticleMesh::depositHalo(EBAMRCellData&              a_meshData,
                               const ParticleContainer<P>& a_particles,
                               const DepositionType        a_depositionType,
                               const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositHalo");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositHalo" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == SpaceDim);

  DataOps::setValue(a_meshData, 0.0);

  // Copy the required particles to a masked particle data holder -- this takes the first strip of cells on the coarse side of the refinemnt boundary.
  constexpr int coarseMaskWidth = 1;

  a_particles.copyMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

  const Interval interv(0, 0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit   = dbl.dataIterator();

    const bool hasCoar = (lvl > 0);

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&        meshData     = (*a_meshData[lvl])[din];
      const ListBox<P>& boxParticles = a_particles[lvl][din];

      // Note: EBParticleMesh construction is cheap.
      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleVectorField>(boxParticles.listItems(), meshData, a_depositionType, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary should end up to the coarse level. Add that mass right now.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
    }

    // 4. The particles on the other side of the refinement boundary should deposit into this level.
    if (hasCoar) {
      // Particles that live "around" the refinement boundary. Must have
      // been filled by the user.
      const ParticleData<P>& coarHaloParticles = *a_particles.getMaskParticles()[lvl - 1];

      // Get the refined coarse grid stuff.
      const int                refRat    = m_refRat[lvl - 1];
      const EBLevelGrid&       eblgFiCo  = m_coarseFinePM[lvl]->getEblgFiCo();
      const DisjointBoxLayout& dblFiCo   = eblgFiCo.getDBL();
      const EBISLayout&        ebislFiCo = eblgFiCo.getEBISL();
      const DataIterator&      ditFiCo   = dblFiCo.dataIterator();

      // Make a buffer we can deposit into.
      LevelData<EBCellFAB> bufferFiCo(dblFiCo, SpaceDim, m_ghost * IntVect::Unit, EBCellFactory(ebislFiCo));

      const int numBoxesFiCo = ditFiCo.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxesFiCo; mybox++) {
        const DataIndex& din = ditFiCo[mybox];

        EBCellFAB&        dataFiCo      = bufferFiCo[din];
        const ListBox<P>& haloParticles = coarHaloParticles[din];

        dataFiCo.setVal(0.0);

        const EBParticleMesh& interp = (*m_ebParticleMeshFiCo[lvl])[din];

        switch (refRat) {
        case 2: {
          interp.deposit2<P, particleVectorField>(haloParticles.listItems(),
                                                  dataFiCo,
                                                  a_depositionType,
                                                  a_forceIrregNGP);

          break;
        }
        case 4: {
          interp.deposit4<P, particleVectorField>(haloParticles.listItems(),
                                                  dataFiCo,
                                                  a_depositionType,
                                                  a_forceIrregNGP);

          break;
        }
        default: {
          MayDay::Error("CD_EBAMRParticleMeshImplem.H - logic bust in EBAMRParticleMesh::depositHalo(RealVect)");

          break;
        }
        }
      }

      // Add the result of the buffer deposition to this level. The buffer above had SpaceDim
      // components but coarseFinePM works with only one component, so we just alias the data.
      for (int dir = 0; dir < SpaceDim; dir++) {
        LevelData<EBCellFAB> meshDataAlias;
        LevelData<EBCellFAB> bufferFiCoAlias;

        aliasLevelData<EBCellFAB>(meshDataAlias, &(*a_meshData[lvl]), Interval(dir, dir));
        aliasLevelData<EBCellFAB>(bufferFiCoAlias, &bufferFiCo, Interval(dir, dir));

        m_coarseFinePM[lvl]->addFiCoDataToFine(meshDataAlias, bufferFiCoAlias);
      }
    }
  }

  a_particles.clearMaskParticles();
}

template <class P, const RealVect& (P::*particleVectorField)() const>
void
EBAMRParticleMesh::depositHaloNGP(EBAMRCellData&              a_meshData,
                                  const ParticleContainer<P>& a_particles,
                                  const DepositionType        a_depositionType,
                                  const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositHaloNGP");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositHaloNGP" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == SpaceDim);

  // TLDR: On each level the particles are deposited as usual. This routine permits particles to live right next to the refinement boundary
  //       so both the coarse-level and fine-level particles have clouds that hang into the fine/coarse level. The way we handle this is as follows:
  //       The fine-level particles will have deposited into ghost cells that lie around the refinement boundary. This mass is put on the coarse grid.
  //       On the coarse level the particles have clouds that extend beneath the fine grid. Instead of interpolating this mass to the fine grid we
  //       run a special deposition procedure which allows us to deposit directly on the fine grid using 2x or 4x the particle width (depending on the
  //       refinement ratio).

  DataOps::setValue(a_meshData, 0.0);

  const Interval interv(0, 0);

  // Copy the required particles to a masked particle data holder -- this takes the first strip of cells on the coarse side of the refinemnt boundary.
  constexpr int coarseMaskWidth = 1;

  // Doing the nasty here...
  ParticleContainer<P>& particles = const_cast<ParticleContainer<P>&>(a_particles);

  particles.transferMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

  // nonHaloParticles will get deposited with the 'a_depositionType' scheme and haloParticles with an NGP scheme
  const AMRParticles<P>& nonHaloParticles = particles.getParticles();
  const AMRParticles<P>& haloParticles    = particles.getMaskParticles();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit   = dbl.dataIterator();

    const bool hasCoar = (lvl > 0);

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&        meshData       = (*a_meshData[lvl])[din];
      const ListBox<P>& particlesCloud = (*nonHaloParticles[lvl])[din];
      const ListBox<P>& particlesNGP   = (*haloParticles[lvl])[din];

      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleVectorField>(particlesCloud.listItems(), meshData, a_depositionType, a_forceIrregNGP);
      interp.deposit<P, particleVectorField>(particlesNGP.listItems(), meshData, DepositionType::NGP, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary should end up to the coarse level. Add that mass right now.
    if (hasCoar) {

      // Because m_coarseFinePM works with one component but we have SpaceDim.
      for (int dir = 0; dir < SpaceDim; dir++) {
        LevelData<EBCellFAB> coarAlias;
        LevelData<EBCellFAB> fineAlias;

        // This lets coarAlias be a 1-component data holder which references the dir-component
        // in a_meshData.
        aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(dir, dir));
        aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(dir, dir));

        // Average data in the fine-level ghost cells and add it to the coarse level.
        m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);
      }
    }
  }

  particles.transferParticles(particles.getMaskParticles());
}

template <class P, RealVect (P::*particleVectorField)() const>
void
EBAMRParticleMesh::depositHaloNGP(EBAMRCellData&              a_meshData,
                                  const ParticleContainer<P>& a_particles,
                                  const DepositionType        a_depositionType,
                                  const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositHaloNGP");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositHaloNGP" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == SpaceDim);

  // TLDR: On each level the particles are deposited as usual. This routine permits particles to live right next to the refinement boundary
  //       so both the coarse-level and fine-level particles have clouds that hang into the fine/coarse level. The way we handle this is as follows:
  //       The fine-level particles will have deposited into ghost cells that lie around the refinement boundary. This mass is put on the coarse grid.
  //       On the coarse level the particles have clouds that extend beneath the fine grid. Instead of interpolating this mass to the fine grid we
  //       run a special deposition procedure which allows us to deposit directly on the fine grid using 2x or 4x the particle width (depending on the
  //       refinement ratio).

  DataOps::setValue(a_meshData, 0.0);

  const Interval interv(0, 0);

  // Copy the required particles to a masked particle data holder -- this takes the first strip of cells on the coarse side of the refinemnt boundary.
  constexpr int coarseMaskWidth = 1;

  // Doing the nasty here...
  ParticleContainer<P>& particles = const_cast<ParticleContainer<P>&>(a_particles);

  particles.transferMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

  // nonHaloParticles will get deposited with the 'a_depositionType' scheme and haloParticles with an NGP scheme
  const AMRParticles<P>& nonHaloParticles = particles.getParticles();
  const AMRParticles<P>& haloParticles    = particles.getMaskParticles();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit   = dbl.dataIterator();

    const bool hasCoar = (lvl > 0);

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&        meshData       = (*a_meshData[lvl])[din];
      const ListBox<P>& particlesCloud = (*nonHaloParticles[lvl])[din];
      const ListBox<P>& particlesNGP   = (*haloParticles[lvl])[din];

      const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

      interp.deposit<P, particleVectorField>(particlesCloud.listItems(), meshData, a_depositionType, a_forceIrregNGP);
      interp.deposit<P, particleVectorField>(particlesNGP.listItems(), meshData, DepositionType::NGP, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level.
    a_meshData[lvl]->exchange(interv, m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary should end up to the coarse level. Add that mass right now.
    if (hasCoar) {

      // Because m_coarseFinePM works with one component but we have SpaceDim.
      for (int dir = 0; dir < SpaceDim; dir++) {
        LevelData<EBCellFAB> coarAlias;
        LevelData<EBCellFAB> fineAlias;

        // This lets coarAlias be a 1-component data holder which references the dir-component
        // in a_meshData.
        aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(dir, dir));
        aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(dir, dir));

        // Average data in the fine-level ghost cells and add it to the coarse level.
        m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);
      }
    }
  }

  particles.transferParticles(particles.getMaskParticles());
}

template <class P, const RealVect& (P::*particleVectorField)() const>
void
EBAMRParticleMesh::depositTransition(EBAMRCellData&              a_meshData,
                                     const ParticleContainer<P>& a_particles,
                                     const DepositionType        a_depositionType,
                                     const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositTransition");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositTransition" << endl;
  }

  MayDay::Abort("EBAMRParticleMesh::depositTransition - not implemented");
}

template <class P, RealVect (P::*particleVectorField)() const>
void
EBAMRParticleMesh::depositTransition(EBAMRCellData&              a_meshData,
                                     const ParticleContainer<P>& a_particles,
                                     const DepositionType        a_depositionType,
                                     const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::depositTransition");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::depositTransition" << endl;
  }

  MayDay::Abort("EBAMRParticleMesh::depositTransition - not implemented");
}

template <class P>
void
EBAMRParticleMesh::transferMaskParticlesTransition(ParticleContainer<P>& a_particles,
                                                   const DepositionType  a_depositionType) const
{
  CH_TIME("EBAMRParticleMesh::transferMaskParticlesTransition");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::transferMaskParticlesTransition" << endl;
  }

  for (int lvl = 0; lvl < m_finestLevel; lvl++) {
    const EBLevelGrid& eblg     = *m_eblgs[lvl];
    const EBLevelGrid& eblgFiCo = m_coarseFinePM[lvl + 1]->getEblgFiCo();

    const DisjointBoxLayout& dbl     = eblg.getDBL();
    const DisjointBoxLayout& dblFiCo = eblgFiCo.getDBL();

    const DataIterator& dit     = dbl.dataIterator();
    const DataIterator& ditFiCo = dblFiCo.dataIterator();

    const int numBoxes     = dit.size();
    const int numBoxesFiCo = ditFiCo.size();

    const Real dx     = m_dx[lvl];
    const Real dxFine = m_dx[lvl + 1];

    if (a_depositionType != DepositionType::NGP) {
      const int maskWidth = this->getTransitionMaskWidth(a_depositionType, m_refRat[lvl]);

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxes; mybox++) {
        const DataIndex& din = dit[mybox];

        const BaseFab<bool>& mask = (*m_transitionMasks.at(maskWidth)[lvl])[din];

        if (mask.isUsable()) {
          List<P>& amrParticles  = (*a_particles.getParticles()[lvl])[din].listItems();
          List<P>& maskParticles = (*a_particles.getMaskParticles()[lvl])[din].listItems();

          const Box cellBox = dblFiCo[din];
          const Box maskBox = mask.box();

          CH_assert(cellBox == maskBox);

          for (ListIterator<P> lit(amrParticles); lit.ok();) {
            const IntVect iv = ParticleOps::getParticleCellIndex(lit().position(), m_probLo, dxFine);

            if (mask(iv)) {
              maskParticles.transfer(lit);
            }
            else {
              ++lit;
            }
          }
        }
      }
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
