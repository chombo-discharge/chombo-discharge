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

template <class P, class Ret, Ret (P::*MemberFunc)() const>
void
EBAMRParticleMesh::deposit(EBAMRCellData&              a_meshData,
                           const ParticleContainer<P>& a_particles,
                           const DepositionType        a_depositionType,
                           const CoarseFineDeposition  a_coarseFineDeposition,
                           const bool                  a_forceIrregNGP)
{
  CH_TIME("EBAMRParticleMesh::deposit");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::deposit" << endl;
  }

  switch (a_coarseFineDeposition) {
  case CoarseFineDeposition::Interp: {
    this->depositInterp<P, Ret, MemberFunc>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::Halo: {
    this->depositHalo<P, Ret, MemberFunc>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::HaloNGP: {
    this->depositHaloNGP<P, Ret, MemberFunc>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  case CoarseFineDeposition::Transition: {
    this->depositTransition<P, Ret, MemberFunc>(a_meshData, a_particles, a_depositionType, a_forceIrregNGP);

    break;
  }
  default: {
    MayDay::Error("EBAMRParticleMesh::deposit - logic bust");

    break;
  }
  }
}

template <class P, class Ret, Ret (P::*MemberFunc)() const>
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

  constexpr int numComp = EBParticleMesh::sanitize<Ret>();

  CH_assert(a_meshData[0]->nComp() == numComp);

  // Reset mesh data
  DataOps::setValue(a_meshData, 0.0);

  // Start deposition loop
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
    const DataIterator&      dit     = dbl.dataIterator();
    const bool               hasCoar = (lvl > 0);

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex&      din       = dit[mybox];
      const EBParticleMesh& interp    = (*m_ebParticleMesh[lvl])[din];
      const List<P>&        particles = a_particles[lvl][din].listItems();

      EBCellFAB& meshData = (*a_meshData[lvl])[din];

      interp.deposit<P, Ret, MemberFunc>(meshData, particles, a_depositionType, 1.0, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level. After this, all the mass should be on the current level.
    a_meshData[lvl]->exchange(Interval(0, numComp - 1), m_levelCopiers[lvl], EBAddOp());

    // 3. If the particles deposited over the coarse-fine boundary, add the mass to the coarse level. Since
    //    the coarse-fine choreography buffer only uses a single component in the memory, we need to alias the data.
    if (hasCoar) {
      for (int comp = 0; comp < numComp; comp++) {
        LevelData<EBCellFAB> coarAlias;
        LevelData<EBCellFAB> fineAlias;

        // This lets coarAlias be a 1-component data holder which references the dir-component
        // in a_meshData.
        aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(comp, comp));
        aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(comp, comp));

        // Add the mass that hangs from the fine level and over the refinement boundary onto the coarse level.
        m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);

        // Likewise, take the particles that deposited mass to underneath the current level
        // and put that mas on the fine level.
        m_coarseFinePM[lvl]->addInvalidCoarseToFine(fineAlias, coarAlias);
      }
    }
  }
}

template <class P, class Ret, Ret (P::*MemberFunc)() const>
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

  constexpr int numComp = EBParticleMesh::sanitize<Ret>();

  CH_assert(a_meshData[0]->nComp() == numComp);

  // Reset mesh data
  DataOps::setValue(a_meshData, 0.0);

  // Copy the required particles to the masked particle data holder.
  constexpr int coarseMaskWidth = 1;

  a_particles.copyMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl = m_eblgs[lvl]->getEBISL();
    const DataIterator&      dit   = dbl.dataIterator();

    const bool hasCoar = (lvl > 0);

    // 1. Deposit particles on this level.
    const int numBoxes = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex&      din       = dit[mybox];
      const EBParticleMesh& interp    = (*m_ebParticleMesh[lvl])[din];
      const List<P>&        particles = a_particles[lvl][din].listItems();

      EBCellFAB& meshData = (*a_meshData[lvl])[din];

      interp.deposit<P, Ret, MemberFunc>(meshData, particles, a_depositionType, 1.0, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level.
    a_meshData[lvl]->exchange(Interval(0, numComp - 1), m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary should end up to the coarse level. Add that mass right now.
    if (hasCoar) {
      m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
    }

    // 4. The particles on the other side of the refinement boundary should deposit into this level and also the coarser level.
    if (hasCoar) {
      const ParticleData<P>& coarHaloParticles = *a_particles.getMaskParticles()[lvl - 1];

      // Get the refined coarse grid stuff.
      const int                refRat    = m_refRat[lvl - 1];
      const EBLevelGrid&       eblgFiCo  = m_coarseFinePM[lvl]->getEblgFiCo();
      const DisjointBoxLayout& dblFiCo   = eblgFiCo.getDBL();
      const EBISLayout&        ebislFiCo = eblgFiCo.getEBISL();
      const DataIterator&      ditFiCo   = dblFiCo.dataIterator();

      // Make a buffer we can deposit into.
      LevelData<EBCellFAB> bufferFiCo(dblFiCo, numComp, m_ghost * IntVect::Unit, EBCellFactory(ebislFiCo));

      const int numBoxesFiCo = ditFiCo.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxesFiCo; mybox++) {
        const DataIndex&  din           = ditFiCo[mybox];
        EBCellFAB&        dataFiCo      = bufferFiCo[din];
        const ListBox<P>& haloParticles = coarHaloParticles[din];

        dataFiCo.setVal(0.0);

        const EBParticleMesh& interp = (*m_ebParticleMeshFiCo[lvl])[din];

        interp.deposit<P, Ret, MemberFunc>(dataFiCo,
                                           haloParticles.listItems(),
                                           a_depositionType,
                                           refRat,
                                           a_forceIrregNGP);
      }

      // Add the result of the buffer deposition to this level. The buffer above could have more than one component, while coarseFinePM works with only one component. The
      // solution is just to alias the damn thing.
      for (int comp = 0; comp < numComp; comp++) {
        LevelData<EBCellFAB> meshDataAlias;
        LevelData<EBCellFAB> bufferFiCoAlias;

        aliasLevelData<EBCellFAB>(meshDataAlias, &(*a_meshData[lvl]), Interval(comp, comp));
        aliasLevelData<EBCellFAB>(bufferFiCoAlias, &bufferFiCo, Interval(comp, comp));

        m_coarseFinePM[lvl]->addFiCoDataToFine(meshDataAlias, bufferFiCoAlias);
      }
    }
  }

  a_particles.clearMaskParticles();
}

template <class P, class Ret, Ret (P::*MemberFunc)() const>
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

  constexpr int numComp = EBParticleMesh::sanitize<Ret>();

  CH_assert(a_meshData[0]->nComp() == numComp);

  // Reset mesh data
  DataOps::setValue(a_meshData, 0.0);

  // Copy the required particles to a masked particle data holder.
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
      const DataIndex&      din            = dit[mybox];
      const EBParticleMesh& interp         = (*m_ebParticleMesh[lvl])[din];
      const List<P>&        particlesCloud = (*nonHaloParticles[lvl])[din].listItems();
      const List<P>&        particlesNGP   = (*haloParticles[lvl])[din].listItems();

      EBCellFAB& meshData = (*a_meshData[lvl])[din];

      interp.deposit<P, Ret, MemberFunc>(meshData, particlesCloud, a_depositionType, 1.0, a_forceIrregNGP);
      interp.deposit<P, Ret, MemberFunc>(meshData, particlesNGP, DepositionType::NGP, 1.0, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level.
    a_meshData[lvl]->exchange(Interval(0, numComp - 1), m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary should end up to the coarse level. Add that mass right now.
    if (hasCoar) {

      // Because m_coarseFinePM works with one component but we have SpaceDim.
      for (int comp = 0; comp < numComp; comp++) {
        LevelData<EBCellFAB> coarAlias;
        LevelData<EBCellFAB> fineAlias;

        // This lets coarAlias be a 1-component data holder which references the comp-component
        // in a_meshData.
        aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(comp, comp));
        aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(comp, comp));

        // Average data in the fine-level ghost cells and add it to the coarse level.
        m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);
      }
    }
  }

  particles.transferParticles(particles.getMaskParticles());
}

template <class P, class Ret, Ret (P::*MemberFunc)() const>
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

  constexpr int numComp = EBParticleMesh::sanitize<Ret>();

  CH_assert(a_meshData[0]->nComp() == numComp);

  // If we're calling with NGP, use a cheaper version.
  if (a_depositionType == DepositionType::NGP) {
    EBAMRParticleMesh::depositInterp<P, Ret, MemberFunc>(a_meshData, a_particles, DepositionType::NGP, a_forceIrregNGP);

    return;
  }

  // Reset mesh data
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
      const DataIndex&      din          = dit[mybox];
      const EBParticleMesh& particleMesh = (*m_ebParticleMesh[lvl])[din];
      const List<P>&        amrParticles = (*particles.getParticles()[lvl])[din].listItems();

      EBCellFAB& meshData = (*a_meshData[lvl])[din];

      particleMesh.deposit<P, Ret, MemberFunc>(meshData, amrParticles, a_depositionType, 1.0, a_forceIrregNGP);
    }

    // 2. Exchange ghost data on this level so that all the mass is on this level.
    a_meshData[lvl]->exchange(Interval(0, numComp - 1), m_levelCopiers[lvl], EBAddOp());

    // 3. Deposition into ghost cells across the refinement boundary must end up on the coarse level.
    if (hasCoar) {
      for (int comp = 0; comp < numComp; comp++) {
        LevelData<EBCellFAB> coarAlias;
        LevelData<EBCellFAB> fineAlias;

        aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(comp, comp));
        aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(comp, comp));

        m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);
      }
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
      LevelData<EBCellFAB> bufferFiCo(dblFiCo, numComp, m_ghost * IntVect::Unit, EBCellFactory(ebislFiCo));

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
          particleMesh.deposit<P, Ret, MemberFunc>(dataFiCo, maskParticles, a_depositionType, 1.0, a_forceIrregNGP);
        }
      }

      // b) Add the data to the fine level. This moves from valid+ghost -> valid
      m_coarseFinePM[lvl + 1]->addFiCoDataToFine(*a_meshData[lvl + 1], bufferFiCo);

      // c) Exchange data on this level
      m_coarseFinePM[lvl + 1]->exchangeAndAddFiCoData(bufferFiCo);      

      // d) Coarsen data from the refined grid to this grid.
      m_coarseFinePM[lvl + 1]->restrictAndAddFiCoDataToCoar(*a_meshData[lvl],
                                                            bufferFiCo,
                                                            EBCoarseFineParticleMesh::Average::Conservative);
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

    const EBISLayout& ebisl = eblg.getEBISL();
    const EBISLayout& ebislFiCo = eblgFiCo.getEBISL();    

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
        const EBISBox& ebisbox = ebislFiCo[din];
        
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
