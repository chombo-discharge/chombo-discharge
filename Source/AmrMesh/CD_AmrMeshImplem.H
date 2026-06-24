/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_AmrMeshImplem.H
  @brief  Implementation of CD_AmrMesh.H
  @author Robert Marskar
*/

#ifndef CD_AMRMESHIMPLEM_H
#define CD_AMRMESHIMPLEM_H

// Our includes
#include <CD_AmrMesh.H>
#include <CD_ParticleOps.H>
#include <CD_NamespaceHeader.H>

template <typename T>
void
AmrMesh::copyData(EBAMRData<T>&       a_dst,
                  const EBAMRData<T>& a_src,
                  const CopyStrategy& a_toRegion,
                  const CopyStrategy& a_fromRegion) const noexcept
{
  CH_TIME("AmrMesh::copyData(EBAMRData, simple)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::copyData(EBAMRData, simple)" << endl;
  }

  const int      nComp    = a_dst[0]->nComp();
  const Interval dstComps = Interval(0, nComp - 1);
  const Interval srcComps = Interval(0, nComp - 1);

  const std::string toRealm   = a_dst.getRealm();
  const std::string fromRealm = a_src.getRealm();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    CH_assert(!(a_dst[lvl].isNull()));
    CH_assert(!(a_src[lvl].isNull()));

    this->copyData(*a_dst[lvl], *a_src[lvl], lvl, toRealm, fromRealm, dstComps, srcComps, a_toRegion, a_fromRegion);
  }
}

template <typename T>
void
AmrMesh::copyData(EBAMRData<T>&       a_dst,
                  const EBAMRData<T>& a_src,
                  const Interval&     a_dstComps,
                  const Interval&     a_srcComps,
                  const CopyStrategy& a_toRegion,
                  const CopyStrategy& a_fromRegion) const noexcept
{
  CH_TIME("AmrMesh::copyData(EBAMRData, full)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::copyData(EBAMRData, full)" << endl;
  }

  const std::string toRealm   = a_dst.getRealm();
  const std::string fromRealm = a_src.getRealm();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    CH_assert(!(a_dst[lvl].isNull()));
    CH_assert(!(a_src[lvl].isNull()));

    this->copyData(*a_dst[lvl], *a_src[lvl], lvl, toRealm, fromRealm, a_dstComps, a_srcComps, a_toRegion, a_fromRegion);
  }
}

template <typename T>
void
AmrMesh::copyData(LevelData<T>&       a_dst,
                  const LevelData<T>& a_src,
                  const int           a_level,
                  const std::string   a_toRealm,
                  const std::string   a_fromRealm,
                  const CopyStrategy& a_toRegion,
                  const CopyStrategy& a_fromRegion) const noexcept
{
  CH_TIME("AmrMesh::copyData(LD<T>_full)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::copyData(LD<T>_full)" << endl;
  }

  const Interval dstComps = Interval(0, a_dst.nComp() - 1);
  const Interval srcComps = Interval(0, a_src.nComp() - 1);

  this->copyData(a_dst, a_src, a_level, a_toRealm, a_fromRealm, dstComps, srcComps, a_toRegion, a_fromRegion);
}

template <typename T>
void
AmrMesh::copyData(LevelData<T>&       a_dst,
                  const LevelData<T>& a_src,
                  const int           a_level,
                  const std::string   a_toRealm,
                  const std::string   a_fromRealm,
                  const Interval&     a_dstComps,
                  const Interval&     a_srcComps,
                  const CopyStrategy& a_toRegion,
                  const CopyStrategy& a_fromRegion) const noexcept
{
  CH_TIME("AmrMesh::copyData(LD<T>_full)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::copyData(LD<T>_full)" << endl;
  }

  CH_assert(a_dstComps.size() == a_srcComps.size());
  CH_assert(a_dst.nComp() > a_dstComps.end());
  CH_assert(a_src.nComp() > a_srcComps.end());

  if (a_toRealm != a_fromRealm) {
    const auto id = std::make_pair(a_fromRealm, a_toRealm);

    if (a_toRegion == CopyStrategy::ValidGhost) {
      CH_assert(a_dst.ghostVect() == m_numGhostCells * IntVect::Unit);
    }
    if (a_fromRegion == CopyStrategy::ValidGhost) {
      CH_assert(a_src.ghostVect() == m_numGhostCells * IntVect::Unit);
    }

    Copier copier;

    if (a_fromRegion == CopyStrategy::Valid) {
      if (a_toRegion == CopyStrategy::Valid) {
        copier = m_validToValidRealmCopiers.at(id)[a_level];
      }
      else if (a_toRegion == CopyStrategy::ValidGhost) {
        copier = m_validToValidGhostRealmCopiers.at(id)[a_level];
      }
      else {
        MayDay::Abort("AmrMesh::copyData - logic bust 1");
      }
    }
    else if (a_fromRegion == CopyStrategy::ValidGhost) {
      if (a_toRegion == CopyStrategy::Valid) {
        copier = m_validGhostToValidRealmCopiers.at(id)[a_level];
      }
      else if (a_toRegion == CopyStrategy::ValidGhost) {
        copier = m_validGhostToValidGhostRealmCopiers.at(id)[a_level];
      }
      else {
        MayDay::Abort("AmrMesh::copyData - logic bust 2");
      }
    }
    else {
      MayDay::Abort("AmrMesh::copyData - logic bust 3");
    }

    a_src.copyTo(a_srcComps, a_dst, a_dstComps, copier);
  }
  else {
    a_src.localCopyTo(a_srcComps, a_dst, a_dstComps);
  }
}

template <typename T>
void
AmrMesh::deallocate(Vector<T*>& a_data) const
{
  CH_TIME("AmrMesh::deallocate(Vector<T*>)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::deallocate(Vector<T*>)" << endl;
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    delete a_data[lvl];
  }
}

template <typename T>
void
AmrMesh::deallocate(Vector<RefCountedPtr<T>>& a_data) const
{
  CH_TIME("AmrMesh::deallocate(Vector<RefCountedPtr<T>)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::deallocate(Vector<RefCountedPtr<T>)" << endl;
  }

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    //    delete &a_data[lvl];
    a_data[lvl] = RefCountedPtr<T>(0);
#if 0
    if(!a_data[lvl].isNull()){
      delete &(*a_data[lvl]);
      a_data[lvl] = RefCountedPtr<T> (NULL);
    }

#endif
  }
}

template <typename T>
void
AmrMesh::deallocate(EBAMRData<T>& a_data) const
{
  CH_TIME("AmrMesh::deallocate(EBAMRData<T>)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::deallocate(EBAMRData<T>)" << endl;
  }

  return this->deallocate(a_data.getData());
}

template <typename T>
void
AmrMesh::alias(Vector<T*>& a_alias, const Vector<RefCountedPtr<T>>& a_data) const
{
  CH_TIME("AmrMesh::alias(Vector<T*>, Vector<RefCountedPtr<T>)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::alias(Vector<T*>, Vector<RefCountedPtr<T>)" << endl;
  }

  a_alias.resize(a_data.size());

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    a_alias[lvl] = &(*a_data[lvl]);
  }
}

template <typename T, typename S>
void
AmrMesh::alias(Vector<T*>& a_alias, const EBAMRData<S>& a_data) const
{
  CH_TIME("AmrMesh::alias(Vector<T*>, EBAMRData<S>)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::alias(Vector<T*>, EBAMRData<S>" << endl;
  }

  return this->alias(a_alias, a_data.getData());
}

template <typename T>
void
AmrMesh::allocate(Vector<RefCountedPtr<ParticleData<T>>>& a_particles, const std::string& a_realm) const
{
  CH_TIME("AmrMesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >, string)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >, string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string
      str = "AmrMesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >, string) - could not find Realm '" + a_realm +
            "'";
    MayDay::Abort(str.c_str());
  }

  if (m_maxBoxSize != m_blockingFactor) {
    MayDay::Abort(
      "AmrMesh::allocate(Vector<RefCountedPtr<ParticleData<T> > >, string) - only constant box sizes supported for particle methods");
  }

  a_particles.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl    = m_realms[a_realm]->getGrids()[lvl];
    const ProblemDomain&     domain = m_realms[a_realm]->getDomains()[lvl];

    const Real dx = m_dx[lvl];

    a_particles[lvl] = RefCountedPtr<ParticleData<T>>(
      new ParticleData<T>(dbl, domain, m_blockingFactor, dx * RealVect::Unit, m_probLo));
  }
}

template <typename T>
void
AmrMesh::allocate(ParticleContainer<T>& a_container, const std::string& a_realm) const
{
  CH_TIME("AmrMesh::allocate(ParticleContainer<T>, int, string)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(ParticleContainer<T>, int, string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::allocate(ParticleContainer<T>, int, string) - could not find Realm '" + a_realm +
                            "'";
    MayDay::Abort(str.c_str());
  }

  if (m_maxBoxSize != m_blockingFactor) {
    MayDay::Abort(
      "AmrMesh::allocate(ParticleContainer<T>, int, string) - only constant box sizes are supported for particle methods");
  }

  a_container.define(m_realms[a_realm]->getGrids(),
                     m_realms[a_realm]->getDomains(),
                     m_realms[a_realm]->getDx(),
                     m_realms[a_realm]->getRefinementRatios(),
                     m_realms[a_realm]->getValidCells(),
                     m_realms[a_realm]->getLevelTiles(),
                     m_probLo,
                     m_blockingFactor,
                     m_finestLevel,
                     a_realm);
}

template <typename P, typename Traits>
void
AmrMesh::allocate(ParticleContainerSoA<P, Traits>& a_container, const std::string& a_realm) const
{
  CH_TIME("AmrMesh::allocate(ParticleContainerSoA<P, Traits>, string)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(ParticleContainerSoA<P, Traits>, string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::allocate(ParticleContainerSoA<P, Traits>, string) - could not find Realm '" +
                            a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  if (m_maxBoxSize != m_blockingFactor) {
    MayDay::Abort(
      "AmrMesh::allocate(ParticleContainerSoA<P, Traits>, string) - only constant box sizes are supported for particle methods");
  }

  a_container.define(m_realms[a_realm]->getGrids(),
                     m_realms[a_realm]->getDomains(),
                     m_realms[a_realm]->getDx(),
                     m_realms[a_realm]->getRefinementRatios(),
                     m_probLo,
                     m_blockingFactor,
                     m_finestLevel,
                     a_realm);
}

template <typename T>
void
AmrMesh::allocatePointer(Vector<RefCountedPtr<T>>& a_data) const
{
  CH_TIME("AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >)" << endl;
  }

  a_data.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    a_data[lvl] = RefCountedPtr<T>(new T());
  }
}

template <typename T>
void
AmrMesh::allocatePointer(Vector<RefCountedPtr<T>>& a_data, const int a_finestLevel) const
{
  CH_TIME("AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocatePointer(Vector<RefCountedPtr<T> >, int)" << endl;
  }

  a_data.resize(1 + a_finestLevel);

  for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
    a_data[lvl] = RefCountedPtr<T>(new T());
  }
}

template <typename T>
void
AmrMesh::allocatePointer(EBAMRData<T>& a_data, const std::string& a_realm) const
{
  CH_TIME("AmrMesh::allocatePointer(EBAMRData<T>, std::string)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocatePointer(EBAMRData<T>, std::string)" << endl;
  }

  this->allocatePointer(a_data.getData());

  a_data.setRealm(a_realm);
}

template <typename T>
void
AmrMesh::allocatePointer(EBAMRData<T>& a_data, const std::string& a_realm, const int a_finestLevel) const
{
  CH_TIME("AmrMesh::allocatePointer(EBAMRData<T>, std::string, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocatePointer(EBAMRData<T>, std::string, int)" << endl;
  }

  this->allocatePointer(a_data.getData(), a_finestLevel);

  a_data.setRealm(a_realm);
}

template <class P>
void
AmrMesh::remapToNewGrids(ParticleContainer<P>& a_particles, const int a_lmin, const int a_newFinestLevel) const noexcept
{
  CH_TIME("AmrMesh::remapToNewGrids");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::remapToNewGrids" << endl;
  }

  const std::string realm = a_particles.getRealm();

  a_particles.regrid(this->getGrids(realm),
                     this->getDomains(),
                     this->getDx(),
                     this->getRefinementRatios(),
                     this->getValidCells(realm),
                     this->getLevelTiles(realm),
                     a_lmin,
                     a_newFinestLevel);
}

template <typename P, typename Traits>
void
AmrMesh::remapToNewGrids(ParticleContainerSoA<P, Traits>& a_particles,
                         const int                        a_lmin,
                         const int                        a_newFinestLevel) const noexcept
{
  CH_TIME("AmrMesh::remapToNewGrids(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::remapToNewGrids(ParticleContainerSoA)" << endl;
  }

  (void)a_lmin; // ParticleContainerSoA::regrid redistributes from the preRegrid cache; no lmin needed.

  const std::string realm = a_particles.getRealm();

  a_particles.regrid(this->getGrids(realm),
                     this->getDomains(),
                     this->getDx(),
                     this->getRefinementRatios(),
                     m_blockingFactor,
                     a_newFinestLevel);
}

template <class P, class Ret, Ret (P::*MemberFunc)() const>
void
AmrMesh::depositParticles(EBAMRCellData&              a_meshData,
                          const std::string&          a_realm,
                          const phase::which_phase&   a_phase,
                          const ParticleContainer<P>& a_particles,
                          const DepositionType        a_depositionType,
                          const CoarseFineDeposition  a_coarseFineDeposition,
                          const bool                  a_forceIrregNGP)
{
  CH_TIME("AmrMesh::depositParticles");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::depositParticles" << endl;
  }

  EBAMRParticleMesh& particleMesh = this->getParticleMesh(a_realm, a_phase);

  particleMesh.deposit<P, Ret, MemberFunc>(a_meshData,
                                           a_particles,
                                           a_depositionType,
                                           a_coarseFineDeposition,
                                           a_forceIrregNGP);
}

template <class P, typename Ret, Ret (P::*MemberFunc)() const>
void
AmrMesh::depositParticles(EBAMRIVData&                a_meshData,
                          const std::string&          a_realm,
                          const phase::which_phase&   a_phase,
                          const ParticleContainer<P>& a_particles) const noexcept
{
  CH_TIME("AmrMesh::depositParticles(surface)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::depositParticles(surface)" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == 1);
  CH_assert(a_meshData.getRealm() == a_particles.getRealm());

  EBAMRSurfaceDeposition& surfaceDeposition = this->getSurfaceDeposition(a_realm, a_phase);

  surfaceDeposition.deposit<P, Ret, MemberFunc>(a_meshData, a_particles);
}

template <typename P, typename Traits>
void
AmrMesh::depositParticles(EBAMRIVData&                           a_meshData,
                          const std::string&                     a_realm,
                          const phase::which_phase&              a_phase,
                          const ParticleContainerSoA<P, Traits>& a_particles) const noexcept
{
  CH_TIME("AmrMesh::depositParticles(surface, SoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::depositParticles(surface, SoA)" << endl;
  }

  CH_assert(a_meshData[0]->nComp() == 1);
  CH_assert(a_meshData.getRealm() == a_particles.getRealm());

  EBAMRSurfaceDeposition& surfaceDeposition = this->getSurfaceDeposition(a_realm, a_phase);

  surfaceDeposition.deposit<P, Traits>(a_meshData, a_particles);
}

template <typename P, typename Traits>
void
AmrMesh::depositWeight(EBAMRCellData&                         a_meshData,
                       const std::string&                     a_realm,
                       const phase::which_phase&              a_phase,
                       const ParticleContainerSoA<P, Traits>& a_particles,
                       const DepositionType                   a_depositionType,
                       const CoarseFineDeposition             a_coarseFineDeposition,
                       const bool                             a_forceIrregNGP)
{
  CH_TIME("AmrMesh::depositWeight(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::depositWeight(ParticleContainerSoA)" << endl;
  }

  EBAMRParticleMeshSoA& particleMesh = this->getParticleMeshSoA(a_realm, a_phase);

  particleMesh.depositWeight(a_meshData, a_particles, a_depositionType, a_coarseFineDeposition, a_forceIrregNGP);
}

template <typename P, typename Traits, typename CellGather>
void
AmrMesh::depositGathered(EBAMRCellData&                         a_meshData,
                         const std::string&                     a_realm,
                         const phase::which_phase&              a_phase,
                         const ParticleContainerSoA<P, Traits>& a_particles,
                         const DepositionType                   a_depositionType,
                         const CoarseFineDeposition             a_coarseFineDeposition,
                         const bool                             a_forceIrregNGP,
                         CellGather                             a_cellGather)
{
  CH_TIME("AmrMesh::depositGathered(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::depositGathered(ParticleContainerSoA)" << endl;
  }

  EBAMRParticleMeshSoA& particleMesh = this->getParticleMeshSoA(a_realm, a_phase);

  particleMesh
    .depositGathered(a_meshData, a_particles, a_depositionType, a_coarseFineDeposition, a_forceIrregNGP, a_cellGather);
}

template <auto... Members, typename P, typename Traits>
void
AmrMesh::depositParticles(EBAMRCellData&                         a_meshData,
                          const std::string&                     a_realm,
                          const phase::which_phase&              a_phase,
                          const ParticleContainerSoA<P, Traits>& a_particles,
                          const DepositionType                   a_depositionType,
                          const CoarseFineDeposition             a_coarseFineDeposition,
                          const bool                             a_forceIrregNGP)
{
  CH_TIME("AmrMesh::depositParticles(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::depositParticles(ParticleContainerSoA)" << endl;
  }

  EBAMRParticleMeshSoA& particleMesh = this->getParticleMeshSoA(a_realm, a_phase);

  particleMesh.template deposit<Members...>(a_meshData,
                                            a_particles,
                                            a_depositionType,
                                            a_coarseFineDeposition,
                                            a_forceIrregNGP);
}

template <class P, class Ret, Ret (P::*MemberFunc)()>
void
AmrMesh::interpolateParticles(ParticleContainer<P>&     a_particles,
                              const std::string&        a_realm,
                              const phase::which_phase& a_phase,
                              const EBAMRCellData&      a_meshScalarField,
                              const DepositionType      a_interpType,
                              const bool                a_forceIrregNGP) const
{
  CH_TIME("AmrMesh::interpolateParticles(scalar)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::interpolateParticles(scalar)" << endl;
  }

  CH_assert(a_realm == a_particles.getRealm());
  CH_assert(a_realm == a_meshScalarField.getRealm());

  EBAMRParticleMesh& particleMesh = this->getParticleMesh(a_realm, a_phase);

  particleMesh.interpolate<P, Ret, MemberFunc>(a_particles, a_meshScalarField, a_interpType, a_forceIrregNGP);
}

template <typename P, typename Traits>
void
AmrMesh::interpolateWeight(ParticleContainerSoA<P, Traits>& a_particles,
                           const std::string&               a_realm,
                           const phase::which_phase&        a_phase,
                           const EBAMRCellData&             a_meshScalarField,
                           const DepositionType             a_interpType,
                           const bool                       a_forceIrregNGP) const
{
  CH_TIME("AmrMesh::interpolateWeight(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::interpolateWeight(ParticleContainerSoA)" << endl;
  }

  EBAMRParticleMeshSoA& particleMesh = this->getParticleMeshSoA(a_realm, a_phase);

  particleMesh.interpolateWeight(a_particles, a_meshScalarField, a_interpType, a_forceIrregNGP);
}

template <auto... Members, typename P, typename Traits>
void
AmrMesh::interpolateParticles(ParticleContainerSoA<P, Traits>& a_particles,
                              const std::string&               a_realm,
                              const phase::which_phase&        a_phase,
                              const EBAMRCellData&             a_meshField,
                              const DepositionType             a_interpType,
                              const bool                       a_forceIrregNGP) const
{
  CH_TIME("AmrMesh::interpolateParticles(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::interpolateParticles(ParticleContainerSoA)" << endl;
  }

  EBAMRParticleMeshSoA& particleMesh = this->getParticleMeshSoA(a_realm, a_phase);

  particleMesh.template interpolate<Members...>(a_particles, a_meshField, a_interpType, a_forceIrregNGP);
}

template <class P>
void
AmrMesh::removeCoveredParticlesIF(ParticleContainer<P>&     a_particles,
                                  const phase::which_phase& a_phase,
                                  const Real                a_tolerance) const
{
  CH_TIME("AmrMesh::removeCoveredParticlesIF");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::removeCoveredParticlesIF" << endl;
  }

  CH_assert(!a_particles.isOrganizedByCell());

  // Figure out the implicit function
  RefCountedPtr<BaseIF> implicitFunction;

  switch (a_phase) {
  case phase::gas: {
    implicitFunction = m_baseif.at(phase::gas);

    break;
  }
  case phase::solid: {
    implicitFunction = m_baseif.at(phase::solid);

    break;
  }
  default: {
    MayDay::Error("AmrMesh::removeCoveredParticlesIF - logic bust");
  }
  }

  // Get the realm where the particles live.
  const std::string whichRealm = a_particles.getRealm();

  // Go through all particles and remove them if they are less than dx*a_tolerance away from the EB.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = this->getGrids(whichRealm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();
    const Real               dx  = this->getDx()[lvl];
    const Real               tol = a_tolerance * dx;

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& particles = a_particles[lvl][din].listItems();

      // Check if particles are outside the implicit function.
      for (ListIterator<P> lit(particles); lit.ok();) {
        const RealVect& pos = lit().position();

        const Real f = implicitFunction->value(pos);

        if (f > tol) {
          particles.remove(lit);
        }
        else {
          lit++;
        }
      }
    }
  }
}

template <typename P, typename Traits>
void
AmrMesh::removeCoveredParticlesIF(ParticleContainerSoA<P, Traits>& a_particles,
                                  const phase::which_phase&        a_phase,
                                  const Real                       a_tolerance) const
{
  CH_TIME("AmrMesh::removeCoveredParticlesIF(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::removeCoveredParticlesIF(ParticleContainerSoA)" << endl;
  }

  // Figure out the implicit function
  RefCountedPtr<BaseIF> implicitFunction;

  switch (a_phase) {
  case phase::gas: {
    implicitFunction = m_baseif.at(phase::gas);

    break;
  }
  case phase::solid: {
    implicitFunction = m_baseif.at(phase::solid);

    break;
  }
  default: {
    MayDay::Error("AmrMesh::removeCoveredParticlesIF(ParticleContainerSoA) - logic bust");
  }
  }

  const std::string whichRealm = a_particles.getRealm();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = this->getGrids(whichRealm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();
    const Real               dx  = this->getDx()[lvl];
    const Real               tol = a_tolerance * dx;

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

      // swap-and-pop removal: do NOT advance i after a removal (a new particle now sits in slot i).
      std::size_t i = 0;
      while (i < leaf.size()) {
        if (implicitFunction->value(leaf.position(i)) > tol) {
          leaf.remove(i);
        }
        else {
          i++;
        }
      }
    }
  }
}

template <class P>
void
AmrMesh::removeCoveredParticlesDiscrete(ParticleContainer<P>&     a_particles,
                                        const phase::which_phase& a_phase,
                                        const Real                a_tolerance) const
{
  CH_TIME("AmrMesh::removeCoveredParticlesDiscrete");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::removeCoveredParticlesDiscrete" << endl;
  }

  CH_assert(!a_particles.isOrganizedByCell());

  // Get the realm where the particles live.
  const std::string whichRealm = a_particles.getRealm();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = this->getGrids(whichRealm)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = this->getEBISLayout(whichRealm, a_phase)[lvl];
    const Real               dx    = this->getDx()[lvl];
    const Real               tol   = a_tolerance * dx;

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din     = dit[mybox];
      const EBISBox&   ebisBox = ebisl[din];
      const Box        region  = ebisBox.getRegion();

      const bool isRegular   = ebisBox.isAllRegular();
      const bool isCovered   = ebisBox.isAllCovered();
      const bool isIrregular = !isRegular && !isCovered;

      List<P>& particles = a_particles[lvl][din].listItems();

      if (isCovered) {
        particles.clear();
      }
      else if (isIrregular) {
        for (ListIterator<P> lit(particles); lit.ok();) {
          const P&        particle = lit();
          const RealVect& pos      = particle.position();

          // Get the cell where the particle lives.
          const IntVect iv = ParticleOps::getParticleCellIndex(pos, m_probLo, dx);

          CH_assert(region.contains(iv));

          if (ebisBox.isCovered(iv)) {
            particles.remove(lit);
          }
          else if (ebisBox.isIrregular(iv)) {

            // There could potentially be multiple degrees of freedom -- we check if the particle is inside at least one
            // of the VoFs in the cell.

            bool insideAtLeastOneVoF = false;

            // For each VoF, compute the projection of the particle position onth the EB face. Because the EB normal points outwards, we have a positive
            // value if the particle lies in the valid region.
            const std::vector<VolIndex> vofs = ebisBox.getVoFs(iv).stdVector();
            for (const auto& vof : vofs) {
              const RealVect ebNormal       = ebisBox.normal(vof);
              const RealVect ebCentroid     = m_probLo + Location::position(Location::Cell::Boundary, vof, ebisBox, dx);
              const Real     faceProjection = ebNormal.dotProduct(pos - ebCentroid);

              if (faceProjection >= -tol) {
                insideAtLeastOneVoF = true;

                break;
              }
            }

            if (!insideAtLeastOneVoF) {
              particles.remove(lit);
            }
            else {
              ++lit;
            }
          }
          else {
            ++lit;
          }
        }
      }
    }
  }
}

template <class P>
void
AmrMesh::removeCoveredParticlesVoxels(ParticleContainer<P>& a_particles, const phase::which_phase& a_phase) const
{
  CH_TIME("AmrMesh::removeCoveredParticlesVoxels");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::removeCoveredParticlesVoxels" << endl;
  }

  CH_assert(!a_particles.isOrganizedByCell());

  // Get the realm where the particles live.
  const std::string whichRealm = a_particles.getRealm();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = this->getGrids(whichRealm)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = this->getEBISLayout(whichRealm, a_phase)[lvl];
    const Real               dx    = this->getDx()[lvl];

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din     = dit[mybox];
      const EBISBox&   ebisBox = ebisl[din];
      const Box        region  = ebisBox.getRegion();

      const bool isRegular   = ebisBox.isAllRegular();
      const bool isCovered   = ebisBox.isAllCovered();
      const bool isIrregular = !isRegular && !isCovered;

      List<P>& particles = a_particles[lvl][din].listItems();

      if (isCovered) {
        particles.clear();
      }
      else if (isIrregular) {
        for (ListIterator<P> lit(particles); lit.ok();) {
          const P&        particle = lit();
          const RealVect& pos      = particle.position();

          // Get the cell where the particle lives.
          const RealVect rv = (pos - m_probLo) / dx;
          const IntVect  iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

          CH_assert(region.contains(iv));

          if (ebisBox.isCovered(iv)) {
            particles.remove(lit);
          }
          else {
            ++lit;
          }
        }
      }
    }
  }
}

template <typename P, typename Traits>
void
AmrMesh::removeCoveredParticlesDiscrete(ParticleContainerSoA<P, Traits>& a_particles,
                                        const phase::which_phase&        a_phase,
                                        const Real                       a_tolerance) const
{
  CH_TIME("AmrMesh::removeCoveredParticlesDiscrete(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::removeCoveredParticlesDiscrete(ParticleContainerSoA)" << endl;
  }

  const std::string whichRealm = a_particles.getRealm();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = this->getGrids(whichRealm)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = this->getEBISLayout(whichRealm, a_phase)[lvl];
    const Real               dx    = this->getDx()[lvl];
    const Real               tol   = a_tolerance * dx;

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din     = dit[mybox];
      const EBISBox&   ebisBox = ebisl[din];
      const Box        region  = ebisBox.getRegion();

      const bool isRegular   = ebisBox.isAllRegular();
      const bool isCovered   = ebisBox.isAllCovered();
      const bool isIrregular = !isRegular && !isCovered;

      ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

      if (isCovered) {
        leaf.clear();
      }
      else if (isIrregular) {
        // swap-and-pop removal: do NOT advance i after a removal.
        std::size_t i = 0;
        while (i < leaf.size()) {
          const RealVect pos = leaf.position(i);
          const IntVect  iv  = ParticleOps::getParticleCellIndex(pos, m_probLo, dx);

          CH_assert(region.contains(iv));

          if (ebisBox.isCovered(iv)) {
            leaf.remove(i);
          }
          else if (ebisBox.isIrregular(iv)) {
            // Inside the valid region if it is on the valid side of at least one of the cell's VoF faces.
            bool insideAtLeastOneVoF = false;

            const std::vector<VolIndex> vofs = ebisBox.getVoFs(iv).stdVector();
            for (const auto& vof : vofs) {
              const RealVect ebNormal       = ebisBox.normal(vof);
              const RealVect ebCentroid     = m_probLo + Location::position(Location::Cell::Boundary, vof, ebisBox, dx);
              const Real     faceProjection = ebNormal.dotProduct(pos - ebCentroid);

              if (faceProjection >= -tol) {
                insideAtLeastOneVoF = true;
                break;
              }
            }

            if (!insideAtLeastOneVoF) {
              leaf.remove(i);
            }
            else {
              i++;
            }
          }
          else {
            i++;
          }
        }
      }
    }
  }
}

template <typename P, typename Traits>
void
AmrMesh::removeCoveredParticlesVoxels(ParticleContainerSoA<P, Traits>& a_particles,
                                      const phase::which_phase&        a_phase) const
{
  CH_TIME("AmrMesh::removeCoveredParticlesVoxels(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::removeCoveredParticlesVoxels(ParticleContainerSoA)" << endl;
  }

  const std::string whichRealm = a_particles.getRealm();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = this->getGrids(whichRealm)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = this->getEBISLayout(whichRealm, a_phase)[lvl];
    const Real               dx    = this->getDx()[lvl];

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din     = dit[mybox];
      const EBISBox&   ebisBox = ebisl[din];
      const Box        region  = ebisBox.getRegion();

      const bool isRegular   = ebisBox.isAllRegular();
      const bool isCovered   = ebisBox.isAllCovered();
      const bool isIrregular = !isRegular && !isCovered;

      ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

      if (isCovered) {
        leaf.clear();
      }
      else if (isIrregular) {
        std::size_t i = 0;
        while (i < leaf.size()) {
          const RealVect pos = leaf.position(i);
          const RealVect rv  = (pos - m_probLo) / dx;
          const IntVect  iv  = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

          CH_assert(region.contains(iv));

          if (ebisBox.isCovered(iv)) {
            leaf.remove(i);
          }
          else {
            i++;
          }
        }
      }
    }
  }
}

template <class P>
void
AmrMesh::transferCoveredParticlesIF(ParticleContainer<P>&     a_particlesFrom,
                                    ParticleContainer<P>&     a_particlesTo,
                                    const phase::which_phase& a_phase,
                                    const Real                a_tolerance) const
{
  CH_TIME("AmrMesh::transferCoveredParticlesIF");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferCoveredParticlesIF" << endl;
  }

  CH_assert(!a_particlesFrom.isOrganizedByCell());
  CH_assert(!a_particlesTo.isOrganizedByCell());

  // Figure out the implicit function
  RefCountedPtr<BaseIF> implicitFunction;

  switch (a_phase) {
  case phase::gas: {
    implicitFunction = m_baseif.at(phase::gas);

    break;
  }
  case phase::solid: {
    implicitFunction = m_baseif.at(phase::solid);

    break;
  }
  default: {
    MayDay::Error("AmrMesh::removeCoveredParticlesIF - logic bust");

    break;
  }
  }

  // Get the realm where the particles live.
  const std::string realmFrom = a_particlesFrom.getRealm();
  const std::string realmTo   = a_particlesTo.getRealm();

  CH_assert(realmFrom == realmTo);

  // Go through all particles and remove them if they are less than dx*a_tolerance away from the EB.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = this->getGrids(realmFrom)[lvl];
    const DataIterator&      dit = dbl.dataIterator();
    const Real               dx  = this->getDx()[lvl];
    const Real               tol = a_tolerance * dx;

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& particlesFrom = a_particlesFrom[lvl][din].listItems();
      List<P>& particlesTo   = a_particlesTo[lvl][din].listItems();

      // Check if particles are outside the implicit function.
      for (ListIterator<P> lit(particlesFrom); lit.ok();) {
        const RealVect& pos = lit().position();

        const Real f = implicitFunction->value(pos);

        if (f > tol) {
          particlesTo.transfer(lit);
        }
        else {
          ++lit;
        }
      }
    }
  }
}

template <class P>
void
AmrMesh::transferCoveredParticlesDiscrete(ParticleContainer<P>&     a_particlesFrom,
                                          ParticleContainer<P>&     a_particlesTo,
                                          const phase::which_phase& a_phase,
                                          const Real                a_tolerance) const
{
  CH_TIME("AmrMesh::transferCoveredParticlesDiscrete");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferCoveredParticlesDiscrete" << endl;
  }

  CH_assert(!a_particlesFrom.isOrganizedByCell());
  CH_assert(!a_particlesTo.isOrganizedByCell());

  // Get the realm where the particles live.
  const std::string realmFrom = a_particlesFrom.getRealm();
  const std::string realmTo   = a_particlesTo.getRealm();

  CH_assert(realmFrom == realmTo);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = this->getGrids(realmFrom)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = this->getEBISLayout(realmFrom, a_phase)[lvl];
    const Real               dx    = this->getDx()[lvl];
    const Real               tol   = a_tolerance * dx;

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din     = dit[mybox];
      const EBISBox&   ebisBox = ebisl[din];
      const Box        region  = ebisBox.getRegion();

      const bool isRegular   = ebisBox.isAllRegular();
      const bool isCovered   = ebisBox.isAllCovered();
      const bool isIrregular = !isRegular && !isCovered;

      List<P>& particlesFrom = a_particlesFrom[lvl][din].listItems();
      List<P>& particlesTo   = a_particlesTo[lvl][din].listItems();

      if (isCovered) {
        particlesTo.catenate(particlesFrom);
      }
      else if (isIrregular) {
        for (ListIterator<P> lit(particlesFrom); lit.ok();) {
          const P&        particle = lit();
          const RealVect& pos      = particle.position();

          // Get the cell where the particle lives.
          const IntVect iv = ParticleOps::getParticleCellIndex(pos, m_probLo, dx);

          CH_assert(region.contains(iv));

          if (ebisBox.isCovered(iv)) {
            particlesTo.transfer(lit);
          }
          else if (ebisBox.isIrregular(iv)) {
            // There could potentially be multiple degrees of freedom -- we check if the particle is inside at least one of the VoFs in the cell.

            bool insideAtLeastOneVoF = false;

            // For each VoF, compute the projection of the particle position onth the EB face. Because the EB normal points outwards, we have a positive
            // value if the particle lies in the valid region.
            const std::vector<VolIndex> vofs = ebisBox.getVoFs(iv).stdVector();
            for (const auto& vof : vofs) {
              const RealVect ebCentroid     = m_probLo + Location::position(Location::Cell::Boundary, vof, ebisBox, dx);
              const RealVect ebNormal       = ebisBox.normal(vof);
              const Real     faceProjection = ebNormal.dotProduct(pos - ebCentroid);

              if (faceProjection >= -tol) {
                insideAtLeastOneVoF = true;

                break;
              }
            }

            if (!insideAtLeastOneVoF) {
              particlesTo.transfer(lit);
            }
            else {
              ++lit;
            }
          }
          else {
            ++lit;
          }
        }
      }
    }
  }
}

template <class P>
void
AmrMesh::transferCoveredParticlesVoxels(ParticleContainer<P>&     a_particlesFrom,
                                        ParticleContainer<P>&     a_particlesTo,
                                        const phase::which_phase& a_phase) const
{
  CH_TIME("AmrMesh::transferCoveredParticlesVoxels");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferCoveredParticlesVoxels" << endl;
  }

  CH_assert(!a_particlesFrom.isOrganizedByCell());
  CH_assert(!a_particlesTo.isOrganizedByCell());

  // Get the realm where the particles live.
  const std::string realmFrom = a_particlesFrom.getRealm();
  const std::string realmTo   = a_particlesTo.getRealm();

  CH_assert(realmFrom == realmTo);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = this->getGrids(realmFrom)[lvl];
    const DataIterator       dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = this->getEBISLayout(realmFrom, a_phase)[lvl];
    const Real               dx    = this->getDx()[lvl];

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din     = dit[mybox];
      const EBISBox&   ebisBox = ebisl[din];
      const Box        region  = ebisBox.getRegion();

      const bool isRegular   = ebisBox.isAllRegular();
      const bool isCovered   = ebisBox.isAllCovered();
      const bool isIrregular = !isRegular && !isCovered;

      List<P>& particlesFrom = a_particlesFrom[lvl][din].listItems();
      List<P>& particlesTo   = a_particlesTo[lvl][din].listItems();

      if (isCovered) {
        particlesTo.catenate(particlesFrom);
      }
      else if (isIrregular) {
        for (ListIterator<P> lit(particlesFrom); lit.ok();) {
          const P&        particle = lit();
          const RealVect& pos      = particle.position();

          // Get the cell where the particle lives.
          const IntVect iv = ParticleOps::getParticleCellIndex(pos, m_probLo, dx);

          CH_assert(region.contains(iv));

          if (ebisBox.isCovered(iv)) {
            particlesTo.transfer(lit);
          }
          else {
            ++lit;
          }
        }
      }
    }
  }
}

template <typename P, typename Traits>
void
AmrMesh::transferCoveredParticlesIF(ParticleContainerSoA<P, Traits>& a_particlesFrom,
                                    ParticleContainerSoA<P, Traits>& a_particlesTo,
                                    const phase::which_phase&        a_phase,
                                    const Real                       a_tolerance) const
{
  CH_TIME("AmrMesh::transferCoveredParticlesIF(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferCoveredParticlesIF(ParticleContainerSoA)" << endl;
  }

  // Figure out the implicit function
  RefCountedPtr<BaseIF> implicitFunction;

  switch (a_phase) {
  case phase::gas: {
    implicitFunction = m_baseif.at(phase::gas);

    break;
  }
  case phase::solid: {
    implicitFunction = m_baseif.at(phase::solid);

    break;
  }
  default: {
    MayDay::Error("AmrMesh::transferCoveredParticlesIF(ParticleContainerSoA) - logic bust");
  }
  }

  const std::string realmFrom = a_particlesFrom.getRealm();
  const std::string realmTo   = a_particlesTo.getRealm();

  CH_assert(realmFrom == realmTo);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = this->getGrids(realmFrom)[lvl];
    const DataIterator&      dit = dbl.dataIterator();
    const Real               dx  = this->getDx()[lvl];
    const Real               tol = a_tolerance * dx;

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P, Traits>& from = a_particlesFrom[lvl][din];
      ParticleSoA<P, Traits>& to   = a_particlesTo[lvl][din];

      // swap-and-pop transfer: do NOT advance i after a transfer (a new particle now sits in slot i).
      std::size_t i = 0;
      while (i < from.size()) {
        if (implicitFunction->value(from.position(i)) > tol) {
          to.appendParticle(from, i);
          from.remove(i);
        }
        else {
          i++;
        }
      }
    }
  }
}

template <typename P, typename Traits>
void
AmrMesh::transferCoveredParticlesDiscrete(ParticleContainerSoA<P, Traits>& a_particlesFrom,
                                          ParticleContainerSoA<P, Traits>& a_particlesTo,
                                          const phase::which_phase&        a_phase,
                                          const Real                       a_tolerance) const
{
  CH_TIME("AmrMesh::transferCoveredParticlesDiscrete(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferCoveredParticlesDiscrete(ParticleContainerSoA)" << endl;
  }

  const std::string realmFrom = a_particlesFrom.getRealm();
  const std::string realmTo   = a_particlesTo.getRealm();

  CH_assert(realmFrom == realmTo);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = this->getGrids(realmFrom)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = this->getEBISLayout(realmFrom, a_phase)[lvl];
    const Real               dx    = this->getDx()[lvl];
    const Real               tol   = a_tolerance * dx;

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din     = dit[mybox];
      const EBISBox&   ebisBox = ebisl[din];
      const Box        region  = ebisBox.getRegion();

      const bool isRegular   = ebisBox.isAllRegular();
      const bool isCovered   = ebisBox.isAllCovered();
      const bool isIrregular = !isRegular && !isCovered;

      ParticleSoA<P, Traits>& from = a_particlesFrom[lvl][din];
      ParticleSoA<P, Traits>& to   = a_particlesTo[lvl][din];

      if (isCovered) {
        to.catenate(from);
      }
      else if (isIrregular) {
        // swap-and-pop transfer: do NOT advance i after a transfer.
        std::size_t i = 0;
        while (i < from.size()) {
          const RealVect pos = from.position(i);
          const IntVect  iv  = ParticleOps::getParticleCellIndex(pos, m_probLo, dx);

          CH_assert(region.contains(iv));

          if (ebisBox.isCovered(iv)) {
            to.appendParticle(from, i);
            from.remove(i);
          }
          else if (ebisBox.isIrregular(iv)) {
            // Inside the valid region if it is on the valid side of at least one of the cell's VoF faces.
            bool insideAtLeastOneVoF = false;

            const std::vector<VolIndex> vofs = ebisBox.getVoFs(iv).stdVector();
            for (const auto& vof : vofs) {
              const RealVect ebCentroid     = m_probLo + Location::position(Location::Cell::Boundary, vof, ebisBox, dx);
              const RealVect ebNormal       = ebisBox.normal(vof);
              const Real     faceProjection = ebNormal.dotProduct(pos - ebCentroid);

              if (faceProjection >= -tol) {
                insideAtLeastOneVoF = true;
                break;
              }
            }

            if (!insideAtLeastOneVoF) {
              to.appendParticle(from, i);
              from.remove(i);
            }
            else {
              i++;
            }
          }
          else {
            i++;
          }
        }
      }
    }
  }
}

template <typename P, typename Traits>
void
AmrMesh::transferCoveredParticlesVoxels(ParticleContainerSoA<P, Traits>& a_particlesFrom,
                                        ParticleContainerSoA<P, Traits>& a_particlesTo,
                                        const phase::which_phase&        a_phase) const
{
  CH_TIME("AmrMesh::transferCoveredParticlesVoxels(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferCoveredParticlesVoxels(ParticleContainerSoA)" << endl;
  }

  const std::string realmFrom = a_particlesFrom.getRealm();
  const std::string realmTo   = a_particlesTo.getRealm();

  CH_assert(realmFrom == realmTo);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = this->getGrids(realmFrom)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = this->getEBISLayout(realmFrom, a_phase)[lvl];
    const Real               dx    = this->getDx()[lvl];

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din     = dit[mybox];
      const EBISBox&   ebisBox = ebisl[din];
      const Box        region  = ebisBox.getRegion();

      const bool isRegular   = ebisBox.isAllRegular();
      const bool isCovered   = ebisBox.isAllCovered();
      const bool isIrregular = !isRegular && !isCovered;

      ParticleSoA<P, Traits>& from = a_particlesFrom[lvl][din];
      ParticleSoA<P, Traits>& to   = a_particlesTo[lvl][din];

      if (isCovered) {
        to.catenate(from);
      }
      else if (isIrregular) {
        std::size_t i = 0;
        while (i < from.size()) {
          const RealVect pos = from.position(i);
          const IntVect  iv  = ParticleOps::getParticleCellIndex(pos, m_probLo, dx);

          CH_assert(region.contains(iv));

          if (ebisBox.isCovered(iv)) {
            to.appendParticle(from, i);
            from.remove(i);
          }
          else {
            i++;
          }
        }
      }
    }
  }
}

template <class P>
void
AmrMesh::transferIrregularParticles(ParticleContainer<P>&          a_dstParticles,
                                    ParticleContainer<P>&          a_srcParticles,
                                    const phase::which_phase       a_phase,
                                    const std::function<void(P&)>& a_transferModifier) const noexcept
{
  CH_TIME("AmrMesh::transferIrregularParticles");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferIrregularParticles" << endl;
  }

  CH_assert(a_dstParticles.getRealm() == a_srcParticles.getRealm());

  const std::string realm = a_dstParticles.getRealm();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl    = this->getGrids(realm)[lvl];
    const DataIterator&      dit    = dbl.dataIterator();
    const EBISLayout&        ebisl  = this->getEBISLayout(realm, a_phase)[lvl];
    const ProblemDomain&     domain = this->getDomains()[lvl];
    const RealVect           probLo = this->getProbLo();
    const Real               dx     = this->getDx()[lvl];

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box      cellBox = dbl[din];
      const EBISBox& ebisbox = ebisl[din];

      List<P>& dstParticles = a_dstParticles[lvl][din].listItems();
      List<P>& srcParticles = a_srcParticles[lvl][din].listItems();

      for (ListIterator<P> lit(srcParticles); lit.ok();) {
        P& p = lit();

        const IntVect iv = ParticleOps::getParticleCellIndex(p.position(), probLo, dx);

        if (!(cellBox.contains(iv))) {
          MayDay::Warning("CD_AmrMeshImplem.H in routine 'transferIrregularParticles' - particle not in box!");

          ++lit;
        }

        if (ebisbox.isIrregular(iv)) {
          bool insideAnyVoF = false;

          const Vector<VolIndex> allVoFs = ebisbox.getVoFs(iv);

          for (int i = 0; i < allVoFs.size(); i++) {
            const VolIndex& vof    = allVoFs[i];
            const RealVect  normal = ebisbox.normal(vof);
            const RealVect  ebPos  = probLo + Location::position(Location::Cell::Boundary, vof, ebisbox, dx);

            // Note: Can have normal = 0!
            if ((p.position() - ebPos).dotProduct(normal) >= 0.0) {
              insideAnyVoF = true;
            }
          }

          if (!insideAnyVoF) {
            const VolIndex& vof   = allVoFs[0];
            const RealVect  ebPos = probLo + Location::position(Location::Cell::Boundary, vof, ebisbox, dx);

            p.position() = ebPos;

            a_transferModifier(p);

            dstParticles.transfer(lit);
          }
          else {
            ++lit;
          }
        }
        else {
          ++lit;
        }
      }
    }
  }
}

template <class P>
void
AmrMesh::intersectParticlesRaycastIF(ParticleContainer<P>&          a_activeParticles,
                                     ParticleContainer<P>&          a_ebParticles,
                                     ParticleContainer<P>&          a_domainParticles,
                                     const phase::which_phase       a_phase,
                                     const Real                     a_tolerance,
                                     const bool                     a_deleteParticles,
                                     const std::function<void(P&)>& a_nonDeletionModifier) const noexcept
{
  CH_TIME("AmrMesh::intersectParticlesRaycastIF");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::intersectParticlesRaycastIF" << endl;
  }

  CH_assert(a_activeParticles.getRealm() == a_ebParticles.getRealm());
  CH_assert(a_activeParticles.getRealm() == a_domainParticles.getRealm());

  a_ebParticles.clearParticles();
  a_domainParticles.clearParticles();

  const std::string whichRealm = a_activeParticles.getRealm();

  // Figure out the implicit function
  RefCountedPtr<BaseIF> implicitFunction;

  switch (a_phase) {
  case phase::gas: {
    implicitFunction = m_baseif.at(phase::gas);

    break;
  }
  case phase::solid: {
    implicitFunction = m_baseif.at(phase::solid);

    break;
  }
  default: {
    MayDay::Error("AmrMesh::intersectParticlesRaycastIF - logic bust");

    break;
  }
  }

  // Safety factor to prevent particles falling off the domain if they intersect the high-side of the domain
  constexpr Real safety = 1.E-12;

  // Level loop -- go through each AMR level
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

    // Handle to various grid stuff.
    const DisjointBoxLayout& dbl = this->getGrids(whichRealm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& activeParticles = a_activeParticles[lvl][din].listItems();
      List<P>& ebParticles     = a_ebParticles[lvl][din].listItems();
      List<P>& domainParticles = a_domainParticles[lvl][din].listItems();

      for (ListIterator<P> lit(activeParticles); lit.ok();) {
        P& particle = lit();

        const RealVect newPos = particle.position();
        const RealVect oldPos = particle.oldPosition();
        const RealVect path   = newPos - oldPos;

        // Cheap initial tests that allow us to skip some intersections tests.
        bool checkEB     = false;
        bool checkDomain = false;

        if (!implicitFunction.isNull()) {
          checkEB = true;
        }
        for (int dir = 0; dir < SpaceDim; dir++) {
          const bool outsideLo = newPos[dir] < m_probLo[dir];
          const bool outsideHi = newPos[dir] > m_probHi[dir];

          if (outsideLo || outsideHi) {
            checkDomain = true;
          }
        }

        // Do the intersection tests.
        if (checkEB || checkDomain) {

          // These are the solution
          Real sDomain = std::numeric_limits<Real>::max();
          Real sEB     = std::numeric_limits<Real>::max();

          bool contactDomain = false;
          bool contactEB     = false;

          // Check if the particle intersected the domain. If it did, we compute sDomain such that the intersection
          // point with the domain is X = X0 + sDomain * (X1-X0) where X1=newPos and X0=oldPos
          if (checkDomain) {
            contactDomain = ParticleOps::domainIntersection(oldPos, newPos, m_probLo, m_probHi, sDomain);
          }

          // Check if the particle intersected the EB. If it did, we compute sEB such that the intersection
          // point with the domain is X = X0 + sEB * (X1-X0) where X1=newPos and X0=oldPos
          if (checkEB) {
            contactEB = ParticleOps::ebIntersectionRaycast(implicitFunction, oldPos, newPos, a_tolerance, sEB);
          }

          // Particle bumped into something.
          if (contactDomain || contactEB) {
            if (sEB <= sDomain) { // Crashed with EB "first".
              const RealVect intersectionPos = oldPos + sEB * path;

              // If we delete the original particles we can just transfer the one we have. Otherwise we create
              // a new particle.
              if (a_deleteParticles) {
                particle.position() = intersectionPos;

                ebParticles.transfer(lit);
              }
              else {
                P p = lit();

                p.position() = intersectionPos;

                ebParticles.add(p);

                a_nonDeletionModifier(particle);

                lit++;
              }
            }
            else {
              // Crashed with domain "first". Safety factor is to prevent particles falling off the high side of the domain
              const Real sSafety = std::max((Real)0.0, sDomain - safety);

              const RealVect intersectionPos = oldPos + sSafety * path;

              if (a_deleteParticles) {
                particle.position() = intersectionPos;

                domainParticles.transfer(lit);
              }
              else {
                P p = lit();

                p.position() = intersectionPos;

                domainParticles.add(p);

                a_nonDeletionModifier(particle);

                ++lit;
              }
            }
          }
          else {
            ++lit;
          }
        }
        else {
          ++lit;
        }
      }
    }
  }

  // These need to be remapped.
  a_ebParticles.remap();
  a_domainParticles.remap();
}

template <class P>
void
AmrMesh::intersectParticlesBisectIF(ParticleContainer<P>&          a_activeParticles,
                                    ParticleContainer<P>&          a_ebParticles,
                                    ParticleContainer<P>&          a_domainParticles,
                                    const phase::which_phase       a_phase,
                                    const Real                     a_bisectionStep,
                                    const bool                     a_deleteParticles,
                                    const std::function<void(P&)>& a_nonDeletionModifier) const noexcept
{
  CH_TIME("AmrMesh::intersectParticlesBisectIF");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::intersectParticlesBisectIF" << endl;
  }

  CH_assert(a_activeParticles.getRealm() == a_ebParticles.getRealm());
  CH_assert(a_activeParticles.getRealm() == a_domainParticles.getRealm());

  // TLDR: This is pretty much a hard-copy of intersectParticlesRaycastIF, with the exception that the EB intersection
  //       test is replaced by a bisection test.

  a_ebParticles.clearParticles();
  a_domainParticles.clearParticles();

  const std::string whichRealm = a_activeParticles.getRealm();

  // Figure out the implicit function
  RefCountedPtr<BaseIF> implicitFunction;

  switch (a_phase) {
  case phase::gas: {
    implicitFunction = m_baseif.at(phase::gas);

    break;
  }
  case phase::solid: {
    implicitFunction = m_baseif.at(phase::solid);

    break;
  }
  default: {
    MayDay::Error("AmrMesh::intersectParticlesBisectIF - logic bust");

    break;
  }
  }

  // Safety factor to prevent particles falling off the domain if they intersect the high-side of the domain
  constexpr Real safety = 1.E-12;

  // Level loop -- go through each AMR level
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

    // Handle to various grid stuff.
    const DisjointBoxLayout& dbl = this->getGrids(whichRealm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& activeParticles = a_activeParticles[lvl][din].listItems();
      List<P>& ebParticles     = a_ebParticles[lvl][din].listItems();
      List<P>& domainParticles = a_domainParticles[lvl][din].listItems();

      for (ListIterator<P> lit(activeParticles); lit.ok();) {
        P& particle = lit();

        const RealVect newPos = particle.position();
        const RealVect oldPos = particle.oldPosition();
        const RealVect path   = newPos - oldPos;

        // Cheap initial tests that allow us to skip some intersections tests.
        bool checkEB     = false;
        bool checkDomain = false;

        if (!implicitFunction.isNull()) {
          checkEB = true;
        }
        for (int dir = 0; dir < SpaceDim; dir++) {
          const bool outsideLo = newPos[dir] < m_probLo[dir];
          const bool outsideHi = newPos[dir] > m_probHi[dir];

          if (outsideLo || outsideHi) {
            checkDomain = true;
          }
        }

        // Do the intersection tests.
        if (checkEB || checkDomain) {

          // These are the solution
          Real sDomain = std::numeric_limits<Real>::max();
          Real sEB     = std::numeric_limits<Real>::max();

          bool contactDomain = false;
          bool contactEB     = false;

          // Check if the particle intersected the domain. If it did, we compute sDomain such that the intersection
          // point with the domain is X = X0 + sDomain * (X1-X0) where X1=newPos and X0=oldPos
          if (checkDomain) {
            contactDomain = ParticleOps::domainIntersection(oldPos, newPos, m_probLo, m_probHi, sDomain);
          }

          // Check if the particle intersected the EB. If it did, we compute sEB such that the intersection
          // point with the domain is X = X0 + sEB * (X1-X0) where X1=newPos and X0=oldPos
          if (checkEB) {
            contactEB = ParticleOps::ebIntersectionBisect(implicitFunction, oldPos, newPos, a_bisectionStep, sEB);
          }

          // Particle bumped into something.
          if (contactDomain || contactEB) {
            if (sEB <= sDomain) {
              // Crashed with EB "first".
              const RealVect intersectionPos = oldPos + sEB * path;

              if (a_deleteParticles) {
                particle.position() = intersectionPos;

                ebParticles.transfer(lit);
              }
              else {
                P p = particle;

                p.position() = intersectionPos;

                ebParticles.add(p);

                ++lit;

                a_nonDeletionModifier(particle);
              }
            }
            else {
              // Crashed with domain "first".
              const Real sSafety = std::max((Real)0.0, sDomain - safety);

              const RealVect intersectionPos = oldPos + sSafety * path;

              if (a_deleteParticles) {
                particle.position() = intersectionPos;

                domainParticles.transfer(lit);
              }
              else {
                P p = particle;

                p.position() = intersectionPos;

                domainParticles.add(p);

                ++lit;

                a_nonDeletionModifier(particle);
              }
            }
          }
          else {
            ++lit;
          }
        }
        else {
          ++lit;
        }
      }
    }
  }

  // These need to be remapped.
  a_ebParticles.remap();
  a_domainParticles.remap();
}

template <auto... OldPosition, typename P, typename Traits>
void
AmrMesh::intersectParticlesRaycastIF(
  ParticleContainerSoA<P, Traits>&                                 a_activeParticles,
  ParticleContainerSoA<P, Traits>&                                 a_ebParticles,
  ParticleContainerSoA<P, Traits>&                                 a_domainParticles,
  const phase::which_phase                                         a_phase,
  const Real                                                       a_tolerance,
  const bool                                                       a_deleteParticles,
  const std::function<void(ParticleSoA<P, Traits>&, std::size_t)>& a_nonDeletionModifier) const noexcept
{
  CH_TIME("AmrMesh::intersectParticlesRaycastIF(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::intersectParticlesRaycastIF(ParticleContainerSoA)" << endl;
  }

  static_assert(
    sizeof...(OldPosition) == SpaceDim,
    "AmrMesh::intersectParticlesRaycastIF(ParticleContainerSoA) - need exactly SpaceDim oldPosition members");

  CH_assert(a_activeParticles.getRealm() == a_ebParticles.getRealm());
  CH_assert(a_activeParticles.getRealm() == a_domainParticles.getRealm());

  a_ebParticles.clearParticles();
  a_domainParticles.clearParticles();

  const std::string whichRealm = a_activeParticles.getRealm();

  // Figure out the implicit function
  RefCountedPtr<BaseIF> implicitFunction;

  switch (a_phase) {
  case phase::gas: {
    implicitFunction = m_baseif.at(phase::gas);

    break;
  }
  case phase::solid: {
    implicitFunction = m_baseif.at(phase::solid);

    break;
  }
  default: {
    MayDay::Error("AmrMesh::intersectParticlesRaycastIF(ParticleContainerSoA) - logic bust");

    break;
  }
  }

  // Safety factor to prevent particles falling off the domain if they intersect the high-side of the domain
  constexpr Real safety = 1.E-12;

  // Level loop -- go through each AMR level
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

    // Handle to various grid stuff.
    const DisjointBoxLayout& dbl = this->getGrids(whichRealm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P, Traits>& active = a_activeParticles[lvl][din];
      ParticleSoA<P, Traits>& eb     = a_ebParticles[lvl][din];
      ParticleSoA<P, Traits>& domain = a_domainParticles[lvl][din];

      // swap-and-pop transfer: do NOT advance i after a deletion (a new particle now sits in slot i).
      std::size_t i = 0;
      while (i < active.size()) {
        const RealVect newPos = active.position(i);

        // Assemble the previous position from the payload columns (left-to-right pack expansion).
        RealVect oldPos;
        int      odir = 0;
        ((oldPos[odir++] = active.template get<OldPosition>(i)), ...);

        const RealVect path = newPos - oldPos;

        // Cheap initial tests that allow us to skip some intersections tests.
        bool checkEB     = false;
        bool checkDomain = false;

        if (!implicitFunction.isNull()) {
          checkEB = true;
        }
        for (int dir = 0; dir < SpaceDim; dir++) {
          const bool outsideLo = newPos[dir] < m_probLo[dir];
          const bool outsideHi = newPos[dir] > m_probHi[dir];

          if (outsideLo || outsideHi) {
            checkDomain = true;
          }
        }

        // Do the intersection tests.
        if (checkEB || checkDomain) {

          // These are the solution
          Real sDomain = std::numeric_limits<Real>::max();
          Real sEB     = std::numeric_limits<Real>::max();

          bool contactDomain = false;
          bool contactEB     = false;

          if (checkDomain) {
            contactDomain = ParticleOps::domainIntersection(oldPos, newPos, m_probLo, m_probHi, sDomain);
          }
          if (checkEB) {
            contactEB = ParticleOps::ebIntersectionRaycast(implicitFunction, oldPos, newPos, a_tolerance, sEB);
          }

          // Particle bumped into something. Move it onto the intersection point and transfer it to the EB/domain
          // container. When the original is kept we restore its end position after copying the intersection point.
          if (contactDomain || contactEB) {
            if (sEB <= sDomain) { // Crashed with EB "first".
              const RealVect intersectionPos = oldPos + sEB * path;

              active.setPosition(i, intersectionPos);
              eb.appendParticle(active, i);

              if (a_deleteParticles) {
                active.remove(i);
              }
              else {
                active.setPosition(i, newPos);
                a_nonDeletionModifier(active, i);
                i++;
              }
            }
            else { // Crashed with domain "first".
              const Real     sSafety         = std::max((Real)0.0, sDomain - safety);
              const RealVect intersectionPos = oldPos + sSafety * path;

              active.setPosition(i, intersectionPos);
              domain.appendParticle(active, i);

              if (a_deleteParticles) {
                active.remove(i);
              }
              else {
                active.setPosition(i, newPos);
                a_nonDeletionModifier(active, i);
                i++;
              }
            }
          }
          else {
            i++;
          }
        }
        else {
          i++;
        }
      }
    }
  }

  // These need to be remapped.
  a_ebParticles.remap();
  a_domainParticles.remap();
}

template <auto... OldPosition, typename P, typename Traits>
void
AmrMesh::intersectParticlesBisectIF(
  ParticleContainerSoA<P, Traits>&                                 a_activeParticles,
  ParticleContainerSoA<P, Traits>&                                 a_ebParticles,
  ParticleContainerSoA<P, Traits>&                                 a_domainParticles,
  const phase::which_phase                                         a_phase,
  const Real                                                       a_bisectionStep,
  const bool                                                       a_deleteParticles,
  const std::function<void(ParticleSoA<P, Traits>&, std::size_t)>& a_nonDeletionModifier) const noexcept
{
  CH_TIME("AmrMesh::intersectParticlesBisectIF(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::intersectParticlesBisectIF(ParticleContainerSoA)" << endl;
  }

  static_assert(
    sizeof...(OldPosition) == SpaceDim,
    "AmrMesh::intersectParticlesBisectIF(ParticleContainerSoA) - need exactly SpaceDim oldPosition members");

  CH_assert(a_activeParticles.getRealm() == a_ebParticles.getRealm());
  CH_assert(a_activeParticles.getRealm() == a_domainParticles.getRealm());

  // TLDR: This is pretty much a hard-copy of intersectParticlesRaycastIF(ParticleContainerSoA), with the exception
  //       that the EB intersection test is replaced by a bisection test.

  a_ebParticles.clearParticles();
  a_domainParticles.clearParticles();

  const std::string whichRealm = a_activeParticles.getRealm();

  // Figure out the implicit function
  RefCountedPtr<BaseIF> implicitFunction;

  switch (a_phase) {
  case phase::gas: {
    implicitFunction = m_baseif.at(phase::gas);

    break;
  }
  case phase::solid: {
    implicitFunction = m_baseif.at(phase::solid);

    break;
  }
  default: {
    MayDay::Error("AmrMesh::intersectParticlesBisectIF(ParticleContainerSoA) - logic bust");

    break;
  }
  }

  // Safety factor to prevent particles falling off the domain if they intersect the high-side of the domain
  constexpr Real safety = 1.E-12;

  // Level loop -- go through each AMR level
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

    // Handle to various grid stuff.
    const DisjointBoxLayout& dbl = this->getGrids(whichRealm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P, Traits>& active = a_activeParticles[lvl][din];
      ParticleSoA<P, Traits>& eb     = a_ebParticles[lvl][din];
      ParticleSoA<P, Traits>& domain = a_domainParticles[lvl][din];

      // swap-and-pop transfer: do NOT advance i after a deletion (a new particle now sits in slot i).
      std::size_t i = 0;
      while (i < active.size()) {
        const RealVect newPos = active.position(i);

        // Assemble the previous position from the payload columns (left-to-right pack expansion).
        RealVect oldPos;
        int      odir = 0;
        ((oldPos[odir++] = active.template get<OldPosition>(i)), ...);

        const RealVect path = newPos - oldPos;

        // Cheap initial tests that allow us to skip some intersections tests.
        bool checkEB     = false;
        bool checkDomain = false;

        if (!implicitFunction.isNull()) {
          checkEB = true;
        }
        for (int dir = 0; dir < SpaceDim; dir++) {
          const bool outsideLo = newPos[dir] < m_probLo[dir];
          const bool outsideHi = newPos[dir] > m_probHi[dir];

          if (outsideLo || outsideHi) {
            checkDomain = true;
          }
        }

        // Do the intersection tests.
        if (checkEB || checkDomain) {

          // These are the solution
          Real sDomain = std::numeric_limits<Real>::max();
          Real sEB     = std::numeric_limits<Real>::max();

          bool contactDomain = false;
          bool contactEB     = false;

          if (checkDomain) {
            contactDomain = ParticleOps::domainIntersection(oldPos, newPos, m_probLo, m_probHi, sDomain);
          }
          if (checkEB) {
            contactEB = ParticleOps::ebIntersectionBisect(implicitFunction, oldPos, newPos, a_bisectionStep, sEB);
          }

          // Particle bumped into something. Move it onto the intersection point and transfer it to the EB/domain
          // container. When the original is kept we restore its end position after copying the intersection point.
          if (contactDomain || contactEB) {
            if (sEB <= sDomain) { // Crashed with EB "first".
              const RealVect intersectionPos = oldPos + sEB * path;

              active.setPosition(i, intersectionPos);
              eb.appendParticle(active, i);

              if (a_deleteParticles) {
                active.remove(i);
              }
              else {
                active.setPosition(i, newPos);
                a_nonDeletionModifier(active, i);
                i++;
              }
            }
            else { // Crashed with domain "first".
              const Real     sSafety         = std::max((Real)0.0, sDomain - safety);
              const RealVect intersectionPos = oldPos + sSafety * path;

              active.setPosition(i, intersectionPos);
              domain.appendParticle(active, i);

              if (a_deleteParticles) {
                active.remove(i);
              }
              else {
                active.setPosition(i, newPos);
                a_nonDeletionModifier(active, i);
                i++;
              }
            }
          }
          else {
            i++;
          }
        }
        else {
          i++;
        }
      }
    }
  }

  // These need to be remapped.
  a_ebParticles.remap();
  a_domainParticles.remap();
}

template <typename P, typename Traits>
void
AmrMesh::transferIrregularParticles(ParticleContainerSoA<P, Traits>& a_dstParticles,
                                    ParticleContainerSoA<P, Traits>& a_srcParticles,
                                    const phase::which_phase         a_phase) const noexcept
{
  CH_TIME("AmrMesh::transferIrregularParticles(ParticleContainerSoA)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferIrregularParticles(ParticleContainerSoA)" << endl;
  }

  CH_assert(a_dstParticles.getRealm() == a_srcParticles.getRealm());

  const std::string realm = a_dstParticles.getRealm();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = this->getGrids(realm)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = this->getEBISLayout(realm, a_phase)[lvl];
    const Real               dx    = m_dx[lvl];

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din     = dit[mybox];
      const Box        cellBox = dbl[din];
      const EBISBox&   ebisbox = ebisl[din];

      ParticleSoA<P, Traits>& dst = a_dstParticles[lvl][din];
      ParticleSoA<P, Traits>& src = a_srcParticles[lvl][din];

      // swap-and-pop transfer: do NOT advance i after a transfer (a new particle now sits in slot i).
      std::size_t i = 0;
      while (i < src.size()) {
        const RealVect pos = src.position(i);
        const IntVect  iv  = ParticleOps::getParticleCellIndex(pos, m_probLo, dx);

        bool moved = false;

        if (cellBox.contains(iv) && ebisbox.isIrregular(iv)) {
          bool insideAnyVoF = false;

          const Vector<VolIndex> allVoFs = ebisbox.getVoFs(iv);
          for (int k = 0; k < allVoFs.size(); k++) {
            const RealVect normal = ebisbox.normal(allVoFs[k]);
            const RealVect ebPos  = m_probLo + Location::position(Location::Cell::Boundary, allVoFs[k], ebisbox, dx);

            if ((pos - ebPos).dotProduct(normal) >= 0.0) {
              insideAnyVoF = true;
            }
          }

          if (!insideAnyVoF) {
            const VolIndex& vof   = allVoFs[0];
            const RealVect  ebPos = m_probLo + Location::position(Location::Cell::Boundary, vof, ebisbox, dx);

            src.setPosition(i, ebPos);
            dst.appendParticle(src, i);
            src.remove(i);

            moved = true;
          }
        }

        if (!moved) {
          i++;
        }
      }
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
