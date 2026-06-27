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

template <typename P, typename Traits>
void
AmrMesh::allocate(ParticleContainer<P, Traits>& a_container, const std::string& a_realm) const
{
  CH_TIME("AmrMesh::allocate(ParticleContainer<P, Traits>, string)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(ParticleContainer<P, Traits>, string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::allocate(ParticleContainer<P, Traits>, string) - could not find Realm '" +
                            a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  a_container.define(m_realms[a_realm]->getGrids(),
                     m_realms[a_realm]->getDomains(),
                     m_realms[a_realm]->getDx(),
                     m_realms[a_realm]->getRefinementRatios(),
                     m_probLo,
                     m_minBlockSize,
                     m_finestLevel,
                     a_realm,
                     &m_realms[a_realm]->getValidCells()); // alias for remap()'s stayer fast path
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

template <typename P, typename Traits>
void
AmrMesh::remapToNewGrids(ParticleContainer<P, Traits>& a_particles,
                         const int                     a_lmin,
                         const int                     a_newFinestLevel) const noexcept
{
  CH_TIME("AmrMesh::remapToNewGrids(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::remapToNewGrids(ParticleContainer)" << endl;
  }

  (void)a_lmin; // ParticleContainer::regrid redistributes from the preRegrid cache; no lmin needed.

  const std::string realm = a_particles.getRealm();

  a_particles.regrid(this->getGrids(realm),
                     this->getDomains(),
                     this->getDx(),
                     this->getRefinementRatios(),
                     m_minBlockSize,
                     a_newFinestLevel);
}

template <typename P, typename Traits>
void
AmrMesh::depositParticles(EBAMRIVData&                        a_meshData,
                          const std::string&                  a_realm,
                          const phase::which_phase&           a_phase,
                          const ParticleContainer<P, Traits>& a_particles) const noexcept
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
AmrMesh::depositWeight(EBAMRCellData&                      a_meshData,
                       const std::string&                  a_realm,
                       const phase::which_phase&           a_phase,
                       const ParticleContainer<P, Traits>& a_particles,
                       const DepositionType                a_depositionType,
                       const CoarseFineDeposition          a_coarseFineDeposition,
                       const bool                          a_forceIrregNGP)
{
  CH_TIME("AmrMesh::depositWeight(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::depositWeight(ParticleContainer)" << endl;
  }

  EBAMRParticleMesh& particleMesh = this->getParticleMesh(a_realm, a_phase);

  particleMesh.depositWeight(a_meshData, a_particles, a_depositionType, a_coarseFineDeposition, a_forceIrregNGP);
}

template <typename P, typename Traits, typename CellGather>
void
AmrMesh::depositGathered(EBAMRCellData&                      a_meshData,
                         const std::string&                  a_realm,
                         const phase::which_phase&           a_phase,
                         const ParticleContainer<P, Traits>& a_particles,
                         const DepositionType                a_depositionType,
                         const CoarseFineDeposition          a_coarseFineDeposition,
                         const bool                          a_forceIrregNGP,
                         CellGather                          a_cellGather)
{
  CH_TIME("AmrMesh::depositGathered(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::depositGathered(ParticleContainer)" << endl;
  }

  EBAMRParticleMesh& particleMesh = this->getParticleMesh(a_realm, a_phase);

  particleMesh
    .depositGathered(a_meshData, a_particles, a_depositionType, a_coarseFineDeposition, a_forceIrregNGP, a_cellGather);
}

template <auto... Members, typename P, typename Traits>
void
AmrMesh::depositParticles(EBAMRCellData&                      a_meshData,
                          const std::string&                  a_realm,
                          const phase::which_phase&           a_phase,
                          const ParticleContainer<P, Traits>& a_particles,
                          const DepositionType                a_depositionType,
                          const CoarseFineDeposition          a_coarseFineDeposition,
                          const bool                          a_forceIrregNGP)
{
  CH_TIME("AmrMesh::depositParticles(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::depositParticles(ParticleContainer)" << endl;
  }

  EBAMRParticleMesh& particleMesh = this->getParticleMesh(a_realm, a_phase);

  particleMesh.template deposit<Members...>(a_meshData,
                                            a_particles,
                                            a_depositionType,
                                            a_coarseFineDeposition,
                                            a_forceIrregNGP);
}

template <typename P, typename Traits>
void
AmrMesh::interpolateWeight(ParticleContainer<P, Traits>& a_particles,
                           const std::string&            a_realm,
                           const phase::which_phase&     a_phase,
                           const EBAMRCellData&          a_meshScalarField,
                           const DepositionType          a_interpType,
                           const bool                    a_forceIrregNGP) const
{
  CH_TIME("AmrMesh::interpolateWeight(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::interpolateWeight(ParticleContainer)" << endl;
  }

  EBAMRParticleMesh& particleMesh = this->getParticleMesh(a_realm, a_phase);

  particleMesh.interpolateWeight(a_particles, a_meshScalarField, a_interpType, a_forceIrregNGP);
}

template <auto... Members, typename P, typename Traits>
void
AmrMesh::interpolateParticles(ParticleContainer<P, Traits>& a_particles,
                              const std::string&            a_realm,
                              const phase::which_phase&     a_phase,
                              const EBAMRCellData&          a_meshField,
                              const DepositionType          a_interpType,
                              const bool                    a_forceIrregNGP) const
{
  CH_TIME("AmrMesh::interpolateParticles(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::interpolateParticles(ParticleContainer)" << endl;
  }

  EBAMRParticleMesh& particleMesh = this->getParticleMesh(a_realm, a_phase);

  particleMesh.template interpolate<Members...>(a_particles, a_meshField, a_interpType, a_forceIrregNGP);
}

template <typename P, typename Traits>
void
AmrMesh::removeCoveredParticlesIF(ParticleContainer<P, Traits>& a_particles,
                                  const phase::which_phase&     a_phase,
                                  const Real                    a_tolerance) const
{
  CH_TIME("AmrMesh::removeCoveredParticlesIF(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::removeCoveredParticlesIF(ParticleContainer)" << endl;
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
    MayDay::Error("AmrMesh::removeCoveredParticlesIF(ParticleContainer) - logic bust");
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

template <typename P, typename Traits>
void
AmrMesh::removeCoveredParticlesDiscrete(ParticleContainer<P, Traits>& a_particles,
                                        const phase::which_phase&     a_phase,
                                        const Real                    a_tolerance) const
{
  CH_TIME("AmrMesh::removeCoveredParticlesDiscrete(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::removeCoveredParticlesDiscrete(ParticleContainer)" << endl;
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
AmrMesh::removeCoveredParticlesVoxels(ParticleContainer<P, Traits>& a_particles,
                                      const phase::which_phase&     a_phase) const
{
  CH_TIME("AmrMesh::removeCoveredParticlesVoxels(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::removeCoveredParticlesVoxels(ParticleContainer)" << endl;
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

template <typename P, typename Traits>
void
AmrMesh::transferCoveredParticlesIF(ParticleContainer<P, Traits>& a_particlesFrom,
                                    ParticleContainer<P, Traits>& a_particlesTo,
                                    const phase::which_phase&     a_phase,
                                    const Real                    a_tolerance) const
{
  CH_TIME("AmrMesh::transferCoveredParticlesIF(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferCoveredParticlesIF(ParticleContainer)" << endl;
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
    MayDay::Error("AmrMesh::transferCoveredParticlesIF(ParticleContainer) - logic bust");
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
AmrMesh::transferCoveredParticlesDiscrete(ParticleContainer<P, Traits>& a_particlesFrom,
                                          ParticleContainer<P, Traits>& a_particlesTo,
                                          const phase::which_phase&     a_phase,
                                          const Real                    a_tolerance) const
{
  CH_TIME("AmrMesh::transferCoveredParticlesDiscrete(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferCoveredParticlesDiscrete(ParticleContainer)" << endl;
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
AmrMesh::transferCoveredParticlesVoxels(ParticleContainer<P, Traits>& a_particlesFrom,
                                        ParticleContainer<P, Traits>& a_particlesTo,
                                        const phase::which_phase&     a_phase) const
{
  CH_TIME("AmrMesh::transferCoveredParticlesVoxels(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferCoveredParticlesVoxels(ParticleContainer)" << endl;
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

template <auto... OldPosition, typename P, typename Traits>
void
AmrMesh::intersectParticlesRaycastIF(
  ParticleContainer<P, Traits>&                                    a_activeParticles,
  ParticleContainer<P, Traits>&                                    a_ebParticles,
  ParticleContainer<P, Traits>&                                    a_domainParticles,
  const phase::which_phase                                         a_phase,
  const Real                                                       a_tolerance,
  const bool                                                       a_deleteParticles,
  const std::function<void(ParticleSoA<P, Traits>&, std::size_t)>& a_nonDeletionModifier) const noexcept
{
  CH_TIME("AmrMesh::intersectParticlesRaycastIF(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::intersectParticlesRaycastIF(ParticleContainer)" << endl;
  }

  static_assert(sizeof...(OldPosition) == SpaceDim,
                "AmrMesh::intersectParticlesRaycastIF(ParticleContainer) - need exactly SpaceDim oldPosition members");

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
    MayDay::Error("AmrMesh::intersectParticlesRaycastIF(ParticleContainer) - logic bust");

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
  ParticleContainer<P, Traits>&                                    a_activeParticles,
  ParticleContainer<P, Traits>&                                    a_ebParticles,
  ParticleContainer<P, Traits>&                                    a_domainParticles,
  const phase::which_phase                                         a_phase,
  const Real                                                       a_bisectionStep,
  const bool                                                       a_deleteParticles,
  const std::function<void(ParticleSoA<P, Traits>&, std::size_t)>& a_nonDeletionModifier) const noexcept
{
  CH_TIME("AmrMesh::intersectParticlesBisectIF(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::intersectParticlesBisectIF(ParticleContainer)" << endl;
  }

  static_assert(sizeof...(OldPosition) == SpaceDim,
                "AmrMesh::intersectParticlesBisectIF(ParticleContainer) - need exactly SpaceDim oldPosition members");

  CH_assert(a_activeParticles.getRealm() == a_ebParticles.getRealm());
  CH_assert(a_activeParticles.getRealm() == a_domainParticles.getRealm());

  // TLDR: This is pretty much a hard-copy of intersectParticlesRaycastIF(ParticleContainer), with the exception
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
    MayDay::Error("AmrMesh::intersectParticlesBisectIF(ParticleContainer) - logic bust");

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
AmrMesh::transferIrregularParticles(ParticleContainer<P, Traits>& a_dstParticles,
                                    ParticleContainer<P, Traits>& a_srcParticles,
                                    const phase::which_phase      a_phase) const noexcept
{
  CH_TIME("AmrMesh::transferIrregularParticles(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::transferIrregularParticles(ParticleContainer)" << endl;
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
