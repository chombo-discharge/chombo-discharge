/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBStateRedistribution.cpp
 * @brief  Implementation of CD_EBStateRedistribution.H
 * @author Robert Marskar
 */

// Std includes
#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

// Chombo includes
#include <CH_Timer.H>
#include <BaseIVFactory.H>
#include <Copier.H>
#include <EBISBox.H>

// Our includes
#include <CD_EBStateRedistribution.H>
#include <CD_VofUtils.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int EBStateRedistribution::s_numGhost;
constexpr int EBStateRedistribution::s_numOffsets;
constexpr int EBStateRedistribution::s_trackerComps;
constexpr int EBStateRedistribution::s_outside;
constexpr int EBStateRedistribution::s_valid;
constexpr int EBStateRedistribution::s_covered;
constexpr int EBStateRedistribution::s_ghost;

EBStateRedistribution::EBStateRedistribution(
  const Vector<RefCountedPtr<EBLevelGrid>>&              a_eblg,
  const Vector<RefCountedPtr<EBLevelGrid>>&              a_eblgCoFi,
  const Vector<int>&                                     a_refinementRatios,
  const int                                              a_finestLevel,
  const int                                              a_order,
  const Real                                             a_targetVolFrac,
  const Vector<RefCountedPtr<EBFluxRedistribution>>&     a_fluxRedistribution,
  const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_coveredByFiner,
  const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_coarseFineGhost) noexcept
{
  CH_TIME("EBStateRedistribution::EBStateRedistribution");

  if (a_order != 0 && a_order != 2) {
    MayDay::Abort("EBStateRedistribution -- the order must be 0 (constant) or 2 (linear)");
  }
  if (!(a_targetVolFrac > 0.0 && a_targetVolFrac <= 1.0)) {
    MayDay::Abort("EBStateRedistribution -- the target volume fraction must be in (0, 1]");
  }

  m_eblg               = a_eblg;
  m_eblgCoFi           = a_eblgCoFi;
  m_coveredByFiner     = a_coveredByFiner;
  m_coarseFineGhost    = a_coarseFineGhost;
  m_refinementRatios   = a_refinementRatios;
  m_fluxRedistribution = a_fluxRedistribution;
  m_finestLevel        = a_finestLevel;
  m_order              = a_order;
  m_targetVolFrac      = a_targetVolFrac;

  if (m_finestLevel > 0) {
    if (m_fluxRedistribution.size() < 1 + m_finestLevel) {
      MayDay::Abort("EBStateRedistribution -- the flux redistribution operators are needed on all levels");
    }

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      if (m_fluxRedistribution[lvl].isNull()) {
        MayDay::Abort("EBStateRedistribution -- the flux redistribution operators are needed on all levels");
      }
    }
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    if (m_coveredByFiner.size() <= lvl || m_coarseFineGhost.size() <= lvl || m_coveredByFiner[lvl].isNull() ||
        m_coarseFineGhost[lvl].isNull()) {
      MayDay::Abort("EBStateRedistribution -- the cell masks are needed on all levels");
    }
    if (!(m_coveredByFiner[lvl]->ghostVect() >= s_numGhost * IntVect::Unit) ||
        !(m_coarseFineGhost[lvl]->ghostVect() >= s_numGhost * IntVect::Unit)) {
      MayDay::Abort("EBStateRedistribution -- the cell masks need s_numGhost ghost cells");
    }
  }

  m_tracker.resize(1 + m_finestLevel);
  m_numNbhd.resize(1 + m_finestLevel);
  m_alpha.resize(1 + m_finestLevel);
  m_nbhdVol.resize(1 + m_finestLevel);
  m_nbhdCentroid.resize(1 + m_finestLevel);
  m_connected.resize(1 + m_finestLevel);
  m_activeBox.resize(1 + m_finestLevel);
  m_nbhdAverage.resize(1 + m_finestLevel);
  m_slopes.resize(1 + m_finestLevel);
  m_syncMass.resize(1 + m_finestLevel);
  m_ghostMass.resize(1 + m_finestLevel);
  m_ghostCopier.resize(1 + m_finestLevel);
  m_coarseCopier.resize(1 + m_finestLevel);
  m_state.resize(1 + m_finestLevel);
  m_massDiff.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    this->defineLevel(lvl);
  }
}

EBStateRedistribution::~EBStateRedistribution() noexcept
{
  CH_TIME("EBStateRedistribution::~EBStateRedistribution");
}

IntVect
EBStateRedistribution::offset(const int a_index) noexcept
{
  CH_assert(a_index >= 0 && a_index < s_numOffsets);

  IntVect ret;

  int index = a_index;

  for (int dir = 0; dir < SpaceDim; dir++) {
    ret[dir] = (index % 3) - 1;
    index /= 3;
  }

  return ret;
}

Real
EBStateRedistribution::volumeFraction(const EBISBox&       a_ebisbox,
                                      const ProblemDomain& a_domain,
                                      const IntVect&       a_iv) noexcept
{
  if (!a_domain.contains(a_iv)) {
    return 0.0;
  }
  if (a_ebisbox.isRegular(a_iv)) {
    return 1.0;
  }
  if (a_ebisbox.isMultiValued(a_iv)) {
    return -1.0;
  }
  if (a_ebisbox.isIrregular(a_iv)) {
    return a_ebisbox.volFrac(VolIndex(a_iv, 0));
  }

  return 0.0;
}

RealVect
EBStateRedistribution::cellCentroid(const EBISBox&       a_ebisbox,
                                    const ProblemDomain& a_domain,
                                    const IntVect&       a_iv) noexcept
{
  if (a_domain.contains(a_iv) && a_ebisbox.isIrregular(a_iv) && !a_ebisbox.isMultiValued(a_iv)) {
    return a_ebisbox.centroid(VolIndex(a_iv, 0));
  }

  return RealVect::Zero;
}

int
EBStateRedistribution::cellType(const BaseFab<bool>& a_coveredByFiner,
                                const BaseFab<bool>& a_coarseFineGhost,
                                const ProblemDomain& a_domain,
                                const IntVect&       a_iv) noexcept
{
  if (a_coarseFineGhost(a_iv, 0)) {
    return s_ghost;
  }
  if (!a_domain.contains(a_iv)) {
    return s_outside;
  }

  return a_coveredByFiner(a_iv, 0) ? s_covered : s_valid;
}

int
EBStateRedistribution::offsetIndex(const IntVect& a_offset) noexcept
{
  int ret    = 0;
  int stride = 1;

  for (int dir = 0; dir < SpaceDim; dir++) {
    CH_assert(std::abs(a_offset[dir]) <= 1);

    ret += (a_offset[dir] + 1) * stride;
    stride *= 3;
  }

  return ret;
}

void
EBStateRedistribution::defineLevel(const int a_level) noexcept
{
  CH_TIME("EBStateRedistribution::defineLevel");

  if (m_eblg[a_level]->getEBISL().getGhost() < s_numGhost) {
    MayDay::Abort("EBStateRedistribution -- the EBISLayout needs at least two ghost cells");
  }

  const EBLevelGrid&       eblg   = *m_eblg[a_level];
  const DisjointBoxLayout& dbl    = eblg.getDBL();
  const EBISLayout&        ebisl  = eblg.getEBISL();
  const ProblemDomain&     domain = eblg.getDomain();
  const DataIterator&      dit    = dbl.dataIterator();
  const int                nbox   = dit.size();
  const IntVect            ghost  = s_numGhost * IntVect::Unit;

  // The domain, grown into the ghost cells across periodic boundaries.
  Box periodicDomain = domain.domainBox();

  for (int dir = 0; dir < SpaceDim; dir++) {
    if (domain.isPeriodic(dir)) {
      periodicDomain.grow(dir, s_numGhost);
    }
  }

  m_tracker[a_level] = RefCountedPtr<LevelData<BaseFab<int>>>(new LevelData<BaseFab<int>>(dbl, s_trackerComps, ghost));
  m_numNbhd[a_level] = RefCountedPtr<LevelData<BaseFab<int>>>(new LevelData<BaseFab<int>>(dbl, 1, ghost));
  m_alpha[a_level]   = RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(dbl, 2, ghost));
  m_nbhdVol[a_level] = RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(dbl, 1, ghost));
  m_nbhdCentroid[a_level] = RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(dbl, SpaceDim, ghost));
  m_connected[a_level]    = RefCountedPtr<LevelData<BaseFab<int>>>(new LevelData<BaseFab<int>>(dbl, 1, IntVect::Zero));
  m_activeBox[a_level]    = RefCountedPtr<LayoutData<bool>>(new LayoutData<bool>(dbl));
  m_nbhdAverage[a_level]  = RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(dbl, 1, ghost));
  m_slopes[a_level]       = RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(dbl, SpaceDim, ghost));
  m_syncMass[a_level]     = RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(dbl, 1, IntVect::Zero));
  m_state[a_level]        = RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(dbl, 1, ghost));

  if (a_level > 0) {
    const EBLevelGrid& eblgCoar = *m_eblg[a_level - 1];

    const DisjointBoxLayout& dblCoarsened = m_eblgCoFi[a_level - 1]->getDBL();

    m_ghostMass[a_level] = RefCountedPtr<LevelData<FArrayBox>>(
      new LevelData<FArrayBox>(dblCoarsened, 1, IntVect::Unit));

    m_ghostCopier[a_level].ghostDefine(dblCoarsened, eblgCoar.getDBL(), eblgCoar.getDomain(), IntVect::Unit);
    m_coarseCopier[a_level].define(eblgCoar.getDBL(), dblCoarsened, eblgCoar.getDomain(), IntVect::Unit);
  }

  LayoutData<IntVectSet> irregCells(dbl);

  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    irregCells[din] = ebisl[din].getIrregIVS(dbl[din]);
  }

  m_massDiff[a_level] = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
    new LevelData<BaseIVFAB<Real>>(dbl, 1, IntVect::Zero, BaseIVFactory<Real>(ebisl, irregCells)));

  LevelData<BaseFab<int>>& tracker      = *m_tracker[a_level];
  LevelData<BaseFab<int>>& numNbhd      = *m_numNbhd[a_level];
  LevelData<FArrayBox>&    alpha        = *m_alpha[a_level];
  LevelData<FArrayBox>&    nbhdVol      = *m_nbhdVol[a_level];
  LevelData<FArrayBox>&    nbhdCentroid = *m_nbhdCentroid[a_level];
  LevelData<BaseFab<int>>& connected    = *m_connected[a_level];
  LayoutData<bool>&        activeBox    = *m_activeBox[a_level];

  // Find the neighbours that the small cells merge with.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din     = dit[mybox];
    const EBISBox&   ebisbox = ebisl[din];

    const BaseFab<bool>& covered = (*m_coveredByFiner[a_level])[din];
    const BaseFab<bool>& cfGhost = (*m_coarseFineGhost[a_level])[din];

    auto kappa = [&](const IntVect& a_iv) -> Real {
      return volumeFraction(ebisbox, domain, a_iv);
    };
    auto type = [&](const IntVect& a_iv) -> int {
      return cellType(covered, cfGhost, domain, a_iv);
    };

    BaseFab<int>& track = tracker[din];

    track.setVal(0);

    std::vector<IntVect> offsets;

    auto kernel = [&](const IntVect& iv) -> void {
      const Real k = kappa(iv);

      if (type(iv) == s_valid && k > 0.0 && k < m_targetVolFrac) {

        // The neighbourhood reaches one cell and the slopes two cells. The cells that they use must be single-valued.
        const Box block = grow(Box(iv, iv), s_numGhost) & periodicDomain;

        for (BoxIterator bit(block); bit.ok(); ++bit) {
          const IntVect& c    = bit();
          const bool     near = Box(iv - IntVect::Unit, iv + IntVect::Unit).contains(c);

          if (kappa(c) < 0.0 && (type(c) == s_valid || (near && type(c) != s_outside))) {
            MayDay::Abort("EBStateRedistribution -- a multi-valued cell lies within two cells of a small cell");
          }
        }

        this->mergeNeighbours(offsets, iv, ebisbox.normal(VolIndex(iv, 0)), ebisbox, covered, cfGhost, domain);

        track(iv, 0) = offsets.size();

        for (int i = 0; i < offsets.size(); i++) {
          track(iv, 1 + i) = offsetIndex(offsets[i]);
        }
      }
    };

    BoxLoops::loop(dbl[din], kernel);
  }

  tracker.exchange();

  // The number of neighbourhoods that each cell is in, and alpha_1 of the merged neighbourhoods. This includes the
  // first layer of ghost cells, which can be in the neighbourhoods of valid cells but have no owner on this level.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const EBISBox&       ebisbox = ebisl[din];
    const BaseFab<bool>& covered = (*m_coveredByFiner[a_level])[din];
    const BaseFab<bool>& cfGhost = (*m_coarseFineGhost[a_level])[din];

    auto kappa = [&](const IntVect& a_iv) -> Real {
      return volumeFraction(ebisbox, domain, a_iv);
    };
    auto type = [&](const IntVect& a_iv) -> int {
      return cellType(covered, cfGhost, domain, a_iv);
    };

    const BaseFab<int>& track = tracker[din];

    BaseFab<int>& nrs = numNbhd[din];
    FArrayBox&    alp = alpha[din];

    nrs.setVal(1);
    alp.setVal(1.0);

    auto kernel = [&](const IntVect& iv) -> void {
      if (type(iv) == s_outside) {
        return;
      }

      for (int idx = 0; idx < s_numOffsets; idx++) {
        const IntVect o = offset(idx);
        const IntVect i = iv + o;

        if (o != IntVect::Zero) {
          for (int n = 1; n <= track(i, 0); n++) {
            if (offset(track(i, n)) == -o) {
              nrs(iv, 0) += 1;
            }
          }
        }
      }

      if (track(iv, 0) > 0) {
        Real volOfNeighbours = 0.0;

        for (int n = 1; n <= track(iv, 0); n++) {
          volOfNeighbours += kappa(iv + offset(track(iv, n)));
        }

        alp(iv, 1) = (m_targetVolFrac - kappa(iv)) / volOfNeighbours;
      }
    };

    BoxLoops::loop(grow(dbl[din], 1) & periodicDomain, kernel);
  }

  numNbhd.exchange();
  alpha.exchange();

  // alpha_0 of each cell.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din     = dit[mybox];
    const EBISBox&   ebisbox = ebisl[din];

    auto kappa = [&](const IntVect& a_iv) -> Real {
      return volumeFraction(ebisbox, domain, a_iv);
    };

    const BaseFab<int>& nrs   = numNbhd[din];
    const BaseFab<int>& track = tracker[din];

    FArrayBox& alp = alpha[din];

    auto kernel = [&](const IntVect& iv) -> void {
      if (kappa(iv) > 0.0) {
        Real a0 = 1.0;

        for (int idx = 0; idx < s_numOffsets; idx++) {
          const IntVect o = offset(idx);
          const IntVect i = iv + o;

          if (o != IntVect::Zero) {
            for (int n = 1; n <= track(i, 0); n++) {
              if (offset(track(i, n)) == -o) {
                a0 -= alp(i, 1) / nrs(iv, 0);
              }
            }
          }
        }

        alp(iv, 0) = a0;
      }
      else {
        alp(iv, 0) = 0.0;
      }
    };

    BoxLoops::loop(dbl[din], kernel);
  }

  alpha.exchange();

  // Volumes and centroids of the neighbourhoods, the connected neighbours of the merged cells, and the boxes where the
  // redistribution changes the state.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din     = dit[mybox];
    const EBISBox&   ebisbox = ebisl[din];

    const BaseFab<bool>& covered = (*m_coveredByFiner[a_level])[din];
    const BaseFab<bool>& cfGhost = (*m_coarseFineGhost[a_level])[din];

    auto kappa = [&](const IntVect& a_iv) -> Real {
      return volumeFraction(ebisbox, domain, a_iv);
    };
    auto cent = [&](const IntVect& a_iv) -> RealVect {
      return cellCentroid(ebisbox, domain, a_iv);
    };
    auto type = [&](const IntVect& a_iv) -> int {
      return cellType(covered, cfGhost, domain, a_iv);
    };

    const BaseFab<int>& nrs   = numNbhd[din];
    const FArrayBox&    alp   = alpha[din];
    const BaseFab<int>& track = tracker[din];

    FArrayBox&    vol     = nbhdVol[din];
    FArrayBox&    centHat = nbhdCentroid[din];
    BaseFab<int>& conn    = connected[din];

    vol.setVal(0.0);
    centHat.setVal(0.0);
    conn.setVal(0);

    bool active = false;

    auto kernel = [&](const IntVect& iv) -> void {
      if (kappa(iv) > 0.0) {
        vol(iv, 0) = alp(iv, 0) * kappa(iv);

        for (int dir = 0; dir < SpaceDim; dir++) {
          centHat(iv, dir) = cent(iv)[dir];
        }

        if (track(iv, 0) > 0) {
          RealVect c = RealVect::Zero;

          for (int dir = 0; dir < SpaceDim; dir++) {
            c[dir] = cent(iv)[dir] * alp(iv, 0) * kappa(iv);
          }

          for (int n = 1; n <= track(iv, 0); n++) {
            const IntVect o = offset(track(iv, n));
            const IntVect r = iv + o;
            const Real    w = alp(iv, 1) * kappa(r) / nrs(r, 0);

            vol(iv, 0) += w;

            for (int dir = 0; dir < SpaceDim; dir++) {
              c[dir] += (cent(r)[dir] + o[dir]) * w;
            }
          }

          for (int dir = 0; dir < SpaceDim; dir++) {
            centHat(iv, dir) = c[dir] / vol(iv, 0);
          }

          // The neighbours in the 3^D block that are connected to this cell, for the slopes. Cells that stand in for
          // another level are included.
          const Vector<VolIndex> neighbours = VofUtils::getVofsInRadius(VolIndex(iv, 0),
                                                                        ebisbox,
                                                                        1,
                                                                        VofUtils::Connectivity::MonotonePath,
                                                                        false);

          for (int i = 0; i < neighbours.size(); i++) {
            const IntVect& nb = neighbours[i].gridIndex();

            if (type(nb) != s_outside) {
              conn(iv, 0) |= (1 << offsetIndex(nb - iv));
            }
          }
        }

        if (track(iv, 0) > 0 || nrs(iv, 0) > 1.0) {
          active = true;
        }
      }
    };

    BoxLoops::loop(dbl[din], kernel);

    activeBox[din] = active;
  }

  nbhdVol.exchange();
  nbhdCentroid.exchange();

  // Ghost cells are in no merged neighbourhood, so their neighbourhood centroid is the cell centroid.
  if (a_level > 0) {
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&     din     = dit[mybox];
      const EBISBox&       ebisbox = ebisl[din];
      const BaseFab<bool>& covered = (*m_coveredByFiner[a_level])[din];
      const BaseFab<bool>& cfGhost = (*m_coarseFineGhost[a_level])[din];

      auto cent = [&](const IntVect& a_iv) -> RealVect {
        return cellCentroid(ebisbox, domain, a_iv);
      };
      auto type = [&](const IntVect& a_iv) -> int {
        return cellType(covered, cfGhost, domain, a_iv);
      };

      FArrayBox& centHat = nbhdCentroid[din];

      auto kernel = [&](const IntVect& iv) -> void {
        if (type(iv) == s_ghost) {
          for (int dir = 0; dir < SpaceDim; dir++) {
            centHat(iv, dir) = cent(iv)[dir];
          }
        }
      };

      BoxLoops::loop(grow(dbl[din], s_numGhost) & periodicDomain, kernel);
    }
  }
}

void
EBStateRedistribution::mergeNeighbours(std::vector<IntVect>& a_offsets,
                                       const IntVect&        a_iv,
                                       const RealVect&       a_normal,
                                       const EBISBox&        a_ebisbox,
                                       const BaseFab<bool>&  a_coveredByFiner,
                                       const BaseFab<bool>&  a_coarseFineGhost,
                                       const ProblemDomain&  a_domain) const noexcept
{
  // Two normal components are taken as equal within this tolerance, as in AMReX.
  constexpr Real smallNormDiff = 1.E-8;

  const Box& domainBox = a_domain.domainBox();

  a_offsets.resize(0);

  auto kappa = [&](const IntVect& a_cell) -> Real {
    return volumeFraction(a_ebisbox, a_domain, a_cell);
  };

  // Whether a cell can be in the neighbourhood: a valid cell, or a cell that stands in for another level.
  auto isUsable = [&](const IntVect& a_cell) -> bool {
    return cellType(a_coveredByFiner, a_coarseFineGhost, a_domain, a_cell) != s_outside;
  };

  bool openLo[SpaceDim];
  bool openHi[SpaceDim];

  for (int dir = 0; dir < SpaceDim; dir++) {
    openLo[dir] = (a_domain.isPeriodic(dir) || a_iv[dir] > domainBox.smallEnd(dir)) && isUsable(a_iv - BASISV(dir));
    openHi[dir] = (a_domain.isPeriodic(dir) || a_iv[dir] < domainBox.bigEnd(dir)) && isUsable(a_iv + BASISV(dir));
  }

  // The directions to merge in, in order, each with the side to merge on. The cell merges along the largest normal
  // component first (equal components: the larger direction index first, as in AMReX), on the side the normal points
  // to. A direction where that side is outside the domain is moved last and uses the other side.
  std::vector<std::pair<int, int>> mergeDirs;
  std::vector<std::pair<int, int>> blockedDirs;

  std::vector<int> byNormal(SpaceDim);

  for (int dir = 0; dir < SpaceDim; dir++) {
    byNormal[dir] = dir;
  }

  std::sort(byNormal.begin(), byNormal.end(), [&](const int a, const int b) -> bool {
    const Real na = std::abs(a_normal[a]);
    const Real nb = std::abs(a_normal[b]);

    return (na != nb) ? (na > nb) : (a > b);
  });

  for (const int dir : byNormal) {
    const int side = (a_normal[dir] > 0.0) ? 1 : -1;

    const bool open      = (side > 0) ? openHi[dir] : openLo[dir];
    const bool otherOpen = (side > 0) ? openLo[dir] : openHi[dir];

    if (open) {
      mergeDirs.emplace_back(dir, side);
    }
    else if (otherOpen) {
      blockedDirs.emplace_back(dir, -side);
    }
  }

  mergeDirs.insert(mergeDirs.end(), blockedDirs.begin(), blockedDirs.end());

  if (mergeDirs.size() == 0) {
    MayDay::Abort("EBStateRedistribution::mergeNeighbours -- a small cell has no neighbour inside the domain");
  }

  // If the first neighbour is covered, the cell instead merges with the cardinal neighbours with the largest volume
  // fractions, as in AMReX: for each direction the side with the larger volume fraction, and the directions in order of
  // decreasing volume fraction. The remaining directions follow in the order above.
  const IntVect firstOffset = mergeDirs[0].second * BASISV(mergeDirs[0].first);

  if (kappa(a_iv + firstOffset) == 0.0) {
    std::vector<std::pair<int, int>> byVolFrac;
    std::vector<Real>                bestVolFrac(SpaceDim, -1.0);

    for (int dir = 0; dir < SpaceDim; dir++) {
      int side = 0;

      if (openHi[dir]) {
        side             = 1;
        bestVolFrac[dir] = kappa(a_iv + BASISV(dir));
      }
      if (openLo[dir] && kappa(a_iv - BASISV(dir)) > bestVolFrac[dir]) {
        side             = -1;
        bestVolFrac[dir] = kappa(a_iv - BASISV(dir));
      }

      if (side != 0 && bestVolFrac[dir] > 0.0) {
        byVolFrac.emplace_back(dir, side);
      }
    }

    std::stable_sort(byVolFrac.begin(),
                     byVolFrac.end(),
                     [&](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                       return bestVolFrac[a.first] > bestVolFrac[b.first];
                     });

    if (byVolFrac.size() == 0) {
      MayDay::Abort("EBStateRedistribution::mergeNeighbours -- all cardinal neighbours of a small cell are covered");
    }

    for (const auto& d : mergeDirs) {
      const bool used = std::any_of(byVolFrac.begin(), byVolFrac.end(), [&](const std::pair<int, int>& u) {
        return u.first == d.first;
      });

      if (!used) {
        byVolFrac.push_back(d);
      }
    }

    mergeDirs = byVolFrac;
  }

  // Merge in one direction, then complete a 2x2 block with the next direction, and in 3D a 2x2x2 block with the last,
  // while the merged volume is below the target or the normal component of the next direction equals one already used
  // (merging in only one of two equal directions would break the symmetry). A direction is skipped if a cell that it
  // adds to the block is outside the domain.
  Real sumVol = kappa(a_iv);

  bool limited = false;

  std::vector<IntVect> block(1, IntVect::Zero);

  for (int m = 0; m < mergeDirs.size(); m++) {
    const int     dir  = mergeDirs[m].first;
    const IntVect step = mergeDirs[m].second * BASISV(dir);

    if (m > 0) {
      bool tie = false;

      for (int u = 0; u < m; u++) {
        if (std::abs(std::abs(a_normal[dir]) - std::abs(a_normal[mergeDirs[u].first])) < smallNormDiff) {
          tie = true;
        }
      }

      if (sumVol >= m_targetVolFrac && !tie) {
        break;
      }
    }

    const int numBlock = block.size();

    bool blockValid = true;

    for (int b = 0; b < numBlock; b++) {
      if (!isUsable(a_iv + block[b] + step)) {
        blockValid = false;
      }
    }

    if (!blockValid) {
      limited = true;

      continue;
    }

    for (int b = 0; b < numBlock; b++) {
      const IntVect o = block[b] + step;

      block.push_back(o);
      a_offsets.push_back(o);

      sumVol += kappa(a_iv + o);
    }
  }

  // Next to a domain boundary the neighbourhood can be smaller than the target. Elsewhere that is an error.
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (!openLo[dir] || !openHi[dir]) {
      limited = true;
    }
  }

  if (sumVol < m_targetVolFrac && !limited) {
    MayDay::Abort("EBStateRedistribution::mergeNeighbours -- could not merge a small cell up to the target volume");
  }
}

RealVect
EBStateRedistribution::computeSlope(const IntVect&       a_iv,
                                    const FArrayBox&     a_average,
                                    const FArrayBox&     a_centroid,
                                    const EBISBox&       a_ebisbox,
                                    const BaseFab<bool>& a_coveredByFiner,
                                    const BaseFab<bool>& a_coarseFineGhost,
                                    const ProblemDomain& a_domain,
                                    const int            a_connected) const noexcept
{
  constexpr Real epsilon       = 1.E-12;
  constexpr Real minSpreadSize = 0.5;

  RealVect slope = RealVect::Zero;

  if (m_order == 0) {
    return slope;
  }

  auto position = [&](const IntVect& a_o) -> RealVect {
    RealVect x;

    for (int dir = 0; dir < SpaceDim; dir++) {
      x[dir] = a_centroid(a_iv + a_o, dir) + a_o[dir] - a_centroid(a_iv, dir);
    }

    return x;
  };

  // If the connected neighbourhood centroids span less than half a cell in a direction, the stencil is widened to two
  // cells in that direction, as in AMReX.
  IntVect width = IntVect::Unit;

  {
    RealVect lo = RealVect::Zero;
    RealVect hi = RealVect::Zero;

    for (int idx = 0; idx < s_numOffsets; idx++) {
      if (a_connected & (1 << idx)) {
        const RealVect x = position(offset(idx));

        for (int dir = 0; dir < SpaceDim; dir++) {
          lo[dir] = std::min(lo[dir], x[dir]);
          hi[dir] = std::max(hi[dir], x[dir]);
        }
      }
    }

    for (int dir = 0; dir < SpaceDim; dir++) {
      if (hi[dir] - lo[dir] < minSpreadSize) {
        width[dir] = 2;
      }
    }
  }

  // Least-squares fit of the neighbourhood averages at the neighbourhood centroids, through the average of this cell.
  Real AtA[SpaceDim][SpaceDim];
  Real Atb[SpaceDim];

  for (int a = 0; a < SpaceDim; a++) {
    Atb[a] = 0.0;
    for (int b = 0; b < SpaceDim; b++) {
      AtA[a][b] = 0.0;
    }
  }

  const bool grown = (width != IntVect::Unit);
  const Box  stencil(-width, width);

  for (BoxIterator bit(stencil); bit.ok(); ++bit) {
    const IntVect o = bit();

    if (o == IntVect::Zero) {
      continue;
    }

    // The 3^D stencil uses the connected neighbours; the widened stencil the cells that are not covered by the EB.
    const bool use = grown ? (cellType(a_coveredByFiner, a_coarseFineGhost, a_domain, a_iv + o) != s_outside &&
                              volumeFraction(a_ebisbox, a_domain, a_iv + o) > 0.0)
                           : static_cast<bool>(a_connected & (1 << offsetIndex(o)));

    if (use) {
      const RealVect x  = position(o);
      const Real     du = a_average(a_iv + o, 0) - a_average(a_iv, 0);

      for (int a = 0; a < SpaceDim; a++) {
        Atb[a] += x[a] * du;
        for (int b = 0; b < SpaceDim; b++) {
          AtA[a][b] += x[a] * x[b];
        }
      }
    }
  }

  // Solve AtA*s = Atb by Gaussian elimination with partial pivoting. A singular system gives no slope.
  for (int col = 0; col < SpaceDim; col++) {
    int pivot = col;

    for (int row = col + 1; row < SpaceDim; row++) {
      if (std::abs(AtA[row][col]) > std::abs(AtA[pivot][col])) {
        pivot = row;
      }
    }

    if (std::abs(AtA[pivot][col]) < epsilon) {
      return RealVect::Zero;
    }

    if (pivot != col) {
      for (int b = 0; b < SpaceDim; b++) {
        std::swap(AtA[col][b], AtA[pivot][b]);
      }
      std::swap(Atb[col], Atb[pivot]);
    }

    for (int row = col + 1; row < SpaceDim; row++) {
      const Real f = AtA[row][col] / AtA[col][col];

      for (int b = col; b < SpaceDim; b++) {
        AtA[row][b] -= f * AtA[col][b];
      }
      Atb[row] -= f * Atb[col];
    }
  }

  for (int row = SpaceDim - 1; row >= 0; row--) {
    Real sum = Atb[row];

    for (int b = row + 1; b < SpaceDim; b++) {
      sum -= AtA[row][b] * slope[b];
    }

    slope[row] = sum / AtA[row][row];
  }

  // Limit the slope so that its values at the centroids of the connected neighbours lie between the averages of this
  // cell and of that neighbour, as AMReX does (amrex_calc_centroid_limiter).
  const Real q = a_average(a_iv, 0);

  RealVect limiter = RealVect::Unit;

  for (int idx = 0; idx < s_numOffsets; idx++) {
    const IntVect o = offset(idx);

    if (o != IntVect::Zero && (a_connected & (1 << idx))) {
      const RealVect x    = position(o);
      const Real     qHat = q + x.dotProduct(slope);
      const Real     qNb  = a_average(a_iv + o, 0);
      const Real     qMax = std::max(qNb, q);
      const Real     qMin = std::min(qNb, q);

      if (qHat - qMax > std::abs(epsilon * qMax) || qHat - qMin < -std::abs(epsilon * qMin)) {
        const Real sml = epsilon * std::max(std::abs(qMax), std::abs(qMin));

        Real lim = 1.0;

        if (qHat - q > sml) {
          lim = std::min(1.0, (qMax - q) / (qHat - q));
        }
        else if (qHat - q < -sml) {
          lim = std::min(1.0, (qMin - q) / (qHat - q));
        }

        for (int dir = 0; dir < SpaceDim; dir++) {
          if (std::abs(x[dir]) > epsilon) {
            limiter[dir] = std::min(limiter[dir], lim);
          }
        }
      }
    }
  }

  for (int dir = 0; dir < SpaceDim; dir++) {
    slope[dir] *= limiter[dir];
  }

  return slope;
}

void
EBStateRedistribution::redistribute(EBAMRCellData& a_divG, const EBAMRCellData& a_phi, const Real a_dt) noexcept
{
  CH_TIME("EBStateRedistribution::redistribute");

  if (a_dt == 0.0) {
    MayDay::Abort("EBStateRedistribution::redistribute -- the time step of the update must not be zero");
  }

  // The provisional state U = phi - dt*D, with D = (kappa*D)/kappa, on all levels before any level is changed. In the
  // ghost cells it is the state of the coarse cell under the ghost cell.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl  = m_eblg[lvl]->getDBL();
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

    CH_assert(a_divG[lvl]->ghostVect() >= s_numGhost * IntVect::Unit);
    CH_assert(a_phi[lvl]->ghostVect() >= s_numGhost * IntVect::Unit);

    if (!(a_divG[lvl]->ghostVect() >= s_numGhost * IntVect::Unit) ||
        !(a_phi[lvl]->ghostVect() >= s_numGhost * IntVect::Unit)) {
      MayDay::Abort("EBStateRedistribution::redistribute -- the state and the divergence need two ghost cells");
    }

    LevelData<FArrayBox>& state = *m_state[lvl];

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const EBCellFAB&     phiEB    = (*a_phi[lvl])[din];
      const EBCellFAB&     divGEB   = (*a_divG[lvl])[din];
      const FArrayBox&     phi      = phiEB.getFArrayBox();
      const FArrayBox&     divG     = divGEB.getFArrayBox();
      const EBISBox&       geometry = m_eblg[lvl]->getEBISL()[din];
      const ProblemDomain& domain   = m_eblg[lvl]->getDomain();

      FArrayBox& U = state[din];

      U.setVal(0.0);

      auto kernel = [&](const IntVect& iv) -> void {
        const Real kappa = volumeFraction(geometry, domain, iv);

        if (kappa > 0.0) {
          U(iv, 0) = phi(iv, 0) - a_dt * divG(iv, 0) / kappa;
        }
        else if (kappa < 0.0) {
          // Multi-valued cell, which can stand in for another level: the volume-weighted average of its VoFs.
          const Vector<VolIndex> vofs = geometry.getVoFs(iv);

          Real mass   = 0.0;
          Real volume = 0.0;

          for (int i = 0; i < vofs.size(); i++) {
            const Real k = geometry.volFrac(vofs[i]);

            mass += k * phiEB(vofs[i], 0) - a_dt * divGEB(vofs[i], 0);
            volume += k;
          }

          U(iv, 0) = (volume > 0.0) ? mass / volume : 0.0;
        }
      };

      BoxLoops::loop(dbl[din], kernel);
    }

    state.exchange();

    if (lvl > 0) {
      const int refToCoar = m_refinementRatios[lvl - 1];

      LevelData<FArrayBox>& coarseState = *m_ghostMass[lvl];

      m_state[lvl - 1]->copyTo(Interval(0, 0), coarseState, Interval(0, 0), m_coarseCopier[lvl]);

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const FArrayBox&     coarU    = coarseState[din];
        const EBISBox&       geometry = m_eblg[lvl]->getEBISL()[din];
        const ProblemDomain& domain   = m_eblg[lvl]->getDomain();
        const BaseFab<bool>& covered  = (*m_coveredByFiner[lvl])[din];
        const BaseFab<bool>& cfGhost  = (*m_coarseFineGhost[lvl])[din];

        auto kappa = [&](const IntVect& a_iv) -> Real {
          return volumeFraction(geometry, domain, a_iv);
        };
        auto type = [&](const IntVect& a_iv) -> int {
          return cellType(covered, cfGhost, domain, a_iv);
        };

        FArrayBox& U = state[din];

        auto kernel = [&](const IntVect& iv) -> void {
          if (type(iv) == s_ghost && kappa(iv) > 0.0) {
            U(iv, 0) = coarU(coarsen(iv, refToCoar), 0);
          }
        };

        BoxLoops::loop(grow(dbl[din], s_numGhost), kernel);
      }
    }
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl  = m_eblg[lvl]->getDBL();
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

    LevelData<FArrayBox>& average = *m_nbhdAverage[lvl];
    LevelData<FArrayBox>& slopes  = *m_slopes[lvl];

    // The neighbourhood averages of the provisional state. Outside the merged neighbourhoods the average is U itself.
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const FArrayBox&     U        = (*m_state[lvl])[din];
      const EBISBox&       geometry = m_eblg[lvl]->getEBISL()[din];
      const ProblemDomain& domain   = m_eblg[lvl]->getDomain();

      auto kappa = [&](const IntVect& a_iv) -> Real {
        return volumeFraction(geometry, domain, a_iv);
      };

      const FArrayBox&    alp   = (*m_alpha[lvl])[din];
      const BaseFab<int>& nrs   = (*m_numNbhd[lvl])[din];
      const FArrayBox&    vol   = (*m_nbhdVol[lvl])[din];
      const BaseFab<int>& track = (*m_tracker[lvl])[din];

      FArrayBox& avg = average[din];

      avg.setVal(0.0);

      auto kernel = [&](const IntVect& iv) -> void {
        if (kappa(iv) > 0.0) {
          if (track(iv, 0) > 0) {
            Real q = U(iv, 0) * alp(iv, 0) * kappa(iv);

            for (int n = 1; n <= track(iv, 0); n++) {
              const IntVect r = iv + offset(track(iv, n));

              if (kappa(r) > 0.0) {
                q += U(r, 0) * alp(iv, 1) * kappa(r) / nrs(r, 0);
              }
            }

            avg(iv, 0) = q / vol(iv, 0);
          }
          else {
            avg(iv, 0) = U(iv, 0);
          }
        }
      };

      BoxLoops::loop(dbl[din], kernel);
    }

    average.exchange();

    // Ghost cells are in no merged neighbourhood, so their neighbourhood average is U.
    if (lvl > 0) {
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const FArrayBox&     U        = (*m_state[lvl])[din];
        const EBISBox&       geometry = m_eblg[lvl]->getEBISL()[din];
        const ProblemDomain& domain   = m_eblg[lvl]->getDomain();
        const BaseFab<bool>& covered  = (*m_coveredByFiner[lvl])[din];
        const BaseFab<bool>& cfGhost  = (*m_coarseFineGhost[lvl])[din];

        auto type = [&](const IntVect& a_iv) -> int {
          return cellType(covered, cfGhost, domain, a_iv);
        };

        FArrayBox& avg = average[din];

        auto kernel = [&](const IntVect& iv) -> void {
          if (type(iv) == s_ghost) {
            avg(iv, 0) = U(iv, 0);
          }
        };

        BoxLoops::loop(grow(dbl[din], s_numGhost), kernel);
      }
    }

    // Limited slopes of the merged neighbourhoods.
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      FArrayBox& slope = slopes[din];

      slope.setVal(0.0);

      if ((*m_activeBox[lvl])[din]) {
        const FArrayBox&     avg      = average[din];
        const FArrayBox&     centHat  = (*m_nbhdCentroid[lvl])[din];
        const EBISBox&       geometry = m_eblg[lvl]->getEBISL()[din];
        const ProblemDomain& domain   = m_eblg[lvl]->getDomain();
        const BaseFab<bool>& covered  = (*m_coveredByFiner[lvl])[din];
        const BaseFab<bool>& cfGhost  = (*m_coarseFineGhost[lvl])[din];
        const BaseFab<int>&  track    = (*m_tracker[lvl])[din];
        const BaseFab<int>&  conn     = (*m_connected[lvl])[din];

        auto kernel = [&](const IntVect& iv) -> void {
          if (track(iv, 0) > 0) {
            const RealVect s = this->computeSlope(iv, avg, centHat, geometry, covered, cfGhost, domain, conn(iv, 0));

            for (int dir = 0; dir < SpaceDim; dir++) {
              slope(iv, dir) = s[dir];
            }
          }
        };

        BoxLoops::loop(dbl[din], kernel);
      }
    }

    slopes.exchange();

    // The corrections on refinement boundaries. The update of a neighbourhood i is split into transfers between pairs
    // of its cells: the cell at offset b from i gives the cell at offset a
    //
    //    kappa_a*w_a * kappa_b*w_b / V_i * (U_b - U_a + s_i.(x_a - x_b)),
    //
    // where w is the weight of the cell in the neighbourhood and s_i the slope. Transfers from valid cells to covered
    // cells are given to the fine cells under the covered cells, and transfers from valid cells to ghost cells to the
    // coarse cells under the ghost cells. The transfers are gathered on the cells that receive them.
    const int  refToCoar = (lvl > 0) ? m_refinementRatios[lvl - 1] : 1;
    const Real volRatio  = std::pow(static_cast<Real>(refToCoar), SpaceDim);

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      FArrayBox& sync = (*m_syncMass[lvl])[din];

      sync.setVal(0.0);

      if (lvl > 0) {
        (*m_ghostMass[lvl])[din].setVal(0.0);
      }

      if (m_finestLevel > 0 && (*m_activeBox[lvl])[din]) {
        const Box box = dbl[din];

        const EBISBox&       geometry = m_eblg[lvl]->getEBISL()[din];
        const ProblemDomain& domain   = m_eblg[lvl]->getDomain();
        const FArrayBox&     U        = (*m_state[lvl])[din];
        const FArrayBox&     alp      = (*m_alpha[lvl])[din];
        const BaseFab<int>&  nrs      = (*m_numNbhd[lvl])[din];
        const FArrayBox&     vol      = (*m_nbhdVol[lvl])[din];
        const FArrayBox&     slope    = slopes[din];
        const BaseFab<int>&  track    = (*m_tracker[lvl])[din];
        const BaseFab<bool>& covered  = (*m_coveredByFiner[lvl])[din];
        const BaseFab<bool>& cfGhost  = (*m_coarseFineGhost[lvl])[din];

        auto kappa = [&](const IntVect& a_iv) -> Real {
          return volumeFraction(geometry, domain, a_iv);
        };
        auto type = [&](const IntVect& a_iv) -> int {
          return cellType(covered, cfGhost, domain, a_iv);
        };

        // Weight of the cell at offset a_o in the neighbourhood of a_i.
        auto weight = [&](const IntVect& a_i, const IntVect& a_o) -> Real {
          return (a_o == IntVect::Zero) ? alp(a_i, 0) : alp(a_i, 1) / nrs(a_i + a_o, 0);
        };

        // Mass that the cell at offset a_from gives the cell at offset a_to in the neighbourhood of a_i.
        auto transfer = [&](const IntVect& a_i, const IntVect& a_from, const IntVect& a_to) -> Real {
          const IntVect from = a_i + a_from;
          const IntVect to   = a_i + a_to;

          const Real kappaFrom = kappa(from);
          const Real kappaTo   = kappa(to);

          // Cells covered by the EB are in no neighbourhood.
          if (kappaFrom <= 0.0 || kappaTo <= 0.0) {
            return 0.0;
          }

          const RealVect centFrom = cellCentroid(geometry, domain, from);
          const RealVect centTo   = cellCentroid(geometry, domain, to);

          Real dq = U(from, 0) - U(to, 0);

          for (int dir = 0; dir < SpaceDim; dir++) {
            dq += slope(a_i, dir) * (centTo[dir] + a_to[dir] - centFrom[dir] - a_from[dir]);
          }

          return kappaTo * weight(a_i, a_to) * kappaFrom * weight(a_i, a_from) * dq / vol(a_i, 0);
        };

        // Whether the cell at offset index a_idx from a_i is in the neighbourhood of a_i.
        auto isMember = [&](const IntVect& a_i, const int a_idx) -> bool {
          for (int n = 1; n <= track(a_i, 0); n++) {
            if (track(a_i, n) == a_idx) {
              return true;
            }
          }

          return false;
        };

        // Transfers from the valid cells to the covered cells of this box.
        auto coveredKernel = [&](const IntVect& iv) -> void {
          if (type(iv) == s_covered && kappa(iv) > 0.0) {
            Real m = 0.0;

            for (int idx = 0; idx < s_numOffsets; idx++) {
              const IntVect o = offset(idx);
              const IntVect i = iv - o;

              if (o != IntVect::Zero && track(i, 0) > 0 && isMember(i, idx)) {
                m += transfer(i, IntVect::Zero, o);

                for (int n = 1; n <= track(i, 0); n++) {
                  const IntVect om = offset(track(i, n));

                  if (type(i + om) == s_valid) {
                    m += transfer(i, om, o);
                  }
                }
              }
            }

            sync(iv, 0) = m;
          }
        };

        BoxLoops::loop(box, coveredKernel);

        // Transfers from the valid cells to the ghost cells, from the neighbourhoods of the cells in this box, given to
        // the coarse cells under the ghost cells.
        if (lvl > 0) {
          FArrayBox& ghostMass = (*m_ghostMass[lvl])[din];

          auto ghostKernel = [&](const IntVect& iv) -> void {
            if (type(iv) == s_ghost && kappa(iv) > 0.0) {
              Real m = 0.0;

              for (int idx = 0; idx < s_numOffsets; idx++) {
                const IntVect o = offset(idx);
                const IntVect i = iv - o;

                if (o != IntVect::Zero && box.contains(i) && track(i, 0) > 0 && isMember(i, idx)) {
                  m += transfer(i, IntVect::Zero, o);

                  for (int n = 1; n <= track(i, 0); n++) {
                    const IntVect om = offset(track(i, n));

                    if (type(i + om) == s_valid) {
                      m += transfer(i, om, o);
                    }
                  }
                }
              }

              ghostMass(coarsen(iv, refToCoar), 0) += m / volRatio;
            }
          };

          BoxLoops::loop(grow(box, 1), ghostKernel);
        }
      }
    }

    // The new state, and from it the right-hand side R = (phi - phi^{n+1})/dt. Cells outside the merged neighbourhoods
    // keep R = D.
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const FArrayBox&     phi      = (*a_phi[lvl])[din].getFArrayBox();
      const EBISBox&       geometry = m_eblg[lvl]->getEBISL()[din];
      const ProblemDomain& domain   = m_eblg[lvl]->getDomain();

      EBCellFAB& divGEB = (*a_divG[lvl])[din];
      FArrayBox& divG   = divGEB.getFArrayBox();

      const bool active = (*m_activeBox[lvl])[din];

      const FArrayBox&     avg     = average[din];
      const FArrayBox&     slope   = slopes[din];
      const FArrayBox&     centHat = (*m_nbhdCentroid[lvl])[din];
      const FArrayBox&     alp     = (*m_alpha[lvl])[din];
      const BaseFab<int>&  nrs     = (*m_numNbhd[lvl])[din];
      const BaseFab<int>&  track   = (*m_tracker[lvl])[din];
      const BaseFab<bool>& covered = (*m_coveredByFiner[lvl])[din];
      const BaseFab<bool>& cfGhost = (*m_coarseFineGhost[lvl])[din];

      auto kappa = [&](const IntVect& a_iv) -> Real {
        return volumeFraction(geometry, domain, a_iv);
      };
      auto type = [&](const IntVect& a_iv) -> int {
        return cellType(covered, cfGhost, domain, a_iv);
      };

      // Value of the reconstruction of the neighbourhood of a_i at the centroid of the cell at offset a_o from a_i.
      auto reconstruct = [&](const IntVect& a_i, const IntVect& a_o) -> Real {
        const RealVect c = cellCentroid(geometry, domain, a_i + a_o);

        Real value = avg(a_i, 0);

        for (int dir = 0; dir < SpaceDim; dir++) {
          value += slope(a_i, dir) * (c[dir] + a_o[dir] - centHat(a_i, dir));
        }

        return value;
      };

      auto kernel = [&](const IntVect& iv) -> void {
        const Real k = kappa(iv);

        if (k > 0.0) {
          if (active && type(iv) == s_valid && (track(iv, 0) > 0 || nrs(iv, 0) > 1.0)) {
            const Real self = (track(iv, 0) > 0) ? reconstruct(iv, IntVect::Zero) : avg(iv, 0);

            Real sum = alp(iv, 0) * nrs(iv, 0) * self;

            for (int idx = 0; idx < s_numOffsets; idx++) {
              const IntVect o = offset(idx);
              const IntVect i = iv - o;

              if (o != IntVect::Zero) {
                for (int n = 1; n <= track(i, 0); n++) {
                  if (track(i, n) == idx) {
                    sum += alp(i, 1) * reconstruct(i, o);
                  }
                }
              }
            }

            const Real phiNew = sum / nrs(iv, 0);

            divG(iv, 0) = (phi(iv, 0) - phiNew) / a_dt;
          }
          else if (k < 1.0) {
            divG(iv, 0) /= k;
          }
        }
        else if (k < 0.0) {
          // Multi-valued cell, which is in no merged neighbourhood: R = D in each of its VoFs.
          const Vector<VolIndex> vofs = geometry.getVoFs(iv);

          for (int i = 0; i < vofs.size(); i++) {
            divGEB(vofs[i], 0) /= geometry.volFrac(vofs[i]);
          }
        }
      };

      BoxLoops::loop(dbl[din], kernel);
    }
  }

  if (m_finestLevel > 0) {
    this->synchronize(a_divG, a_dt);
  }
}

void
EBStateRedistribution::synchronize(EBAMRCellData& a_divG, const Real a_dt) noexcept
{
  CH_TIME("EBStateRedistribution::synchronize");

  // The corrections from the ghost cells go to the coarse cells under them.
  for (int lvl = 1; lvl <= m_finestLevel; lvl++) {
    m_ghostMass[lvl]->copyTo(Interval(0, 0),
                             *m_syncMass[lvl - 1],
                             Interval(0, 0),
                             m_ghostCopier[lvl],
                             LDaddOp<FArrayBox>());
  }

  // The corrections of the covered cells go to the fine cells under them, as the same change of the state in each fine
  // cell. A fine cell that is itself covered passes its part on to the next level. This reuses m_ghostMass, which holds
  // the coarse data on the coarsening of the fine grids.
  for (int lvl = 0; lvl < m_finestLevel; lvl++) {
    const int  refToFine = m_refinementRatios[lvl];
    const Real volRatio  = std::pow(static_cast<Real>(refToFine), SpaceDim);

    LevelData<FArrayBox>& coarseMass = *m_ghostMass[lvl + 1];

    m_syncMass[lvl]->copyTo(Interval(0, 0), coarseMass, Interval(0, 0), m_coarseCopier[lvl + 1]);

    const DisjointBoxLayout& dblFine  = m_eblg[lvl + 1]->getDBL();
    const DataIterator&      ditFine  = dblFine.dataIterator();
    const int                nboxFine = ditFine.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nboxFine; mybox++) {
      const DataIndex& din = ditFine[mybox];

      const FArrayBox&     mass     = coarseMass[din];
      const EBISBox&       geometry = m_eblg[lvl + 1]->getEBISL()[din];
      const ProblemDomain& domain   = m_eblg[lvl + 1]->getDomain();
      const BaseFab<bool>& covered  = (*m_coveredByFiner[lvl + 1])[din];
      const BaseFab<bool>& cfGhost  = (*m_coarseFineGhost[lvl + 1])[din];

      auto kappa = [&](const IntVect& a_iv) -> Real {
        return volumeFraction(geometry, domain, a_iv);
      };
      auto type = [&](const IntVect& a_iv) -> int {
        return cellType(covered, cfGhost, domain, a_iv);
      };

      EBCellFAB& divGEB = (*a_divG[lvl + 1])[din];
      FArrayBox& divG   = divGEB.getFArrayBox();
      FArrayBox& sync   = (*m_syncMass[lvl + 1])[din];

      // Volume fraction of a fine cell, summed over the VoFs of a multi-valued cell.
      auto cellVolume = [&](const IntVect& a_iv) -> Real {
        if (kappa(a_iv) >= 0.0) {
          return kappa(a_iv);
        }

        const Vector<VolIndex> vofs = geometry.getVoFs(a_iv);

        Real volume = 0.0;

        for (int i = 0; i < vofs.size(); i++) {
          volume += geometry.volFrac(vofs[i]);
        }

        return volume;
      };

      auto kernel = [&](const IntVect& ivCoar) -> void {
        const Real m = mass(ivCoar, 0);

        if (m != 0.0) {
          const Box fineCells = refine(Box(ivCoar, ivCoar), refToFine);

          Real sumKappa = 0.0;

          for (BoxIterator bit(fineCells); bit.ok(); ++bit) {
            sumKappa += cellVolume(bit());
          }

          // The mass in units of the fine cell volume, spread as the same change of the state.
          const Real dphi = m * volRatio / sumKappa;

          for (BoxIterator bit(fineCells); bit.ok(); ++bit) {
            const IntVect& iv = bit();

            if (type(iv) == s_covered) {
              sync(iv, 0) += cellVolume(iv) * dphi;
            }
            else if (kappa(iv) > 0.0) {
              divG(iv, 0) -= dphi / a_dt;
            }
            else if (kappa(iv) < 0.0) {
              const Vector<VolIndex> vofs = geometry.getVoFs(iv);

              for (int i = 0; i < vofs.size(); i++) {
                divGEB(vofs[i], 0) -= dphi / a_dt;
              }
            }
          }
        }
      };

      BoxLoops::loop(coarsen(dblFine[din], refToFine), kernel);
    }
  }

  // A cell with volume fraction kappa takes the fraction kappa of its correction M, which changes the state by M/kappa
  // * kappa = M. The rest is redistributed to its neighbours.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl  = m_eblg[lvl]->getDBL();
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const FArrayBox&     sync     = (*m_syncMass[lvl])[din];
      const EBISBox&       geometry = m_eblg[lvl]->getEBISL()[din];
      const ProblemDomain& domain   = m_eblg[lvl]->getDomain();
      const BaseFab<bool>& covered  = (*m_coveredByFiner[lvl])[din];
      const BaseFab<bool>& cfGhost  = (*m_coarseFineGhost[lvl])[din];

      auto kappa = [&](const IntVect& a_iv) -> Real {
        return volumeFraction(geometry, domain, a_iv);
      };
      auto type = [&](const IntVect& a_iv) -> int {
        return cellType(covered, cfGhost, domain, a_iv);
      };

      EBCellFAB&       divGEB   = (*a_divG[lvl])[din];
      FArrayBox&       divG     = divGEB.getFArrayBox();
      BaseIVFAB<Real>& massDiff = (*m_massDiff[lvl])[din];

      massDiff.setVal(0.0);

      auto kernel = [&](const IntVect& iv) -> void {
        const Real m = sync(iv, 0);

        if (type(iv) == s_valid && m != 0.0) {
          const Real k = kappa(iv);

          // A multi-valued cell takes all of the correction, as the same change of the state in each VoF.
          if (k < 0.0) {
            const Vector<VolIndex> vofs = geometry.getVoFs(iv);

            Real volume = 0.0;

            for (int i = 0; i < vofs.size(); i++) {
              volume += geometry.volFrac(vofs[i]);
            }

            for (int i = 0; i < vofs.size(); i++) {
              divGEB(vofs[i], 0) -= m / (volume * a_dt);
            }

            return;
          }

          divG(iv, 0) -= m / a_dt;

          if (k < 1.0) {
            massDiff(VolIndex(iv, 0), 0) = -(1.0 - k) * m / a_dt;
          }
        }
      };

      BoxLoops::loop(dbl[din], kernel);
    }
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const Real     scale     = 1.0;
    const Interval variables = Interval(0, 0);

    if (lvl > 0) {
      m_fluxRedistribution[lvl]->redistributeCoar(*a_divG[lvl - 1], *m_massDiff[lvl], scale, variables);
    }

    m_fluxRedistribution[lvl]->redistributeLevel(*a_divG[lvl], *m_massDiff[lvl], scale, variables);

    if (lvl < m_finestLevel) {
      m_fluxRedistribution[lvl]->redistributeFine(*a_divG[lvl + 1], *m_massDiff[lvl], scale, variables);
    }
  }
}

#include <CD_NamespaceFooter.H>
