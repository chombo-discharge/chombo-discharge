/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ParticleMergeStepper.cpp
 * @brief  Implementation of CD_ParticleMergeStepper.H
 * @author Robert Marskar
 */

// Std includes
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ParticleMergeStepper.H>
#include <CD_ParallelOps.H>
#include <CD_Random.H>
#include <CD_NamespaceHeader.H>

namespace {

/**
 * @brief Reduction of the merged particle population, accumulated over the valid cells of one level.
 * @details Everything here is a histogram or a scalar reduction -- there is no per-cell storage. The
 * per-cell quantities are read straight out of the container's cell-sorted CSR mapping.
 */
struct Stats
{
  static constexpr int s_numSubBins = 50;  ///< Bins for the sub-cell coordinate histograms.
  static constexpr int s_maxCount   = 96;  ///< Largest per-cell particle count that is binned.
  static constexpr int s_maxWeight  = 600; ///< Largest particle weight that is binned.
  static constexpr int s_jointBins  = 6;   ///< Bins per direction in the joint sub-cell histogram.
  static constexpr int s_numModes   = 8;   ///< Sub-cell Fourier modes k = 1..s_numModes.

  /**
   * @brief Construct the histograms, sized from the patch dimension.
   * @param[in] a_blockSize Patch size in cells (for the patch-relative histograms).
   */
  Stats(const int a_blockSize) : m_blockSize(a_blockSize)
  {
    for (int dir = 0; dir < SpaceDim; dir++) {
      m_subCount[dir].assign(s_numSubBins, 0.0);
      m_subWeight[dir].assign(s_numSubBins, 0.0);
      m_patchWeight[dir].assign(a_blockSize, 0.0);
      m_patchCount[dir].assign(a_blockSize, 0.0);
      m_blockW[dir].assign(a_blockSize, 0.0);
      m_blockW2[dir].assign(a_blockSize, 0.0);
      m_blockN[dir].assign(a_blockSize, 0.0);
      m_cosK[dir].assign(s_numModes, 0.0);
      m_sinK[dir].assign(s_numModes, 0.0);
    }
    m_cellCount.assign(s_maxCount + 1, 0.0);
    m_weightHist.assign(s_maxWeight + 1, 0.0);
    std::size_t jointSize = 1;
    for (int dir = 0; dir < SpaceDim; dir++) {
      jointSize *= s_jointBins;
    }
    m_joint.assign(jointSize, 0.0);
  }

  int                 m_blockSize;             ///< Patch size in cells.
  std::vector<double> m_subCount[SpaceDim];    ///< Sub-cell coordinate histogram, per particle.
  std::vector<double> m_subWeight[SpaceDim];   ///< Sub-cell coordinate histogram, weighted.
  std::vector<double> m_patchWeight[SpaceDim]; ///< Weight vs patch-relative cell index.
  std::vector<double> m_patchCount[SpaceDim];  ///< Particle count vs patch-relative cell index.
  std::vector<double> m_blockW[SpaceDim];      ///< Sum of per-cell weight vs block phase.
  std::vector<double> m_blockW2[SpaceDim];     ///< Sum of per-cell weight^2 vs block phase.
  std::vector<double> m_blockN[SpaceDim];      ///< Cell count vs block phase.
  std::vector<double> m_cosK[SpaceDim];        ///< Sum cos(2*pi*k*u) over particles, k = 1..s_numModes.
  std::vector<double> m_sinK[SpaceDim];        ///< Sum sin(2*pi*k*u) over particles, k = 1..s_numModes.
  std::vector<double> m_cellCount;             ///< Histogram of particles per cell.
  std::vector<double> m_weightHist;            ///< Histogram of per-particle weight (0..s_maxWeight).
  std::vector<double> m_joint;                 ///< Joint sub-cell position histogram (s_jointBins^SpaceDim).
  double              m_sumU2         = 0.0;   ///< Sum of (u-0.5)^2 over particles and directions.
  double              m_sumNeff       = 0.0;   ///< Sum over cells of (sum w)^2 / (sum w^2).
  double              m_sumTopFrac    = 0.0;   ///< Sum over cells of (heaviest weight) / (cell weight).
  double              m_maxWeight     = 0.0;   ///< Largest single particle weight seen.
  double              m_minWeight     = std::numeric_limits<double>::max(); ///< Smallest particle weight seen.
  double              m_totalWeight   = 0.0;                                ///< Total weight on this level.
  double              m_totalParticle = 0.0;                                ///< Total particles on this level.
  double              m_numCells      = 0.0;                                ///< Number of valid cells visited.
  double              m_sumW          = 0.0; ///< Sum of per-cell weight (== m_totalWeight).
  double              m_sumW2         = 0.0; ///< Sum of per-cell weight squared.
  double              m_absDev        = 0.0; ///< Sum |cellWeight - nominal| over valid cells.

  /**
   * @brief Reduce every accumulator across MPI ranks.
   */
  void
  reduce()
  {
    for (int dir = 0; dir < SpaceDim; dir++) {
      for (auto& v : m_subCount[dir]) {
        v = ParallelOps::sum(v);
      }
      for (auto& v : m_subWeight[dir]) {
        v = ParallelOps::sum(v);
      }
      for (auto& v : m_patchWeight[dir]) {
        v = ParallelOps::sum(v);
      }
      for (auto& v : m_patchCount[dir]) {
        v = ParallelOps::sum(v);
      }
      for (auto& v : m_blockW[dir]) {
        v = ParallelOps::sum(v);
      }
      for (auto& v : m_blockW2[dir]) {
        v = ParallelOps::sum(v);
      }
      for (auto& v : m_blockN[dir]) {
        v = ParallelOps::sum(v);
      }
      for (auto& v : m_cosK[dir]) {
        v = ParallelOps::sum(v);
      }
      for (auto& v : m_sinK[dir]) {
        v = ParallelOps::sum(v);
      }
    }
    for (auto& v : m_cellCount) {
      v = ParallelOps::sum(v);
    }
    for (auto& v : m_weightHist) {
      v = ParallelOps::sum(v);
    }
    for (auto& v : m_joint) {
      v = ParallelOps::sum(v);
    }
    m_sumU2      = ParallelOps::sum(m_sumU2);
    m_sumNeff    = ParallelOps::sum(m_sumNeff);
    m_sumTopFrac = ParallelOps::sum(m_sumTopFrac);
    m_maxWeight  = ParallelOps::max(m_maxWeight);
    m_minWeight  = ParallelOps::min(m_minWeight);

    m_totalWeight   = ParallelOps::sum(m_totalWeight);
    m_totalParticle = ParallelOps::sum(m_totalParticle);
    m_numCells      = ParallelOps::sum(m_numCells);
    m_sumW          = ParallelOps::sum(m_sumW);
    m_sumW2         = ParallelOps::sum(m_sumW2);
    m_absDev        = ParallelOps::sum(m_absDev);
  }
};

/**
 * @brief Walk the cell-sorted container and fill the reduction.
 * @param[out]    a_stats     Reduction to fill.
 * @param[in,out] a_particles Particle container (cell-sorted in place).
 * @param[in]     a_amr       AMR mesh.
 */
void
gatherStats(Stats&                          a_stats,
            ParticleContainer<ItoParticle>& a_particles,
            const RefCountedPtr<AmrMesh>&   a_amr,
            const std::string&              a_realm,
            const int                       a_level)
{
  a_particles.organizeParticlesByCell();

  const RealVect probLo = a_amr->getProbLo();

  {
    const int                lvl = a_level;
    const DisjointBoxLayout& dbl = a_amr->getGrids(a_realm)[lvl];
    const Real               dx  = a_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const Box                       box        = dbl[dit()];
      const ParticleSoA<ItoParticle>& leaf       = a_particles[lvl][dit()];
      const BaseFab<bool>&            validCells = (*a_amr->getValidCells(a_realm)[lvl])[dit()];

      for (std::size_t c = 0; c < leaf.numCells(); c++) {
        // Recover the cell index from the Fortran-ordered CSR cell key.
        const IntVect boxSize = box.size();
        std::size_t   rem     = c;
        IntVect       iv      = box.smallEnd();
        for (int dir = 0; dir < SpaceDim; dir++) {
          iv[dir] += static_cast<int>(rem % static_cast<std::size_t>(boxSize[dir]));
          rem /= static_cast<std::size_t>(boxSize[dir]);
        }

        // Skip cells the domain margin excludes, and cells covered by a finer level. A covered coarse
        // cell holds no particles -- they live on the finest level that covers them -- so counting it
        // would report an empty cell that is not empty, and would make the per-level weight look as
        // though the merge had lost it.
        // Skip cells covered by a finer level. A covered coarse cell holds no particles -- they live on
        // the finest level that covers them -- so counting it would report an empty cell that is not
        // empty, and would make the per-level weight look as though the merge had lost it.
        if (!validCells(iv, 0)) {
          continue;
        }

        const std::size_t lo = leaf.cellStart(c);
        const std::size_t hi = lo + leaf.particlesInCell(c);

        double cellWeight        = 0.0;
        double cellWeightSquared = 0.0;
        double cellTopWeight     = 0.0;

        for (std::size_t i = lo; i < hi; i++) {
          const RealVect x = leaf.position(i);
          const Real     w = leaf.weight(i);

          int jbin = 0;

          cellWeight += w;
          cellWeightSquared += w * w;
          cellTopWeight = std::max(cellTopWeight, static_cast<double>(w));

          a_stats.m_maxWeight = std::max(a_stats.m_maxWeight, static_cast<double>(w));
          a_stats.m_minWeight = std::min(a_stats.m_minWeight, static_cast<double>(w));

          a_stats.m_weightHist[std::max(0, std::min(Stats::s_maxWeight, static_cast<int>(std::lround(w))))] += 1.0;

          for (int dir = 0; dir < SpaceDim; dir++) {
            const double u   = (x[dir] - probLo[dir]) / dx - static_cast<double>(iv[dir]);
            int          bin = static_cast<int>(u * Stats::s_numSubBins);
            bin              = std::max(0, std::min(Stats::s_numSubBins - 1, bin));

            a_stats.m_subCount[dir][bin] += 1.0;
            a_stats.m_subWeight[dir][bin] += w;

            const int pb = ((iv[dir] % a_stats.m_blockSize) + a_stats.m_blockSize) % a_stats.m_blockSize;
            a_stats.m_patchWeight[dir][pb] += w;
            a_stats.m_patchCount[dir][pb] += 1.0;

            a_stats.m_sumU2 += (u - 0.5) * (u - 0.5);

            // Sub-cell Fourier amplitudes. |Z_k| = |mean exp(2*pi*i*k*u)| is 0 for a uniform sub-cell
            // distribution (floor ~ 1/sqrt(N)) and O(1) if the particles pile onto a lattice of period
            // 1/k. Two-term recurrence rather than s_numModes trig calls per particle per direction.
            {
              const double th = 2.0 * M_PI * u;
              const double c1 = std::cos(th);
              const double s1 = std::sin(th);

              double cPrev = 1.0;
              double sPrev = 0.0;
              double cCur  = c1;
              double sCur  = s1;

              for (int k = 0; k < Stats::s_numModes; k++) {
                a_stats.m_cosK[dir][k] += cCur;
                a_stats.m_sinK[dir][k] += sCur;

                const double cNext = 2.0 * c1 * cCur - cPrev;
                const double sNext = 2.0 * c1 * sCur - sPrev;

                cPrev = cCur;
                sPrev = sCur;
                cCur  = cNext;
                sCur  = sNext;
              }
            }

            int jb = static_cast<int>(u * Stats::s_jointBins);
            jb     = std::max(0, std::min(Stats::s_jointBins - 1, jb));
            jbin   = jbin * Stats::s_jointBins + jb;
          }

          a_stats.m_joint[jbin] += 1.0;
        }

        const int nInCell = static_cast<int>(hi - lo);

        for (int dir = 0; dir < SpaceDim; dir++) {
          const int pb = ((iv[dir] % a_stats.m_blockSize) + a_stats.m_blockSize) % a_stats.m_blockSize;

          a_stats.m_blockW[dir][pb] += cellWeight;
          a_stats.m_blockW2[dir][pb] += cellWeight * cellWeight;
          a_stats.m_blockN[dir][pb] += 1.0;
        }

        // Effective particles per cell: (sum w)^2 / (sum w^2) equals the particle count exactly when the
        // cell's mass is spread evenly over its particles, and falls below it as the weights spread out.
        if (cellWeightSquared > 0.0) {
          a_stats.m_sumNeff += cellWeight * cellWeight / cellWeightSquared;
          a_stats.m_sumTopFrac += cellTopWeight / cellWeight;
        }

        a_stats.m_cellCount[std::min(nInCell, Stats::s_maxCount)] += 1.0;
        a_stats.m_totalWeight += cellWeight;
        a_stats.m_totalParticle += static_cast<double>(nInCell);
        a_stats.m_numCells += 1.0;
        a_stats.m_sumW += cellWeight;
        a_stats.m_sumW2 += cellWeight * cellWeight;
      }
    }
  }
}

/**
 * @brief Seed a_ppc unit-weight particles per valid cell, uniformly distributed.
 * @details One Poisson draw per patch with positions uniform over the patch, rather than a fixed count
 * placed inside each cell. Both give the same mean density, but the unstratified draw carries no imprint
 * of the cell structure -- a fixed count per cell makes a count-median split land exactly on cell faces,
 * which is a degenerate case rather than a representative one. Positions landing in cells a finer level
 * covers are rejected: those particles would be discarded on the next remap, so keeping them would make
 * the seeded population differ from the one the merge actually sees.
 * @param[in,out] a_particles Container to seed into.
 * @param[in]     a_amr       AMR mesh.
 * @param[in]     a_realm     Realm the container lives on.
 * @param[in]     a_ppc       Particles per cell.
 */
void
fillUniform(ParticleContainer<ItoParticle>& a_particles,
            const RefCountedPtr<AmrMesh>&   a_amr,
            const std::string&              a_realm,
            const int                       a_ppc)
{
  const RealVect probLo = a_amr->getProbLo();

  for (int lvl = 0; lvl <= a_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_amr->getGrids(a_realm)[lvl];
    const Real               dx  = a_amr->getDx()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&          din        = dit[mybox];
      const Box                 box        = dbl[din];
      ParticleSoA<ItoParticle>& leaf       = a_particles[lvl][din];
      const BaseFab<bool>&      validCells = (*a_amr->getValidCells(a_realm)[lvl])[din];

      const RealVect patchLo = probLo + dx * RealVect(box.smallEnd());
      const RealVect patchHi = probLo + dx * RealVect(box.bigEnd() + IntVect::Unit);

      const long long numDraw = Random::getPoisson<long long>(static_cast<Real>(box.numPts()) * a_ppc);

      leaf.reserve(leaf.size() + static_cast<std::size_t>(numDraw));

      for (long long p = 0; p < numDraw; p++) {
        RealVect x = RealVect::Zero;

        for (int dir = 0; dir < SpaceDim; dir++) {
          x[dir] = patchLo[dir] + (patchHi[dir] - patchLo[dir]) * Random::getUniformReal01();
        }

        IntVect iv;
        for (int dir = 0; dir < SpaceDim; dir++) {
          iv[dir] = static_cast<int>(std::floor((x[dir] - probLo[dir]) / dx));
        }

        if (validCells(iv, 0)) {
          leaf.append(x, 1.0);
        }
      }
    }
  }
}

/**
 * @brief Print one histogram as "bin value ratio-to-flat".
 * @param[in] a_name Histogram label.
 * @param[in] a_hist Histogram contents.
 */
void
report(const std::string& a_name, const std::vector<double>& a_hist)
{
  double sum = 0.0;
  for (const double v : a_hist) {
    sum += v;
  }

  if (sum <= 0.0) {
    return;
  }

  const double flat = sum / static_cast<double>(a_hist.size());

  double maxDev = 0.0;
  double chi2   = 0.0;
  for (const double v : a_hist) {
    maxDev = std::max(maxDev, std::abs(v / flat - 1.0));
    chi2 += (v - flat) * (v - flat) / flat;
  }

  pout() << a_name << ": bins=" << a_hist.size() << " max|dev|=" << 100.0 * maxDev << "%"
         << " chi2/dof=" << chi2 / static_cast<double>(a_hist.size() - 1) << endl;

  pout() << "   ";
  for (const double v : a_hist) {
    pout() << std::fixed << std::setprecision(4) << v / flat << " ";
  }
  pout() << endl;
}

/**
 * @brief Print the before/after report for one AMR level, laid out to be read.
 * @param[in] a_pre       Statistics gathered before the merge.
 * @param[in] a_post      Statistics gathered after the merge.
 * @param[in] a_level     AMR level the statistics belong to.
 * @param[in] a_dx        Grid spacing on that level.
 * @param[in] a_targetPPC Merge target, for the population and N_eff comparisons.
 * @param[in] a_blockSize Patch size in cells, for the block-phase comparison.
 */
void
reportLevel(const Stats& a_pre,
            const Stats& a_post,
            const int    a_level,
            const Real   a_dx,
            const int    a_targetPPC,
            const int    a_blockSize)
{
  if (a_post.m_numCells <= 0.0 || a_pre.m_numCells <= 0.0 || a_post.m_totalParticle <= 0.0) {
    pout() << "   Level " << a_level << ": no valid cells on this rank's share of the level" << endl;

    return;
  }

  const double meanPre   = a_pre.m_sumW / a_pre.m_numCells;
  const double meanPost  = a_post.m_sumW / a_post.m_numCells;
  const double sigmaPre  = std::sqrt(std::max(0.0, a_pre.m_sumW2 / a_pre.m_numCells - meanPre * meanPre));
  const double sigmaPost = std::sqrt(std::max(0.0, a_post.m_sumW2 / a_post.m_numCells - meanPost * meanPost));

  int minN = Stats::s_maxCount;
  int maxN = 0;
  for (std::size_t i = 0; i < a_post.m_cellCount.size(); i++) {
    if (a_post.m_cellCount[i] > 0.0) {
      minN = std::min(minN, static_cast<int>(i));
      maxN = std::max(maxN, static_cast<int>(i));
    }
  }

  // Strongest sub-cell Fourier mode, over every direction. A centroid placement drives the particles
  // onto a lattice of leaf centres, which shows up here as a large amplitude at k = ppc^(1/D).
  double bestZ = 0.0;
  int    bestK = 0;
  int    bestD = 0;
  for (int dir = 0; dir < SpaceDim; dir++) {
    for (int k = 0; k < Stats::s_numModes; k++) {
      const double c = a_post.m_cosK[dir][k] / a_post.m_totalParticle;
      const double s = a_post.m_sinK[dir][k] / a_post.m_totalParticle;
      const double z = std::sqrt(c * c + s * s);

      if (z > bestZ) {
        bestZ = z;
        bestK = k + 1;
        bestD = dir;
      }
    }
  }

  double blockLo = 1.E30;
  double blockHi = 0.0;
  for (int b = 0; b < a_blockSize; b++) {
    const double v = (a_post.m_blockW[0][b] / a_post.m_blockN[0][b]) / meanPost;

    blockLo = std::min(blockLo, v);
    blockHi = std::max(blockHi, v);
  }

  const double e2      = a_post.m_sumU2 / (SpaceDim * a_post.m_totalParticle);
  const double uniform = 1.0 / 12.0;
  const double neff    = a_post.m_sumNeff / a_post.m_numCells;
  const double floorZ  = 1.0 / std::sqrt(a_post.m_totalParticle);

  pout() << "   Level " << a_level << "   dx = " << a_dx << "   cells = " << static_cast<long long>(a_post.m_numCells)
         << endl;
  pout() << "     population        : " << std::fixed << std::setprecision(2)
         << a_post.m_totalParticle / a_post.m_numCells << " particles/cell (target " << a_targetPPC << ", range "
         << minN << " .. " << maxN << ")" << endl;
  pout() << "     weight            : " << std::scientific << std::setprecision(3) << a_pre.m_totalWeight << " -> "
         << a_post.m_totalWeight << "   (change " << (a_post.m_totalWeight - a_pre.m_totalWeight) / a_pre.m_totalWeight
         << ")" << endl;
  pout() << "     cell density      : sigma " << std::fixed << std::setprecision(3) << sigmaPre << " -> " << sigmaPost
         << "   (ratio " << sigmaPost / std::max(sigmaPre, std::numeric_limits<double>::min())
         << ", spans transport AND merge)" << endl;
  pout() << "     weight spread     : N_eff " << std::setprecision(2) << neff << " of " << a_targetPPC
         << "   (heaviest particle holds " << std::setprecision(1) << 100.0 * a_post.m_sumTopFrac / a_post.m_numCells
         << "% of its cell, weights " << std::setprecision(0) << a_post.m_minWeight << " .. " << a_post.m_maxWeight
         << ")" << endl;
  pout() << "     sub-cell spread   : E[(u-0.5)^2] = " << std::setprecision(6) << e2 << "   (uniform " << uniform
         << ", " << std::showpos << std::setprecision(1) << 100.0 * (e2 / uniform - 1.0) << "%)" << std::noshowpos
         << endl;
  pout() << "     sub-cell lattice  : max |Z_k| = " << std::setprecision(4) << bestZ << " at k = " << bestK << " (dir "
         << bestD << ")   (shot floor " << floorZ << ")" << endl;
  pout() << "     block frequency   : " << std::setprecision(2) << 100.0 * (blockHi - blockLo)
         << "% peak-to-peak across the " << a_blockSize << "-cell patch" << endl;
  pout() << std::defaultfloat;
}

} // namespace
ParticleMergeStepper::ParticleMergeStepper(RefCountedPtr<ItoSolver>& a_solver)
  : Physics::BrownianWalker::BrownianWalkerStepper(a_solver)
{
  CH_TIME("ParticleMergeStepper::ParticleMergeStepper");

  ParmParse pp("ParticleMergeStepper");

  pp.get("new_particles_per_cell", m_newParticlesPerCell);
  pp.get("print_stats", m_printStats);
}

ParticleMergeStepper::~ParticleMergeStepper() = default;

void
ParticleMergeStepper::seedParticles() noexcept
{
  CH_TIME("ParticleMergeStepper::seedParticles");

  ParticleContainer<ItoParticle>& particles = m_solver->getParticles(ItoSolver::WhichContainer::Bulk);

  particles.organizeParticlesByPatch();

  fillUniform(particles, m_amr, m_realm, m_newParticlesPerCell);
}

void
ParticleMergeStepper::initialData()
{
  CH_TIME("ParticleMergeStepper::initialData");
  if (m_verbosity > 5) {
    pout() << "ParticleMergeStepper::initialData" << endl;
  }

  this->seedParticles();

  m_solver->remap();
  m_solver->depositParticles();

  m_solver->setParticleDiffusion(m_diffCo);
  m_solver->setParticleMobility(m_mobility);
  m_solver->interpolateVelocities();
}

Real
ParticleMergeStepper::advance(const Real a_dt)
{
  CH_TIME("ParticleMergeStepper::advance");
  if (m_verbosity > 5) {
    pout() << "ParticleMergeStepper::advance" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_solver->getParticles(ItoSolver::WhichContainer::Bulk);

  this->seedParticles();

  std::vector<Stats> pre;

  if (m_printStats) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      pre.emplace_back(m_amr->getMaxBoxSize());
      gatherStats(pre.back(), particles, m_amr, m_realm, lvl);
      pre.back().reduce();
    }
  }

  particles.organizeParticlesByPatch();

  // Transport and merge. The merge down to ItoSolver.particles_per_cell happens inside here.
  const Real dt = Physics::BrownianWalker::BrownianWalkerStepper::advance(a_dt);

  if (m_printStats) {
    std::vector<Stats> post;

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      post.emplace_back(m_amr->getMaxBoxSize());
      gatherStats(post.back(), particles, m_amr, m_realm, lvl);
      post.back().reduce();
    }

    pout() << endl
           << "-------------------------------------------------------------------------------" << endl
           << " ParticleMergeStepper -- merge statistics, step " << m_timeStep + 1 << endl
           << "-------------------------------------------------------------------------------" << endl;

    for (int lvl = 0; lvl <= m_amr->getFinestLevel() && lvl < static_cast<int>(post.size()); lvl++) {
      reportLevel(pre[lvl],
                  post[lvl],
                  lvl,
                  m_amr->getDx()[lvl],
                  m_solver->getParticlesPerCell()[0],
                  m_amr->getMaxBoxSize());
    }

    pout() << "-------------------------------------------------------------------------------" << endl << endl;
  }

  return dt;
}

#include <CD_NamespaceFooter.H>
