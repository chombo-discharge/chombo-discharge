/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Example which measures what ItoSolver super-particle merging does to a uniform particle population
 * @author Robert Marskar
 */

// Std includes
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_AmrMesh.H>
#include <CD_ItoSolver.H>
#include <CD_ItoSpecies.H>
#include <CD_ComputationalGeometry.H>
#include <CD_Random.H>
#include <CD_ParallelOps.H>
#include <CD_DischargeIO.H>

using namespace ChomboDischarge;

namespace {

/**
 * @brief Reduction of the merged particle population, accumulated over the interior cells only.
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
  double              m_totalWeight   = 0.0;                                ///< Total weight in the interior.
  double              m_totalParticle = 0.0;                                ///< Total particles in the interior.
  double              m_numCells      = 0.0;                                ///< Number of interior cells visited.
  double              m_sumW          = 0.0; ///< Sum of per-cell weight (== m_totalWeight).
  double              m_sumW2         = 0.0; ///< Sum of per-cell weight squared.
  double              m_absDev        = 0.0; ///< Sum |cellWeight - nominal| over interior cells.

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
 * @param[in]     a_margin    Cells excluded next to each domain face.
 */
void
gatherStats(Stats&                          a_stats,
            ParticleContainer<ItoParticle>& a_particles,
            const RefCountedPtr<AmrMesh>&   a_amr,
            const int                       a_margin)
{
  a_particles.organizeParticlesByCell();

  const RealVect probLo = a_amr->getProbLo();

  for (int lvl = 0; lvl <= a_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl    = a_amr->getGrids(Realm::Primal)[lvl];
    const Real               dx     = a_amr->getDx()[lvl];
    const Box                domBox = a_amr->getDomains()[lvl].domainBox();

    // The analysis region: the domain shrunk by a_margin cells on every face. Everything outside it is
    // discarded, so the cells that ARE analysed all sit in a fully populated neighbourhood and never see
    // the domain edge. This is what keeps the edge from contaminating the aliasing signal.
    Box interior = domBox;
    interior.grow(-a_margin);

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const Box                       box  = dbl[dit()];
      const ParticleSoA<ItoParticle>& leaf = a_particles[lvl][dit()];

      for (std::size_t c = 0; c < leaf.numCells(); c++) {
        // Recover the cell index from the Fortran-ordered CSR cell key.
        const IntVect boxSize = box.size();
        std::size_t   rem     = c;
        IntVect       iv      = box.smallEnd();
        for (int dir = 0; dir < SpaceDim; dir++) {
          iv[dir] += static_cast<int>(rem % static_cast<std::size_t>(boxSize[dir]));
          rem /= static_cast<std::size_t>(boxSize[dir]);
        }

        if (!interior.contains(iv)) {
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
 * @brief Fill the container with a_ppc unit-weight particles per cell, uniform inside each cell.
 * @param[in,out] a_particles Container to fill.
 * @param[in]     a_amr       AMR mesh.
 * @param[in]     a_ppc       Particles per cell.
 */
void
fillUniform(ParticleContainer<ItoParticle>& a_particles,
            const RefCountedPtr<AmrMesh>&   a_amr,
            const int                       a_ppc,
            const bool                      a_stratified)
{
  const RealVect probLo = a_amr->getProbLo();

  for (int lvl = 0; lvl <= a_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_amr->getGrids(Realm::Primal)[lvl];
    const Real               dx  = a_amr->getDx()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&          din  = dit[mybox];
      const Box                 box  = dbl[din];
      ParticleSoA<ItoParticle>& leaf = a_particles[lvl][din];

      leaf.reserve(leaf.size() + static_cast<std::size_t>(box.numPts()) * a_ppc);

      // Unstratified: one Poisson draw for the whole patch, positions uniform over the patch. This is
      // the same mean density, but the per-cell counts are no longer identical and the sample carries
      // no imprint of the cell structure -- which is exactly what a count-median kd split needs in
      // order NOT to land on cell faces.
      if (!a_stratified) {
        const RealVect patchLo = probLo + dx * RealVect(box.smallEnd());
        const RealVect patchHi = probLo + dx * RealVect(box.bigEnd() + IntVect::Unit);

        const long long numDraw = Random::getPoisson<long long>(static_cast<Real>(box.numPts()) * a_ppc);

        for (long long p = 0; p < numDraw; p++) {
          RealVect x = RealVect::Zero;
          for (int dir = 0; dir < SpaceDim; dir++) {
            x[dir] = patchLo[dir] + (patchHi[dir] - patchLo[dir]) * Random::getUniformReal01();
          }

          leaf.append(x, 1.0);
        }

        continue;
      }

      for (BoxIterator bit(box); bit.ok(); ++bit) {
        const IntVect  iv = bit();
        const RealVect lo = probLo + dx * RealVect(iv);

        for (int p = 0; p < a_ppc; p++) {
          RealVect x = RealVect::Zero;
          for (int dir = 0; dir < SpaceDim; dir++) {
            x[dir] = lo[dir] + dx * Random::getUniformReal01();
          }

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

} // namespace

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  {
    ParmParse pp("ParticleMerge");

    int initPPC   = 32;
    int targetPPC = 16;
    int margin    = 16;

    pp.get("init_ppc", initPPC);
    pp.get("target_ppc", targetPPC);
    pp.get("margin", margin);

    int         numRounds       = 1;
    int         stratifiedInput = 1;
    int         writeParts      = 0;
    std::string tag             = "run";
    pp.get("rounds", numRounds);
    pp.get("stratified", stratifiedInput);
    pp.get("write_particles", writeParts);
    pp.get("tag", tag);
    const bool stratified = (stratifiedInput != 0);

    // ---- Geometry-free AMR setup (no EB anywhere; the merge is the only thing under test) ----
    auto compgeom = RefCountedPtr<ComputationalGeometry>(new ComputationalGeometry());
    auto amr      = RefCountedPtr<AmrMesh>(new AmrMesh());

    amr->setMultifluidIndexSpace(compgeom->getMfIndexSpace());
    amr->sanityCheck();
    amr->buildDomains();

    Random::seed();

    auto species = RefCountedPtr<ItoSpecies>(new ItoSpecies("aliasing", 0, false, false));
    auto solver  = RefCountedPtr<ItoSolver>(new ItoSolver());

    solver->setVerbosity(-1);
    solver->parseOptions();
    solver->setAmr(amr);
    solver->setSpecies(species);
    solver->setPhase(phase::gas);
    solver->setComputationalGeometry(compgeom);
    solver->setRealm(Realm::Primal);

    amr->registerRealm(Realm::Primal);
    solver->registerOperators();

    compgeom->buildGeometries(amr->getFinestDomain(),
                              amr->getProbLo(),
                              amr->getFinestDx(),
                              amr->getMaxEbisBoxSize(),
                              amr->getNumberOfEbGhostCells(),
                              -1);

    amr->setBaseImplicitFunction(phase::gas, compgeom->getGasImplicitFunction());
    amr->setBaseImplicitFunction(phase::solid, compgeom->getSolidImplicitFunction());

    Vector<IntVectSet> tags(1 + amr->getMaxAmrDepth());
    amr->regridAmr(tags, 0);
    amr->regridOperators(0);

    solver->allocate();

    const int blockSize = amr->getMaxBoxSize();

    pout() << "===== ItoSolver particle-merge diagnostic =====" << endl;
    pout() << "domain      = " << amr->getDomains()[0].domainBox() << endl;
    pout() << "block size  = " << blockSize << endl;
    pout() << "init ppc    = " << initPPC << endl;
    pout() << "target ppc  = " << targetPPC << endl;
    pout() << "margin      = " << margin << " cells" << endl;
    pout() << "rounds      = " << numRounds << endl;
    pout() << "fill        = "
           << (stratified ? "stratified (exactly ppc per cell)" : "unstratified (Poisson over patch)") << endl;

    ParticleContainer<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

    for (int round = 1; round <= numRounds; round++) {

      // Each round adds a FRESH uniform population on top of whatever survived the previous merge,
      // then merges the union back down to the target.
      particles.organizeParticlesByPatch();
      fillUniform(particles, amr, initPPC, stratified);

      // Dump the population the merge is about to see (first round only) so the before/after pair is
      // inspectable in the same viewer.
      if (writeParts != 0 && round == 1) {
        DischargeIO::writeH5Part("parts_" + tag + "_initial.h5part", particles, amr->getProbLo(), 0.0);
      }

      Stats pre(blockSize);
      gatherStats(pre, particles, amr, margin);
      pre.reduce();

      particles.organizeParticlesByPatch();
      solver->makeSuperparticles(ItoSolver::WhichContainer::Bulk, targetPPC);

      Stats post(blockSize);
      gatherStats(post, particles, amr, margin);
      post.reduce();

      const double meanPre  = pre.m_sumW / pre.m_numCells;
      const double meanPost = post.m_sumW / post.m_numCells;
      const double varPre   = pre.m_sumW2 / pre.m_numCells - meanPre * meanPre;
      const double varPost  = post.m_sumW2 / post.m_numCells - meanPost * meanPost;

      int minN = Stats::s_maxCount;
      int maxN = 0;
      for (std::size_t i = 0; i < post.m_cellCount.size(); i++) {
        if (post.m_cellCount[i] > 0.0) {
          minN = std::min(minN, static_cast<int>(i));
          maxN = std::max(maxN, static_cast<int>(i));
        }
      }

      pout() << endl << "########## ROUND " << round << " ##########" << endl;
      pout() << "  particles  pre -> post = " << pre.m_totalParticle << " -> " << post.m_totalParticle << endl;
      pout() << "  weight     pre -> post = " << std::setprecision(12) << pre.m_totalWeight << " -> "
             << post.m_totalWeight << std::setprecision(6) << endl;
      pout() << "  regional weight change = " << (post.m_totalWeight - pre.m_totalWeight) / pre.m_totalWeight << endl;
      pout() << "  cell weight mean       = " << meanPre << " -> " << meanPost << endl;
      pout() << "  cell weight stdev      = " << std::sqrt(std::max(0.0, varPre)) << " -> "
             << std::sqrt(std::max(0.0, varPost)) << "   (rel " << std::sqrt(std::max(0.0, varPost)) / meanPost << ")"
             << endl;
      pout() << "  E[(u-0.5)^2] pre->post = " << pre.m_sumU2 / (SpaceDim * pre.m_totalParticle) << " -> "
             << post.m_sumU2 / (SpaceDim * post.m_totalParticle) << "   (uniform 0.0833333)" << endl;
      pout() << "  particles/cell range   = " << minN << " .. " << maxN << endl;

      // Block-frequency probe: per-cell density and its spread as a function of the cell's phase
      // within its 8-cell patch. A patch-local algorithm that treats patch-edge cells differently
      // from interior cells shows up here and nowhere else.
      for (int dir = 0; dir < SpaceDim; dir++) {
        pout() << "  block-phase dir " << dir << " mean:";
        for (int b = 0; b < blockSize; b++) {
          pout() << " " << std::fixed << std::setprecision(5)
                 << (post.m_blockW[dir][b] / post.m_blockN[dir][b]) / meanPost;
        }
        pout() << endl;

        pout() << "  block-phase dir " << dir << " sdev:";
        for (int b = 0; b < blockSize; b++) {
          const double m = post.m_blockW[dir][b] / post.m_blockN[dir][b];
          const double v = post.m_blockW2[dir][b] / post.m_blockN[dir][b] - m * m;
          pout() << " " << std::fixed << std::setprecision(5) << std::sqrt(std::max(0.0, v));
        }
        pout() << endl;
      }

      {
        double unmerged = 0.0;
        double meanW    = 0.0;
        for (std::size_t i = 0; i < post.m_weightHist.size(); i++) {
          meanW += static_cast<double>(i) * post.m_weightHist[i];
        }
        unmerged = post.m_weightHist[1];
        pout() << "  unmerged (w==1) frac   = " << unmerged / post.m_totalParticle
               << "   mean weight = " << meanW / post.m_totalParticle << endl;
        pout() << "  weight histogram (w:frac):";
        for (std::size_t i = 0; i < post.m_weightHist.size(); i++) {
          if (post.m_weightHist[i] / post.m_totalParticle > 0.005) {
            pout() << " " << i << ":" << std::fixed << std::setprecision(3)
                   << post.m_weightHist[i] / post.m_totalParticle;
          }
        }
        pout() << endl;
      }

      pout() << "  N_eff per cell         = " << post.m_sumNeff / post.m_numCells << " of " << targetPPC
             << "   heaviest share = " << post.m_sumTopFrac / post.m_numCells << "   weight range " << post.m_minWeight
             << " .. " << post.m_maxWeight << endl;

      pout() << "  sub-cell lattice |Z_k|, k = 1..8   (uniform floor ~ " << 1.0 / std::sqrt(post.m_totalParticle) << ")"
             << endl;
      for (int dir = 0; dir < SpaceDim; dir++) {
        pout() << "    dir " << dir << ":";
        for (int k = 0; k < Stats::s_numModes; k++) {
          const double c  = post.m_cosK[dir][k] / post.m_totalParticle;
          const double sn = post.m_sinK[dir][k] / post.m_totalParticle;
          pout() << " " << std::fixed << std::setprecision(5) << std::sqrt(c * c + sn * sn);
        }
        pout() << endl;
      }
      {
        double jsum = 0.0;
        for (const double v : post.m_joint) {
          jsum += v;
        }
        const double jflat = jsum / static_cast<double>(post.m_joint.size());
        double       jmin  = 1.E30;
        double       jmax  = 0.0;
        for (const double v : post.m_joint) {
          jmin = std::min(jmin, v / jflat);
          jmax = std::max(jmax, v / jflat);
        }
        pout() << "  joint 6^D sub-cell min/max = " << jmin << " / " << jmax << endl;
      }

      report("  subcell PRE  dir 0", pre.m_subCount[0]);

      for (int dir = 0; dir < SpaceDim; dir++) {
        report("  subcell POST dir " + std::to_string(dir), post.m_subCount[dir]);
      }

      report("  patch-mass   dir 0", post.m_patchWeight[0]);

      if (writeParts != 0 && round == numRounds) {
        DischargeIO::writeH5Part("parts_" + tag + "_final.h5part", particles, amr->getProbLo(), 1.0 * round);
      }
    }

    pout() << "===== done =====" << endl;
  }

  ChomboDischarge::finalize();
}
