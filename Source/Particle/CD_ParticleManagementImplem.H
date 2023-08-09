/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ParticleManagementImplem.H
  @brief  Implementation of CD_ParticleManagement.H
  @author Robert Marskar
*/

#ifndef CD_ParticleManagementImplem_H
#define CD_ParticleManagementImplem_H

// Std includes
#include <utility>
#include <type_traits>

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_Random.H>
#include <CD_ParticleManagement.H>
#include <CD_NamespaceHeader.H>

namespace ParticleManagement {

  template <class P, Real& (P::*weight)(), const RealVect& (P::*position)() const>
  static inline void
  partitionAndSplitEqualWeightKD(KDNode<P>& a_node, const BinaryParticleReconcile<P> a_particleReconcile) noexcept
  {
    CH_assert(!(a_node.isInteriorNode()));
    CH_assert(a_node.weight() > 2.0 - std::numeric_limits<Real>::min());

    constexpr Real splitThresh = 2.0 - std::numeric_limits<Real>::min();

    // Particles and node weight.
    typename KDNode<P>::ParticleList& particles = a_node.getParticles();

    const Real W = a_node.weight();

    // A. Figure out which coordinate direction we should partition and sort
    //    the particles.
    RealVect loCorner = +std::numeric_limits<Real>::max() * RealVect::Unit;
    RealVect hiCorner = -std::numeric_limits<Real>::max() * RealVect::Unit;

    for (auto& p : particles) {
      const RealVect& pos = (p.*position)();
      for (int dir = 0; dir < SpaceDim; dir++) {
        loCorner[dir] = std::min(pos[dir], loCorner[dir]);
        hiCorner[dir] = std::max(pos[dir], hiCorner[dir]);
      }
    }

    const int splitDir = (hiCorner - loCorner).maxDir(true);

    auto sortCrit = [splitDir](const P& p1, const P& p2) -> bool {
      return (p1.*position)()[splitDir] < (p2.*position)()[splitDir];
    };

    std::sort(particles.begin(), particles.end(), sortCrit);

    // B. Determine the "median particle" and start computing the weight in the
    //    two halves.
    size_t id = 0;
    Real   wl = 0.0;
    Real   wr = W - (particles[id].*weight)();

    for (size_t i = 1; i < particles.size(); i++) {
      const Real& w = (particles[id].*weight)();

      if (wl + w < wr) {
        id = i;
        wl += w;
        wr = W - wl - (particles[id].*weight)();
      }
      else {
        break;
      }
    }

    // C. Copy the two particle halves to each subnode.
    P p = particles[id];

    typename KDNode<P>::ParticleList pl;
    typename KDNode<P>::ParticleList pr;

    std::move(particles.begin(), particles.begin() + id, std::back_inserter(pl));
    std::move(particles.begin() + id + 1, particles.end(), std::back_inserter(pr));

    const Real& pw = (p.*weight)();
    const Real  dw = wr - wl;

    CH_assert(wl + wr + pw == W);

    // D. Assign the median particle; split the particle if we can.
    if (pw >= splitThresh && pw >= std::abs(dw)) {
      Real dwl = dw;
      Real dwr = 0.0;
      Real ddw = pw - dw;

      const long long N = (long long)ddw;

      if (N > 0LL) {

        const long long Nr = N / 2;
        const long long Nl = N - Nr;

        dwl += (ddw / N) * Nl;
        dwr += (ddw / N) * Nr;
      }

      if (dwl > 0.0 && dwr > 0.0) {
        // Splitting particle.

        P il(p);
        P ir(p);

        CH_assert(dwl >= 1.0);
        CH_assert(dwr >= 1.0);

        wl += dwl;
        wr += dwr;

        (il.*weight)() = dwl;
        (ir.*weight)() = dwr;

        // User can reconcile other particle properties.
        a_particleReconcile(il, ir, p);

        pl.emplace_back(std::move(il));
        pr.emplace_back(std::move(ir));
      }
      else if (dwl > 0.0 && dwr == 0.0) {
        // Particle assigned to left node.

        P il(p);

        CH_assert(dwl >= 1.0);

        wl += dwl;
        (il.*weight)() = dwl;
        pl.emplace_back(std::move(il));
      }
      else if (dwl == 0.0 && dwr > 0.0) {
        // Particle assigned to right node.

        P ir(p);

        CH_assert(dwr >= 1.0);

        wr += dwr;
        (ir.*weight)() = dwr;
        pr.emplace_back(std::move(ir));
      }
      else {
        // Should not happen.

        MayDay::Abort("ParticleManagement::partitionAndSplitEqualWeightKD - logic bust");
      }
    }
    else {
      if (wl <= wr) {
        wl += pw;
        pl.emplace_back(std::move(p));
      }
      else {
        wr += pw;
        pr.emplace_back(std::move(p));
      }
    }

    // E. If this breaks, weight is not conserved or we broke the median particle splitting; the weight difference
    //    between the left/right node should be at most one physical particle.
    CH_assert(wl + wr == W);
    CH_assert(std::abs(wl - wr) <= 1.0);

    // F. Instantiate the child nodes.
    particles.resize(0);

    a_node.getLeft()  = std::make_shared<KDNode<P>>(pl);
    a_node.getRight() = std::make_shared<KDNode<P>>(pr);

    a_node.getLeft()->weight()  = wl;
    a_node.getRight()->weight() = wr;

    // G. Debug code; make sure particle weights make sense
#ifndef NDEBUG
    Real WL = 0.0;
    Real WR = 0.0;

    // Note: Not auto& l/r : pl/pr because the particles were into the child node.
    for (auto& l : a_node.getLeft()->getParticles()) {
      CH_assert((l.*weight)() >= 1.0);

      WL += (l.*weight)();
    }

    for (auto& r : a_node.getRight()->getParticles()) {
      CH_assert((r.*weight)() >= 1.0);

      WR += (r.*weight)();
    }

    CH_assert(WL == wl);
    CH_assert(WR == wr);
#endif
  }

  template <class P, Real& (P::*weight)(), const RealVect& (P::*position)() const>
  inline std::vector<std::shared_ptr<KDNode<P>>>
  recursivePartitionAndSplitEqualWeightKD(typename KDNode<P>::ParticleList& a_inputParticles,
                                          const int                         a_maxLeaves,
                                          const BinaryParticleReconcile<P>  a_particleReconcile) noexcept
  {
    CH_TIME("ParticleManagement::recursivePartitionAndSplitEqualWeightKD");

    Real W = 0.0;

    for (auto& p : a_inputParticles) {
      W += (p.*weight)();
    }

    std::vector<std::shared_ptr<KDNode<P>>> leaves;

    leaves.emplace_back(std::make_shared<KDNode<P>>(a_inputParticles));
    leaves[0]->weight() = W;

    bool keepGoing = true;

    while (keepGoing && leaves.size() < a_maxLeaves) {
      keepGoing = false;

      std::vector<std::shared_ptr<KDNode<P>>> newLeaves;

      for (const auto& l : leaves) {
        if (l->weight() > 2.0 - std::numeric_limits<Real>::min()) {
          ParticleManagement::partitionAndSplitEqualWeightKD<P, weight, position>(*l, a_particleReconcile);

          newLeaves.emplace_back(l->getLeft());
          newLeaves.emplace_back(l->getRight());

          keepGoing = true;
        }
        else {
          newLeaves.emplace_back(l);
        }

        // Break out if we have sufficient leaf nodes.
        if (newLeaves.size() >= a_maxLeaves) {
          break;
        }
      }

      leaves = newLeaves;
    }

    return leaves;
  }

  template <typename P, typename T, typename>
  inline void
  removePhysicalParticles(List<P>& a_particles, const T a_numPhysPartToRemove) noexcept
  {
    CH_TIME("ParticleManagement::removePhysicalParticles");

    constexpr T zero = (T)0;

    // Obviously an error that the user should catch.
    if (a_numPhysPartToRemove < zero) {
      MayDay::Error("ParticleManagement::removePhysicalParticles - 'a_numPhysPartoToRemove < 0'");
    }

    if (a_particles.length() > 0) {
      ListIterator<P> lit(a_particles);

      T numRemoved = zero;

      // 1. Compute the minimum particle weight.
      T minWeight = std::numeric_limits<T>::max();
      for (lit.begin(); lit.ok(); ++lit) {
        minWeight = std::min(minWeight, (T)lit().weight());
      }

      // 2. Trim particle weights down to minWeight.
      for (lit.begin(); lit.ok(); ++lit) {
        const T diff1 = (T)lit().weight() - minWeight;
        const T diff2 = a_numPhysPartToRemove - numRemoved;

        CH_assert(diff1 >= zero);
        CH_assert(diff2 >= zero);

        const T r = std::max(0LL, std::min(diff1, diff2));

        lit().weight() -= 1.0 * r;
        numRemoved += r;
      }

      // 3. "Uniformly" subtract the particle weights.
      if (a_numPhysPartToRemove - numRemoved > zero) {
        const T numCompParticles = (T)a_particles.length();
        const T uniformWeight    = (a_numPhysPartToRemove - numRemoved) / numCompParticles;
        const T uniformRemainder = (a_numPhysPartToRemove - numRemoved) % numCompParticles;

        // Uniformly remove weight from each particle.
        if (uniformWeight > zero) {
          for (lit.begin(); lit.ok(); ++lit) {
            lit().weight() -= 1.0 * uniformWeight;

            numRemoved += uniformWeight;
          }
        }

        // May have to remove remainder from multiple particles because their weights might be less
        // then the actual remainder.
        if (uniformRemainder > zero) {
          T W = 0;

          for (lit.begin(); lit.ok(); ++lit) {

            // Never remove so that weight is negative.
            const T w = std::min((T)lit().weight(), uniformRemainder - W);

            lit().weight() -= 1.0 * w;

            W += w;
            numRemoved += w;

            if (W == uniformRemainder) {
              break;
            }
          }
        }
      }

      // Debug code.
      CH_assert(numRemoved == a_numPhysPartToRemove);
#ifdef NDEBUG
      for (lit.begin(); lit.ok(); ++lit) {
        CH_assert(lit().weight() >= 0.0);
      }
#endif
    }
  }

  template <typename P>
  inline void
  deleteParticles(List<P>& a_particles, const Real a_weightThresh) noexcept
  {
    CH_TIME("ParticleManagement::deleteParticles");

    for (ListIterator<P> lit(a_particles); lit.ok();) {
      if (lit().weight() < a_weightThresh) {
        a_particles.remove(lit);
      }
      else {
        ++lit;
      }
    }
  }

  template <typename T, typename>
  inline std::vector<T>
  partitionParticleWeights(const T a_numPhysicalParticles, const T a_maxCompParticles) noexcept
  {
    std::vector<T> ret(0);

    constexpr T zero = (T)0;
    constexpr T one  = (T)1;

    if (a_maxCompParticles > zero) {
      if (a_numPhysicalParticles <= a_maxCompParticles) {
        ret.resize(a_numPhysicalParticles, one);
      }
      else {
        const T W = a_numPhysicalParticles / a_maxCompParticles;
        const T r = a_numPhysicalParticles % a_maxCompParticles;

        if (W > zero) {
          ret.resize(a_maxCompParticles, W);

          ret.back() += r;
        }
        else {
          ret.resize(1, r);
        }
      }
    }

    return ret;
  }

  template <typename T, typename>
  inline T
  partitionParticles(const T a_numParticles)
  {
#ifdef CH_MPI
    const T quotient  = a_numParticles / numProc();
    const T remainder = a_numParticles % numProc();

    Vector<T> particlesPerRank(numProc(), quotient);

    for (int i = 0; i < remainder; i++) {
      particlesPerRank[i]++;
    }

    return particlesPerRank[procID()];
#else
    return a_numParticles;
#endif
  }

  template <typename P, typename T, typename>
  inline void
  drawRandomParticles(List<P>& a_particles, const T a_numParticles, const std::function<RealVect()>& a_distribution)
  {

    a_particles.clear();

    const T numParticles = partitionParticles(a_numParticles);

    for (T t = 0; t < numParticles; t++) {
      P p;
      p.position() = a_distribution();

      a_particles.add(p);
    }
  }

  template <typename P, typename T, typename>
  inline void
  drawSphereParticles(List<P>&       a_particles,
                      const T        a_numParticles,
                      const RealVect a_center,
                      const Real     a_radius) noexcept
  {
    CH_TIME("ParticleManagement::drawSphereParticles");

    a_particles.clear();

    const T numParticles = partitionParticles(a_numParticles);

    for (T t = 0; t < numParticles; t++) {
      P p;

      RealVect& x = p.position();

      x = std::numeric_limits<Real>::max() * RealVect::Unit;

      while (x.vectorLength() > a_radius) {
        for (int d = 0; d < SpaceDim; d++) {
          x[d] = a_radius * Random::getUniformReal11();
        }
      }

      x += a_center;

      a_particles.add(p);
    }
  }

  template <typename P, typename T, typename>
  inline void
  drawBoxParticles(List<P>&       a_particles,
                   const T        a_numParticles,
                   const RealVect a_loCorner,
                   const RealVect a_hiCorner) noexcept
  {
    CH_TIME("ParticleManagement::drawBoxParticles");

    CH_assert(a_hiCorner >= a_loCorner);

    auto ranBox = [&]() -> RealVect {
      return RealVect(D_DECL(a_loCorner[0] + (a_hiCorner[0] - a_loCorner[0]) * Random::getUniformReal01(),
                             a_loCorner[1] + (a_hiCorner[1] - a_loCorner[1]) * Random::getUniformReal01(),
                             a_loCorner[2] + (a_hiCorner[2] - a_loCorner[2]) * Random::getUniformReal01()));
    };

    drawRandomParticles(a_particles, a_numParticles, ranBox);
  }

  template <typename P, typename T, typename>
  inline void
  drawGaussianParticles(List<P>&       a_particles,
                        const T        a_numParticles,
                        const RealVect a_center,
                        const Real     a_radius) noexcept
  {
    CH_TIME("ParticleManagement::drawGaussianParticles");

    std::normal_distribution<Real> gauss(0.0, a_radius);

    auto ranGauss = [&]() -> RealVect {
      return a_center + Random::get(gauss) * Random::getDirection();
    };

    drawRandomParticles(a_particles, a_numParticles, ranGauss);
  }
} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#endif
