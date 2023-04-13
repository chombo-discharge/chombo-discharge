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

  template <class P>
  inline KDNode<P>::KDNode() : m_left(nullptr), m_right(nullptr), m_weight(0.0)
  {}

  template <class P>
  inline KDNode<P>::KDNode(ParticleList& a_particles)
    : m_left(nullptr), m_right(nullptr), m_weight(0.0), m_particles(std::move(a_particles))
  {}

  template <class P>
  inline KDNode<P>::~KDNode()
  {}

  template <class P>
  inline const typename KDNode<P>::ParticleList&
  KDNode<P>::getParticles() const noexcept
  {
    return m_particles;
  }

  template <class P>
  inline typename KDNode<P>::ParticleList&
  KDNode<P>::getParticles() noexcept
  {
    return m_particles;
  }

  template <class P>
  inline const Real&
  KDNode<P>::weight() const noexcept
  {
    return m_weight;
  }

  template <class P>
  inline Real&
  KDNode<P>::weight() noexcept
  {
    return m_weight;
  }

  template <class P>
  inline bool
  KDNode<P>::isLeafNode() const noexcept
  {
    return (m_left == nullptr) && (m_right == nullptr);
  }

  template <class P>
  inline bool
  KDNode<P>::isInteriorNode() const noexcept
  {
    return !(this->isLeafNode());
  }

  template <class P>
  inline typename KDNode<P>::ParticleList
  KDNode<P>::gatherParticles() const noexcept
  {
    ParticleList primitives;

    this->gatherParticles(primitives);

    return primitives;
  }

  template <class P>
  inline void
  KDNode<P>::gatherParticles(ParticleList& a_particles) const noexcept
  {
    if (this->isLeafNode()) {
      a_particles.reserve(a_particles.size() + m_particles.size());
      a_particles.insert(a_particles.end(), m_particles.begin(), m_particles.end());
    }
    else {
      m_left->gatherParticles(a_particles);
      m_right->gatherParticles(a_particles);
    }
  }

  template <class P>
  inline typename KDNode<P>::ParticleList
  KDNode<P>::moveParticles() noexcept
  {
    ParticleList primitives;

    this->moveParticles(primitives);

    return primitives;
  }

  template <class P>
  inline void
  KDNode<P>::moveParticles(ParticleList& a_particles) noexcept
  {
    if (this->isLeafNode()) {
      a_particles.reserve(a_particles.size() + m_particles.size());

      std::move(m_particles.begin(), m_particles.end(), std::back_inserter(a_particles));
    }
    else {
      m_left->moveParticles(a_particles);
      m_right->moveParticles(a_particles);
    }
  }

  template <class P>
  inline std::shared_ptr<KDNode<P>>&
  KDNode<P>::getLeft() noexcept
  {
    return m_left;
  }

  template <class P>
  inline std::shared_ptr<KDNode<P>>&
  KDNode<P>::getRight() noexcept
  {
    return m_right;
  }

  /*!
    @brief Partitioning method for KDNode<P> which equilibriates the weight across the two halves when splitting a node
    @note This only adjusts the weight of the median particle. The particle class P should have a copy constructor for the other
    fields; note that these fields are simply copied over when the particle is split and are not adjusted. 
  */
  template <class P, Real& (P::*weight)(), const RealVect& (P::*position)() const>
  typename KDNode<P>::Partitioner PartitionEqualWeight = [](KDNode<P>& a_node) -> void {
    CH_assert(!(a_node.isInteriorNode()));

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

      if (dwl > 0.0) {
        P il(p);

        CH_assert(dwl >= 1.0);

        wl += dwl;
        (il.*weight)() = dwl;
        pl.emplace_back(std::move(il));
      }

      if (dwr > 0.0) {
        P ir(p);

        CH_assert(dwr >= 1.0);

        wr += dwr;
        (ir.*weight)() = dwr;
        pr.emplace_back(std::move(ir));
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
  };

  /*!
    @brief Stop function which returns true if weight in node < 2
    @note The particle type is P and P::*weight is a pointer to a member function providing const-reference access to the particle weight.
  */
  template <class P, const Real& (P::*weight)() const>
  typename KDNode<P>::StopFunction InsufficientWeightForSplitting =
    [](const typename KDNode<P>::ParticleList& a_particles) {
      constexpr Real splitThresh = 2.0 - std::numeric_limits<Real>::min();

      Real W = 0.0;
      for (const auto& p : a_particles) {
        W += (p.*weight)();
      }

      return (W < splitThresh) ? true : false;
    };

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

    const T numParticles = ParticleManagement::partitionParticles(a_numParticles);

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

    const T numParticles = ParticleManagement::partitionParticles(a_numParticles);

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

    ParticleManagement::drawRandomParticles(a_particles, a_numParticles, ranBox);
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

    ParticleManagement::drawRandomParticles(a_particles, a_numParticles, ranGauss);
  }
} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#endif
