/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SuperParticlesImplem.H
  @brief  Implementation of CD_SuperParticles.H
  @author Robert Marskar
*/

#ifndef CD_SuperParticlesImplem_H
#define CD_SuperParticlesImplem_H

// Std includes
#include <utility>

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_SuperParticles.H>
#include <CD_NamespaceHeader.H>

namespace SuperParticles {

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

      const long long N = llround(ddw);

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
} // namespace SuperParticles

#include <CD_NamespaceFooter.H>

#endif
