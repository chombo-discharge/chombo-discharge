/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleManagementImplem.H
  @brief  Implementation of CD_ParticleManagement.H
  @author Robert Marskar
*/

#ifndef CD_PARTICLEMANAGEMENTIMPLEM_H
#define CD_PARTICLEMANAGEMENTIMPLEM_H

// Std includes
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <limits>
#include <utility>
#include <type_traits>

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_Random.H>
#include <CD_ParticleLoops.H>
#include <CD_LoadBalancing.H>
#include <CD_ParticleManagement.H>
#include <CD_NamespaceHeader.H>

namespace ParticleManagement {

  template <class P, Real& (P::*weight)(), const RealVect& (P::*position)() const>
  inline void
  buildEqualWeightKDLeaves(const std::vector<P>&                       a_particles,
                           const int                                   a_maxLeaves,
                           const BinaryParticleReconcile<P>&           a_particleReconcile,
                           std::vector<std::pair<const P*, const P*>>& a_leaves) noexcept
  {
    CH_TIME("ParticleManagement::buildEqualWeightKDLeaves");

    // Half-open particle range [lo, hi) into the flat working buffer, plus the node weight. This replaces
    // the old heap-allocated, shared_ptr-linked KDNode<P> (one alloc + atomic refcount per node).
    struct NodeRange
    {
      std::size_t lo;
      std::size_t hi;
      Real        w;
    };

    // Reusable per-thread scratch: capacity is retained across calls (cleared, never freed), so after the
    // first cell the merge does no heap allocation. The two particle buffers are ping-ponged level by level,
    // so each node's particles are always a contiguous span -- no per-node particle list.
    thread_local std::vector<P>         s_cur;
    thread_local std::vector<P>         s_nxt;
    thread_local std::vector<NodeRange> s_curNodes;
    thread_local std::vector<NodeRange> s_nxtNodes;

    a_leaves.clear();

    const std::size_t numInput = a_particles.size();
    if (numInput == 0 || a_maxLeaves <= 0) {
      return;
    }

    constexpr Real    splitThresh = 2.0 - std::numeric_limits<Real>::min();
    const std::size_t maxLeaves   = static_cast<std::size_t>(a_maxLeaves);

    // Seed the root level. Reserve room for split products: the median split adds at most one particle per
    // internal node, i.e. fewer than maxLeaves extra particles over the whole build.
    s_cur.clear();
    s_cur.reserve(numInput + maxLeaves);
    s_cur.insert(s_cur.end(), a_particles.begin(), a_particles.end());

    // Total weight (read from the non-const working copy; the weight() member-pointer is non-const-qualified).
    Real W = 0.0;
    for (P& p : s_cur) {
      W += (p.*weight)();
    }

    s_curNodes.clear();
    s_curNodes.push_back({static_cast<std::size_t>(0), numInput, W});

    // Breadth-first: split every splittable leaf of the current level, streaming the children contiguously
    // into the next buffer, until we reach maxLeaves leaves or nothing can be split further.
    bool keepGoing = true;
    while (keepGoing && s_curNodes.size() < maxLeaves) {
      keepGoing = false;

      s_nxt.clear();
      s_nxtNodes.clear();

      // Each split turns one leaf into two (+1 leaf), so this level may split at most this many nodes
      // before reaching maxLeaves. Nodes beyond the budget are carried forward unchanged -- they must NOT
      // be dropped (doing so would discard their particles and break weight conservation).
      std::size_t splitBudget = maxLeaves - s_curNodes.size();

      for (std::size_t ni = 0; ni < s_curNodes.size(); ni++) {
        const NodeRange node = s_curNodes[ni];

        if (node.w > splitThresh && splitBudget > 0) {
          splitBudget--;

          // ---- split this node (equal-weight median split) into s_nxt ----
          P* const          beg = s_cur.data() + node.lo;
          P* const          end = s_cur.data() + node.hi;
          const std::size_t n   = node.hi - node.lo;
          const Real        Wn  = node.w;

          // A. Split along the longest bounding-box extent.
          RealVect loCorner = +std::numeric_limits<Real>::max() * RealVect::Unit;
          RealVect hiCorner = -std::numeric_limits<Real>::max() * RealVect::Unit;
          for (P* p = beg; p != end; ++p) {
            const RealVect& pos = ((*p).*position)();
            for (int dir = 0; dir < SpaceDim; dir++) {
              loCorner[dir] = std::min(pos[dir], loCorner[dir]);
              hiCorner[dir] = std::max(pos[dir], hiCorner[dir]);
            }
          }
          const int splitDir = (hiCorner - loCorner).maxDir(true);

          std::sort(beg, end, [splitDir](const P& p1, const P& p2) -> bool {
            return ((p1).*position)()[splitDir] < ((p2).*position)()[splitDir];
          });

          // B. Weight-balanced median particle.
          std::size_t id = 0;
          Real        wl = 0.0;
          Real        wr = Wn - (beg[0].*weight)();
          for (std::size_t i = 1; i < n; i++) {
            const Real& w = (beg[id].*weight)();
            if (wl + w < wr) {
              id = i;
              wl += w;
              wr = Wn - wl - (beg[id].*weight)();
            }
            else {
              break;
            }
          }

          P          med = beg[id]; // non-const: weight() is a non-const-qualified member function
          const Real pw  = (med.*weight)();
          const Real dw  = wr - wl;

          CH_assert(wl + wr + pw == Wn);

          // C. Decide how the median is distributed; produce at most one extra particle per side.
          bool hasExtraL = false;
          bool hasExtraR = false;
          P    extraL;
          P    extraR;

          if (pw >= splitThresh && pw >= std::abs(dw)) {
            Real            dwl    = dw;
            Real            dwr    = 0.0;
            const Real      ddw    = pw - dw;
            const long long Nsplit = (long long)ddw;

            if (Nsplit > 0LL) {
              const long long Nr = Nsplit / 2;
              const long long Nl = Nsplit - Nr;

              dwl += (ddw / Nsplit) * Nl;
              dwr += (ddw / Nsplit) * Nr;
            }

            if (dwl > 0.0 && dwr > 0.0) {
              // Split the median particle across both children.
              extraL = med;
              extraR = med;

              CH_assert(dwl >= 1.0);
              CH_assert(dwr >= 1.0);

              wl += dwl;
              wr += dwr;

              (extraL.*weight)() = dwl;
              (extraR.*weight)() = dwr;

              a_particleReconcile(extraL, extraR, med);

              hasExtraL = true;
              hasExtraR = true;
            }
            else if (dwl > 0.0 && dwr == 0.0) {
              extraL = med;
              CH_assert(dwl >= 1.0);
              wl += dwl;
              (extraL.*weight)() = dwl;
              hasExtraL          = true;
            }
            else if (dwl == 0.0 && dwr > 0.0) {
              extraR = med;
              CH_assert(dwr >= 1.0);
              wr += dwr;
              (extraR.*weight)() = dwr;
              hasExtraR          = true;
            }
            else {
              MayDay::Abort("ParticleManagement::buildEqualWeightKDLeaves - logic bust");
            }
          }
          else {
            // Median assigned whole to the lighter side (weight unchanged).
            if (wl <= wr) {
              wl += pw;
              extraL    = med;
              hasExtraL = true;
            }
            else {
              wr += pw;
              extraR    = med;
              hasExtraR = true;
            }
          }

          CH_assert(wl + wr == Wn);
          CH_assert(std::abs(wl - wr) <= 1.0);

          // D. Stream the two children contiguously into s_nxt. Particle order within each child matches the
          //    old implementation (base half first, then the median's piece), so the leaf reductions are
          //    bit-for-bit identical. NOTE: reading from s_cur (stable) and writing to s_nxt (separate buffer).
          const std::size_t leftLo = s_nxt.size();
          s_nxt.insert(s_nxt.end(), beg, beg + id);
          if (hasExtraL) {
            s_nxt.push_back(std::move(extraL));
          }
          s_nxtNodes.push_back({leftLo, s_nxt.size(), wl});

          const std::size_t rightLo = s_nxt.size();
          s_nxt.insert(s_nxt.end(), beg + id + 1, end);
          if (hasExtraR) {
            s_nxt.push_back(std::move(extraR));
          }
          s_nxtNodes.push_back({rightLo, s_nxt.size(), wr});

          keepGoing = true;
        }
        else {
          // Cannot split further -- carry this leaf forward unchanged.
          const std::size_t lo = s_nxt.size();
          s_nxt.insert(s_nxt.end(), s_cur.begin() + node.lo, s_cur.begin() + node.hi);
          s_nxtNodes.push_back({lo, s_nxt.size(), node.w});
        }
      }

      s_cur.swap(s_nxt);
      s_curNodes.swap(s_nxtNodes);
    }

    // Emit the leaves as contiguous ranges into the (now stable) working buffer.
    a_leaves.reserve(s_curNodes.size());
    for (const NodeRange& r : s_curNodes) {
      const P* const base = s_cur.data() + r.lo;
      a_leaves.emplace_back(base, base + (r.hi - r.lo));
    }
  }

  template <class P, const RealVect& (P::*position)() const>
  inline void
  mergeAdjacentNearest(std::vector<P>&                          a_particles,
                       const std::size_t                        a_target,
                       const std::function<void(P&, const P&)>& a_combine) noexcept
  {
    CH_TIME("ParticleManagement::mergeAdjacentNearest");

    const std::size_t n = a_particles.size();
    if (a_target == 0 || n <= a_target) {
      return;
    }

    auto dist2 = [&a_particles](const std::size_t i, const std::size_t j) -> Real {
      const RealVect& pi = (a_particles[i].*position)();
      const RealVect& pj = (a_particles[j].*position)();
      Real            s  = 0.0;
      for (int d = 0; d < SpaceDim; d++) {
        const Real dd = pi[d] - pj[d];
        s += dd * dd;
      }
      return s;
    };

    // Index-based doubly-linked chain over the (spatially-ordered) list + liveness + version tags. No
    // pointers and no per-call allocation after warmup.
    thread_local std::vector<int>      prev, next;
    thread_local std::vector<uint8_t>  alive;
    thread_local std::vector<uint32_t> ver;
    prev.assign(n, -1);
    next.assign(n, -1);
    alive.assign(n, 1);
    ver.assign(n, 0);
    for (std::size_t i = 0; i < n; i++) {
      prev[i] = (i == 0) ? -1 : static_cast<int>(i - 1);
      next[i] = (i + 1 == n) ? -1 : static_cast<int>(i + 1);
    }

    struct Edge
    {
      Real     d;
      int      a;
      int      b;
      uint32_t va;
      uint32_t vb;
    };
    const auto cmp = [](const Edge& x, const Edge& y) {
      return x.d > y.d;
    }; // min-heap
    thread_local std::vector<Edge> heap;
    heap.clear();
    for (std::size_t i = 0; i + 1 < n; i++) {
      heap.push_back({dist2(i, i + 1), static_cast<int>(i), static_cast<int>(i + 1), 0, 0});
    }
    std::make_heap(heap.begin(), heap.end(), cmp);

    std::size_t count = n;
    while (count > a_target && !heap.empty()) {
      std::pop_heap(heap.begin(), heap.end(), cmp);
      const Edge edge = heap.back();
      heap.pop_back();

      const int a = edge.a;
      const int b = edge.b;
      if (!alive[a] || !alive[b] || ver[a] != edge.va || ver[b] != edge.vb || next[a] != b) {
        continue; // stale entry
      }

      // Combine b into a; a's centroid (and one neighbour) changes, so bump a's version.
      a_combine(a_particles[a], a_particles[b]);
      alive[b]     = 0;
      const int nb = next[b];
      next[a]      = nb;
      if (nb >= 0) {
        prev[nb] = a;
      }
      count--;
      ver[a]++;

      if (prev[a] >= 0) {
        heap.push_back(
          {dist2(static_cast<std::size_t>(prev[a]), static_cast<std::size_t>(a)), prev[a], a, ver[prev[a]], ver[a]});
        std::push_heap(heap.begin(), heap.end(), cmp);
      }
      if (next[a] >= 0) {
        heap.push_back(
          {dist2(static_cast<std::size_t>(a), static_cast<std::size_t>(next[a])), a, next[a], ver[a], ver[next[a]]});
        std::push_heap(heap.begin(), heap.end(), cmp);
      }
    }

    // Compact the survivors back into a_particles (preserving chain order).
    thread_local std::vector<P> survivors;
    survivors.clear();
    survivors.reserve(count);
    for (std::size_t i = 0; i < n; i++) {
      if (alive[i]) {
        survivors.push_back(a_particles[i]);
      }
    }
    a_particles.swap(survivors);
  }

  template <class Packed,
            Real& (Packed::*packWeight)(),
            const RealVect& (Packed::*packPosition)() const,
            class P,
            class Traits>
  inline ParticleMerger<P, Traits>
  makeSfcNearestNeighborMerger(std::function<Packed(const ParticleSoA<P, Traits>&, std::size_t)> a_gather,
                               std::function<void(Packed&, const Packed&)>                       a_combine,
                               std::function<void(ParticleSoA<P, Traits>&, const Packed&)>       a_scatter) noexcept
  {
    return [=](ParticleSoA<P, Traits>& a_particles, const CellInfo& /*a_cellInfo*/, int a_ppc) noexcept {
      CH_TIMERS("ParticleManagement::makeSfcNearestNeighborMerger");
      CH_TIMER("ParticleManagement::makeSfcNearestNeighborMerger::populate", t1);
      CH_TIMER("ParticleManagement::makeSfcNearestNeighborMerger::sort_merge", t2);
      CH_TIMER("ParticleManagement::makeSfcNearestNeighborMerger::emit", t3);

      CH_START(t1);
      thread_local std::vector<Packed> particles;
      particles.clear();
      particles.reserve(a_particles.size());
      Real W = 0.0;
      for (std::size_t i = 0; i < a_particles.size(); i++) {
        Packed p = a_gather(a_particles, i);
        W += (p.*packWeight)();
        particles.emplace_back(std::move(p));
      }
      CH_STOP(t1);

      if (W < 2.0 || a_ppc <= 0) {
        return;
      }

      const std::size_t target = static_cast<std::size_t>(a_ppc);

      CH_START(t2);
      if (particles.size() > target) {
        RealVect loCorner = +std::numeric_limits<Real>::max() * RealVect::Unit;
        RealVect hiCorner = -std::numeric_limits<Real>::max() * RealVect::Unit;
        for (const Packed& p : particles) {
          const RealVect& x = (p.*packPosition)();
          for (int dir = 0; dir < SpaceDim; dir++) {
            loCorner[dir] = std::min(loCorner[dir], x[dir]);
            hiCorner[dir] = std::max(hiCorner[dir], x[dir]);
          }
        }

        RealVect invExtent = RealVect::Zero;
        for (int dir = 0; dir < SpaceDim; dir++) {
          invExtent[dir] = (hiCorner[dir] > loCorner[dir]) ? 1.0 / (hiCorner[dir] - loCorner[dir]) : 0.0;
        }

        constexpr uint32_t                                         maxCoord = (1U << 21) - 1U;
        thread_local std::vector<std::pair<uint64_t, std::size_t>> keyed;
        keyed.clear();
        keyed.reserve(particles.size());
        for (std::size_t i = 0; i < particles.size(); i++) {
          const RealVect&                x = (particles[i].*packPosition)();
          std::array<uint32_t, SpaceDim> coords;
          for (int dir = 0; dir < SpaceDim; dir++) {
            Real t      = (x[dir] - loCorner[dir]) * invExtent[dir];
            t           = std::max(0.0, std::min(1.0, t));
            coords[dir] = static_cast<uint32_t>(t * maxCoord);
          }
          keyed.emplace_back(LoadBalancing::hilbertIndex<SpaceDim>(coords), i);
        }
        std::sort(keyed.begin(), keyed.end(), [](const auto& a, const auto& b) {
          return a.first < b.first;
        });

        thread_local std::vector<Packed> sorted;
        sorted.clear();
        sorted.reserve(particles.size());
        for (const auto& k : keyed) {
          sorted.push_back(particles[k.second]);
        }
        particles.swap(sorted);

        mergeAdjacentNearest<Packed, packPosition>(particles, target, a_combine);
      }
      else if (particles.size() < target) {
        while (particles.size() < target) {
          std::size_t hi = 0;
          for (std::size_t i = 1; i < particles.size(); i++) {
            if ((particles[i].*packWeight)() > (particles[hi].*packWeight)()) {
              hi = i;
            }
          }
          if ((particles[hi].*packWeight)() < 2.0) {
            break;
          }
          Packed half = particles[hi];
          (particles[hi].*packWeight)() *= 0.5;
          (half.*packWeight)() *= 0.5;
          particles.push_back(half);
        }
      }
      CH_STOP(t2);

      CH_START(t3);
      a_particles.clear();
      for (const Packed& p : particles) {
        a_scatter(a_particles, p);
      }
      CH_STOP(t3);
    };
  }

  template <typename P, typename Traits, typename T, typename>
  inline void
  removePhysicalParticles(ParticleSoA<P, Traits>& a_particles, const T a_numPhysPartToRemove) noexcept
  {
    CH_TIME("ParticleManagement::removePhysicalParticles(SoA)");

    constexpr T zero = (T)0;

    if (a_numPhysPartToRemove < zero) {
      MayDay::Error("ParticleManagement::removePhysicalParticles(SoA) - 'a_numPhysPartoToRemove < 0'");
    }

    const std::size_t numComp = a_particles.size();

    if (numComp > 0) {
      T numRemoved = zero;

      // 1. Compute the minimum particle weight.
      T minWeight = std::numeric_limits<T>::max();
      for (std::size_t i = 0; i < numComp; i++) {
        minWeight = std::min(minWeight, (T)a_particles.weight(i));
      }

      // 2. Trim particle weights down to minWeight.
      for (std::size_t i = 0; i < numComp; i++) {
        const T diff1 = (T)a_particles.weight(i) - minWeight;
        const T diff2 = a_numPhysPartToRemove - numRemoved;

        CH_assert(diff1 >= zero);
        CH_assert(diff2 >= zero);

        const T r = std::max(0LL, std::min(diff1, diff2));

        a_particles.weight(i) -= 1.0 * r;
        numRemoved += r;
      }

      // 3. "Uniformly" subtract the particle weights.
      if (a_numPhysPartToRemove - numRemoved > zero) {
        const T numCompParticles = (T)numComp;
        const T uniformWeight    = (a_numPhysPartToRemove - numRemoved) / numCompParticles;
        const T uniformRemainder = (a_numPhysPartToRemove - numRemoved) % numCompParticles;

        if (uniformWeight > zero) {
          // uniformWeight is constant over the loop, so the running accumulation equals
          // uniformWeight * numComp -- hoist it out so the body is a pure elementwise update.
          double* const w = a_particles.weightColumn();

          ParticleLoops::loop(a_particles, [&](const std::size_t i) {
            w[i] -= 1.0 * uniformWeight;
          });

          numRemoved += uniformWeight * static_cast<T>(numComp);
        }

        if (uniformRemainder > zero) {
          T W = 0;

          for (std::size_t i = 0; i < numComp; i++) {

            // Never remove so that weight is negative.
            const T w = std::min((T)a_particles.weight(i), uniformRemainder - W);

            a_particles.weight(i) -= 1.0 * w;

            W += w;
            numRemoved += w;

            if (W == uniformRemainder) {
              break;
            }
          }
        }
      }

      CH_assert(numRemoved == a_numPhysPartToRemove);
    }
  }

  template <typename P, typename Traits>
  inline void
  deleteParticles(ParticleSoA<P, Traits>& a_particles, const Real a_weightThresh) noexcept
  {
    CH_TIME("ParticleManagement::deleteParticles(SoA)");

    // Swap-and-pop: do not advance the index after a removal (a new particle now occupies slot i).
    std::size_t i = 0;
    while (i < a_particles.size()) {
      if (a_particles.weight(i) < a_weightThresh) {
        a_particles.remove(i);
      }
      else {
        i++;
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
        T       r = a_numPhysicalParticles % a_maxCompParticles;

        if (W > zero) {
          ret.resize(a_maxCompParticles, W);

          for (int i = 0; i < ret.size() && r > zero; i++) {
            ret[i] += one;
            r--;
          }
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

  template <typename P, typename Traits, typename T, typename>
  inline void
  drawRandomParticles(ParticleSoA<P, Traits>&          a_particles,
                      const T                          a_numParticles,
                      const std::function<RealVect()>& a_distribution)
  {
    a_particles.clear();

    // Each rank draws its own share (unit weight) into the buffer. Routing to owners happens later when
    // the buffer is added to a ParticleContainer.
    const T numParticles = partitionParticles(a_numParticles);

    for (T t = 0; t < numParticles; t++) {
      a_particles.append(a_distribution(), 1.0);
    }
  }

  template <typename P, typename Traits, typename T, typename>
  inline void
  drawGaussianParticles(ParticleSoA<P, Traits>& a_particles,
                        const T                 a_numParticles,
                        const RealVect&         a_center,
                        const Real              a_radius) noexcept
  {
    CH_TIME("ParticleManagement::drawGaussianParticles(SoA)");

    std::normal_distribution<Real> gauss(0.0, a_radius);

    auto ranGauss = [&]() -> RealVect {
      return a_center + Random::get(gauss) * Random::getDirection();
    };

    drawRandomParticles(a_particles, a_numParticles, ranGauss);
  }

  template <typename P, typename Traits, typename T, typename>
  inline void
  drawBoxParticles(ParticleSoA<P, Traits>& a_particles,
                   const T                 a_numParticles,
                   const RealVect&         a_loCorner,
                   const RealVect&         a_hiCorner) noexcept
  {
    CH_TIME("ParticleManagement::drawBoxParticles(SoA)");

    CH_assert(a_hiCorner >= a_loCorner);

    auto ranBox = [&]() -> RealVect {
      return RealVect(D_DECL(a_loCorner[0] + (a_hiCorner[0] - a_loCorner[0]) * Random::getUniformReal01(),
                             a_loCorner[1] + (a_hiCorner[1] - a_loCorner[1]) * Random::getUniformReal01(),
                             a_loCorner[2] + (a_hiCorner[2] - a_loCorner[2]) * Random::getUniformReal01()));
    };

    drawRandomParticles(a_particles, a_numParticles, ranBox);
  }

  template <typename P, typename Traits, typename T, typename>
  inline void
  drawSphereParticles(ParticleSoA<P, Traits>& a_particles,
                      const T                 a_numParticles,
                      const RealVect&         a_center,
                      const Real              a_radius) noexcept
  {
    CH_TIME("ParticleManagement::drawSphereParticles(SoA)");

    auto ranSphere = [&]() -> RealVect {
      RealVect x = std::numeric_limits<Real>::max() * RealVect::Unit;

      while (x.vectorLength() > a_radius) {
        for (int d = 0; d < SpaceDim; d++) {
          x[d] = a_radius * Random::getUniformReal11();
        }
      }

      return x + a_center;
    };

    drawRandomParticles(a_particles, a_numParticles, ranSphere);
  }
} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#endif
