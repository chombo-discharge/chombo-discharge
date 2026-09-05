/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ParticleManagementImplem.H
 * @brief  Implementation of CD_ParticleManagement.H
 * @author Robert Marskar
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
#include <CD_ParticleManagement.H>
#include <CD_NamespaceHeader.H>

namespace ParticleManagement {

inline ParticleMergeMethod
mergeMethodFromString(const std::string& a_str) noexcept
{
  // Legacy selectors. Each one is now a method plus a specifier, or is gone; abort with the replacement
  // spelled out rather than silently picking one of the axes for the user.
  const auto replaced = [&a_str](const std::string& a_replacement) -> void {
    const std::string msg = "ParticleManagement::mergeMethodFromString - '" + a_str + "' has been replaced by " +
                            a_replacement;

    MayDay::Abort(msg.c_str());
  };

  if (a_str == "none") {
    return ParticleMergeMethod::None;
  }
  else if (a_str == "kd_cell") {
    return ParticleMergeMethod::KdCell;
  }
  else if (a_str == "kd_patch") {
    return ParticleMergeMethod::KdPatch;
  }
  else if (a_str == "kd_amr") {
    return ParticleMergeMethod::KdAmr;
  }
  else if (a_str == "nn_amr") {
    return ParticleMergeMethod::NnAmr;
  }
  else if (a_str == "reinitialize") {
    return ParticleMergeMethod::Reinitialize;
  }
  else if (a_str == "external") {
    return ParticleMergeMethod::External;
  }
  else if (a_str == "kd_carve") {
    replaced("'kd_amr' with ItoSolver.kd_amr_boundary = carve");

    return ParticleMergeMethod::KdAmr;
  }
  else if (a_str == "kd_skin_nn") {
    replaced("'kd_amr' with ItoSolver.kd_amr_boundary = nn");

    return ParticleMergeMethod::KdAmr;
  }
  else if (a_str == "nn_pair_tree") {
    replaced("'nn_amr' with ItoSolver.nn_search = tree");

    return ParticleMergeMethod::NnAmr;
  }
  else if (a_str == "nn_pair_hash") {
    replaced("'nn_amr' with ItoSolver.nn_search = hash");

    return ParticleMergeMethod::NnAmr;
  }
  else if (a_str == "nn_pair_onecell") {
    replaced("'nn_amr' with ItoSolver.nn_search = onecell");

    return ParticleMergeMethod::NnAmr;
  }
  else if (a_str == "nn_sfc") {
    // Removed: a Hilbert sort followed by repeated nearest-adjacent-pair merging inside one cell. It
    // was strictly worse than the per-cell kd merge on every axis that was measured, so it is gone
    // rather than deprecated.
    MayDay::Abort("ParticleManagement::mergeMethodFromString - 'nn_sfc' has been removed; use 'kd_cell'");

    return ParticleMergeMethod::None;
  }
  else if (a_str == "equal_weight_kd" || a_str == "equal_weight_kd_sampled" || a_str == "reinitialize_bvh") {
    replaced("'kd_cell' with ItoSolver.kd_partition = weight and ItoSolver.kd_placement = " +
             std::string(a_str == "equal_weight_kd" ? "centroid"
                                                    : (a_str == "equal_weight_kd_sampled" ? "sample" : "random")));

    return ParticleMergeMethod::KdCell;
  }
  else {
    MayDay::Abort(("ParticleManagement::mergeMethodFromString - unknown merge method '" + a_str + "'").c_str());

    return ParticleMergeMethod::None;
  }
}

inline KDPartition
kdPartitionFromString(const std::string& a_str) noexcept
{
  if (a_str == "weight") {
    return KDPartition::Weight;
  }
  else if (a_str == "count") {
    return KDPartition::Count;
  }
  else if (a_str == "hybrid") {
    return KDPartition::Hybrid;
  }
  else if (a_str == "weight_capped") {
    return KDPartition::WeightCapped;
  }
  else {
    const std::string msg = "ParticleManagement::kdPartitionFromString - unknown partition '" + a_str +
                            "' (expected 'weight', 'count', 'hybrid' or 'weight_capped')";
    MayDay::Abort(msg.c_str());
  }

  return KDPartition::Weight;
}

inline KDPlacement
kdPlacementFromString(const std::string& a_str) noexcept
{
  if (a_str == "centroid") {
    return KDPlacement::Centroid;
  }
  else if (a_str == "sample") {
    return KDPlacement::Sample;
  }
  else if (a_str == "random") {
    return KDPlacement::Random;
  }
  else {
    const std::string msg = "ParticleManagement::kdPlacementFromString - unknown placement '" + a_str +
                            "' (expected 'centroid', 'sample' or 'random')";
    MayDay::Abort(msg.c_str());
  }

  return KDPlacement::Centroid;
}

inline KDSplitPlacement
kdSplitPlacementFromString(const std::string& a_str) noexcept
{
  if (a_str == "center") {
    return KDSplitPlacement::Center;
  }
  else if (a_str == "jitter") {
    return KDSplitPlacement::Jitter;
  }
  else if (a_str == "cell") {
    return KDSplitPlacement::Cell;
  }
  else {
    const std::string msg = "ParticleManagement::kdSplitPlacementFromString - unknown split placement '" + a_str +
                            "' (expected 'center', 'jitter' or 'cell')";
    MayDay::Abort(msg.c_str());
  }

  return KDSplitPlacement::Center;
}

inline KDAmrBoundary
kdAmrBoundaryFromString(const std::string& a_str) noexcept
{
  if (a_str == "carve") {
    return KDAmrBoundary::Carve;
  }
  else if (a_str == "nn") {
    return KDAmrBoundary::Nn;
  }
  else {
    const std::string msg = "ParticleManagement::kdAmrBoundaryFromString - unknown boundary rule '" + a_str +
                            "' (expected 'carve' or 'nn')";
    MayDay::Abort(msg.c_str());
  }

  return KDAmrBoundary::Carve;
}

inline NNSearch
nnSearchFromString(const std::string& a_str) noexcept
{
  if (a_str == "tree") {
    return NNSearch::Tree;
  }
  else if (a_str == "hash") {
    return NNSearch::Hash;
  }
  else if (a_str == "onecell") {
    return NNSearch::OneCell;
  }
  else {
    const std::string msg = "ParticleManagement::nnSearchFromString - unknown search backend '" + a_str +
                            "' (expected 'tree', 'hash' or 'onecell')";
    MayDay::Abort(msg.c_str());
  }

  return NNSearch::Tree;
}

namespace detail {

/**
 * @brief Place one piece of a split particle.
 * @details Shared by the per-cell and whole-patch builds so both spread split pieces the same way.
 *
 * The jitter kernel is TRUNCATED to the cell, not clamped and not reflected. Clamping stacks every
 * out-of-range draw onto the wall. Reflecting folds the overhanging part of the kernel back on top of the
 * part that stayed, doubling the density within half a bandwidth of each wall -- and since the fold is
 * per-direction, a parent near a corner is folded once per direction and the corner gets 2^D times the
 * density it should. Truncation is the correct conditional distribution: draw uniformly from the part of
 * the kernel that lies inside the cell. Its only cost is a mild inward bias for parents within half a
 * bandwidth of a wall, which is second order in the bandwidth and has no directional structure.
 * @param[in] a_placement Which rule to apply. Center returns the parent's position unchanged.
 * @param[in] a_parent    The parent particle's position.
 * @param[in] a_cellLo    Low corner of the cell the parent lives in.
 * @param[in] a_dx        Cell spacing, per direction.
 * @param[in] a_spacing   Jitter bandwidth -- the local mean interparticle spacing.
 * @return The piece's position, always inside the parent's cell.
 */
inline RealVect
kdSplitDisplace(const KDSplitPlacement a_placement,
                const RealVect&        a_parent,
                const RealVect&        a_cellLo,
                const RealVect&        a_dx,
                const Real             a_spacing) noexcept
{
  if (a_placement == KDSplitPlacement::Center) {
    return a_parent;
  }

  RealVect x;

  for (int dir = 0; dir < SpaceDim; dir++) {
    const Real cellLo = a_cellLo[dir];
    const Real cellHi = cellLo + a_dx[dir] * (1.0 - 1.0e-12);

    Real lo = cellLo;
    Real hi = cellHi;

    if (a_placement == KDSplitPlacement::Jitter) {
      // Truncate the kernel to the cell rather than folding it back in -- see the details above.
      lo = std::max(cellLo, a_parent[dir] - 0.5 * a_spacing);
      hi = std::min(cellHi, a_parent[dir] + 0.5 * a_spacing);
    }

    x[dir] = (hi > lo) ? (lo + (hi - lo) * Random::getUniformReal01()) : lo;
  }

  return x;
}

template <class P, Real P::*weight, RealVect P::*position>
inline void
buildKDCellLeaves(const std::vector<P>&                           a_particles,
                  const int                                       a_maxLeaves,
                  const Real                                      a_weightMedianLength,
                  const bool                                      a_capWeights,
                  const std::function<RealVect(const RealVect&)>& a_splitPosition,
                  const BinaryParticleReconcile<P>&               a_particleReconcile,
                  std::vector<std::pair<const P*, const P*>>&     a_leaves) noexcept
{
  CH_TIME("ParticleManagement::buildKDCellLeaves");

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
  // internal node, i.e. fewer than maxLeaves extra particles over the whole build. The optional weight
  // cap below adds at most maxLeaves more, since every piece it makes carries at least the cap.
  s_cur.clear();
  s_cur.reserve(numInput + 2 * maxLeaves);
  s_cur.insert(s_cur.end(), a_particles.begin(), a_particles.end());

  // Total weight, read from the working copy.
  Real W = 0.0;
  for (const P& p : s_cur) {
    W += p.*weight;
  }

  // Optional weight cap. A particle heavier than a leaf's target weight cannot be equalized away by any
  // partition -- it alone sets a floor on the leaf it lands in -- so divide it into co-located integer
  // pieces first. The division is exactly conservative and moves no mass, so unlike a median split it
  // cannot perturb sub-cell position; it only removes the floor. Bounded growth: every piece carries at
  // least one cap's worth of weight, so at most maxLeaves pieces are ever created.
  if (a_capWeights && a_maxLeaves > 1) {
    const long long cap = std::max(1LL, static_cast<long long>(std::ceil(W / static_cast<Real>(a_maxLeaves))));

    const std::size_t numOriginal = s_cur.size();

    for (std::size_t i = 0; i < numOriginal; i++) {
      const long long w = static_cast<long long>(std::llround(s_cur[i].*weight));

      if (w > cap) {
        // Integer pieces, each no heavier than the cap, summing exactly to w.
        const long long              numPieces = (w + cap - 1) / cap;
        const std::vector<long long> pieces    = partitionParticleWeights<long long>(w, numPieces);

        if (!pieces.empty()) {
          const P parent = s_cur[i];

          s_cur[i].*weight = static_cast<Real>(pieces[0]);

          for (std::size_t k = 1; k < pieces.size(); k++) {
            P daughter         = parent;
            daughter.*weight   = static_cast<Real>(pieces[k]);
            daughter.*position = a_splitPosition(parent.*position);

            // Same contract as the median split: the caller fixes up the non-weight members.
            a_particleReconcile(s_cur[i], daughter, parent);

            s_cur.push_back(std::move(daughter));
          }
        }
      }
    }
  }

  s_curNodes.clear();
  s_curNodes.push_back({static_cast<std::size_t>(0), s_cur.size(), W});

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

      // ---- split this node into s_nxt ----
      P* const          beg = s_cur.data() + node.lo;
      P* const          end = s_cur.data() + node.hi;
      const std::size_t n   = node.hi - node.lo;
      const Real        Wn  = node.w;

      // A. Split along the longest bounding-box extent. The box is also what selects the split rule, so
      // it is computed before the splittability test rather than after it -- the two questions cannot be
      // separated once the rule depends on the node's size.
      RealVect loCorner = +std::numeric_limits<Real>::max() * RealVect::Unit;
      RealVect hiCorner = -std::numeric_limits<Real>::max() * RealVect::Unit;

      for (P* p = beg; p != end; ++p) {
        const RealVect& pos = (*p).*position;

        for (int dir = 0; dir < SpaceDim; dir++) {
          loCorner[dir] = std::min(pos[dir], loCorner[dir]);
          hiCorner[dir] = std::max(pos[dir], hiCorner[dir]);
        }
      }

      const int splitDir = (hiCorner - loCorner).maxDir(true);

      // A node narrower than the crossover splits at the weight median, a wider one at the count median.
      // The caller passes a length, not a rule: 'weight' is a crossover at infinity and 'count' one below
      // zero, so a fixed rule and the size-dependent one are the same knob (see KDPartition).
      const bool weightMedian = (hiCorner[splitDir] - loCorner[splitDir]) <= a_weightMedianLength;

      // What makes a node splittable depends on the rule that was just chosen for it. The weight median
      // needs at least two units of weight to divide; the count median needs at least two particles, and
      // never divides a particle, so it can drive a node down to a single member however heavy that
      // member is.
      const bool splittable = weightMedian ? (node.w > splitThresh) : (n >= 2);

      if (splittable && splitBudget > 0) {
        splitBudget--;

        std::sort(beg, end, [splitDir](const P& p1, const P& p2) -> bool {
          return (p1.*position)[splitDir] < (p2.*position)[splitDir];
        });

        // Count median. The two halves get as nearly equal a particle COUNT as possible, and no particle
        // is ever divided, so this path creates no new particles and needs no reconciliation -- the leaf
        // count can therefore never overshoot the target during the build, unlike the weight median.
        // The weights of the halves are whatever the partition happens to give.
        if (!weightMedian) {
          const std::size_t mid = n / 2;

          Real wl = 0.0;
          for (std::size_t i = 0; i < mid; i++) {
            wl += beg[i].*weight;
          }

          const std::size_t leftLo = s_nxt.size();
          s_nxt.insert(s_nxt.end(), beg, beg + mid);
          s_nxtNodes.push_back({leftLo, s_nxt.size(), wl});

          const std::size_t rightLo = s_nxt.size();
          s_nxt.insert(s_nxt.end(), beg + mid, end);
          s_nxtNodes.push_back({rightLo, s_nxt.size(), Wn - wl});

          keepGoing = true;

          continue;
        }

        // B. Weight-balanced median particle.
        std::size_t id = 0;
        Real        wl = 0.0;
        Real        wr = Wn - beg[0].*weight;

        for (std::size_t i = 1; i < n; i++) {
          const Real& w = beg[id].*weight;
          if (wl + w < wr) {
            id = i;
            wl += w;
            wr = Wn - wl - beg[id].*weight;
          }
          else {
            break;
          }
        }

        const P    med = beg[id];
        const Real pw  = med.*weight;
        const Real dw  = wr - wl;

        CH_assert(wl + wr + pw == Wn);

        // C. Decide how the median is distributed; produce at most one extra particle per side.
        bool hasExtraL = false;
        bool hasExtraR = false;
        P    extraL;
        P    extraR;

        if (pw >= splitThresh && pw >= std::abs(dw)) {
          Real dwl = dw;
          Real dwr = 0.0;

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

            extraL.*position = a_splitPosition(med.*position);
            extraR.*position = a_splitPosition(med.*position);

            CH_assert(dwl >= 1.0);
            CH_assert(dwr >= 1.0);

            wl += dwl;
            wr += dwr;

            extraL.*weight = dwl;
            extraR.*weight = dwr;

            a_particleReconcile(extraL, extraR, med);

            hasExtraL = true;
            hasExtraR = true;
          }
          else if (dwl > 0.0 && dwr == 0.0) {
            extraL           = med;
            extraL.*position = a_splitPosition(med.*position);
            CH_assert(dwl >= 1.0);
            wl += dwl;
            extraL.*weight = dwl;
            hasExtraL      = true;
          }
          else if (dwl == 0.0 && dwr > 0.0) {
            extraR           = med;
            extraR.*position = a_splitPosition(med.*position);
            CH_assert(dwr >= 1.0);
            wr += dwr;
            extraR.*weight = dwr;
            hasExtraR      = true;
          }
          else {
            MayDay::Abort("ParticleManagement::buildKDCellLeaves - logic bust");
          }
        }
        else {
          // Median assigned whole to the lighter side (weight unchanged).
          if (wl <= wr) {
            wl        = wl + pw;
            extraL    = med;
            hasExtraL = true;
          }
          else {
            wr        = wr + pw;
            extraR    = med;
            hasExtraR = true;
          }
        }

        CH_assert(std::abs(wl + wr - Wn) <= Wn * std::numeric_limits<Real>::epsilon() * 16);
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

} // namespace detail

template <class Packed, Real Packed::*packWeight, RealVect Packed::*packPosition, class P, class Traits>
inline ParticleMerger<P, Traits>
makeKDCellMerger(
  const Real                                                                                  a_weightMedianCellWidths,
  const bool                                                                                  a_capWeights,
  const KDSplitPlacement                                                                      a_splitPlacement,
  std::function<RealVect()>                                                                   a_probLo,
  std::function<Packed(const ParticleSoA<P, Traits>&, std::size_t)>                           a_gather,
  BinaryParticleReconcile<Packed>                                                             a_reconcile,
  std::function<void(ParticleSoA<P, Traits>&, const Packed*, const Packed*, const CellInfo&)> a_scatterLeaf) noexcept
{
  return [=](ParticleSoA<P, Traits>& a_particles, const CellInfo& a_cellInfo, int a_ppc) noexcept {
    CH_TIMERS("ParticleManagement::makeKDCellMerger");
    CH_TIMER("ParticleManagement::makeKDCellMerger::populate", t1);
    CH_TIMER("ParticleManagement::makeKDCellMerger::build_kd", t2);
    CH_TIMER("ParticleManagement::makeKDCellMerger::scatter", t3);

    CH_START(t1);
    thread_local std::vector<Packed> particles;

    particles.clear();
    particles.reserve(a_particles.size());

    Real W = 0.0;

    for (std::size_t i = 0; i < a_particles.size(); i++) {
      Packed p = a_gather(a_particles, i);
      W += p.*packWeight;
      particles.emplace_back(std::move(p));
    }
    CH_STOP(t1);

    if (W < 2.0 || a_ppc <= 0) {
      return;
    }

    // Where a split particle's pieces go. Center is the identity, so the default costs nothing. In a cut
    // cell every rule falls back to Center: a displaced piece could otherwise land inside the boundary.
    const Real     dxCell  = a_cellInfo.getDx();
    const RealVect cellLo  = a_probLo() + dxCell * RealVect(a_cellInfo.getGridIndex());
    const bool     regular = a_cellInfo.getVolFrac() >= 1.0;

    // Local mean interparticle spacing -- the finest scale this cell's population can resolve, and so the
    // natural jitter bandwidth. No free parameter.
    const Real spacing = dxCell / std::pow(static_cast<Real>(std::max(1, a_ppc)), 1.0 / SpaceDim);

    const std::function<RealVect(const RealVect&)> splitPosition = [=](const RealVect& a_parent) -> RealVect {
      if (!regular || a_splitPlacement == KDSplitPlacement::Center) {
        return a_parent;
      }

      return detail::kdSplitDisplace(a_splitPlacement, a_parent, cellLo, dxCell * RealVect::Unit, spacing);
    };

    CH_START(t2);
    thread_local std::vector<std::pair<const Packed*, const Packed*>> leaves;

    // The crossover is given in cell widths, so it becomes a length only once the cell is known.
    detail::buildKDCellLeaves<Packed, packWeight, packPosition>(particles,
                                                                a_ppc,
                                                                a_weightMedianCellWidths * a_cellInfo.getDx(),
                                                                a_capWeights,
                                                                splitPosition,
                                                                a_reconcile,
                                                                leaves);
    CH_STOP(t2);

    CH_START(t3);
    a_particles.clear();
    for (const auto& leaf : leaves) {
      a_scatterLeaf(a_particles, leaf.first, leaf.second, a_cellInfo);
    }
    CH_STOP(t3);
  };
}

template <class Context, class P, class Traits>
inline ParticleMerger<P, Traits>
makeReinitializeMerger(std::function<std::pair<long long, Context>(const ParticleSoA<P, Traits>&)> a_aggregate,
                       std::function<void(ParticleSoA<P, Traits>&, const RealVect&, long long, const Context&)> a_emit,
                       std::function<RealVect()> a_probLo) noexcept
{
  return [=](ParticleSoA<P, Traits>& a_particles, const CellInfo& a_cellInfo, int a_ppc) noexcept {
    CH_TIME("ParticleManagement::makeReinitializeMerger");

    if (a_ppc <= 0) {
      return;
    }

    const auto [numPhysical, context] = a_aggregate(a_particles);

    if (numPhysical <= 0LL) {
      return;
    }

    const std::vector<long long> weights = partitionParticleWeights(numPhysical, (long long)a_ppc);

    const Real      dx            = a_cellInfo.getDx();
    const Real      kappa         = a_cellInfo.getVolFrac();
    const RealVect  cellPos       = a_probLo() + dx * (a_cellInfo.getGridIndex() + 0.5 * RealVect::Unit);
    const RealVect& validLo       = a_cellInfo.getValidLo();
    const RealVect& validHi       = a_cellInfo.getValidHi();
    const RealVect& bndryCentroid = a_cellInfo.getBndryCentroid();
    const RealVect& bndryNormal   = a_cellInfo.getBndryNormal();

    a_particles.clear();

    for (const long long wt : weights) {
      const RealVect x = Random::randomPosition(cellPos, validLo, validHi, bndryCentroid, bndryNormal, dx, kappa);
      a_emit(a_particles, x, wt, context);
    }
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

        for (std::size_t i = 0; i < ret.size() && r > zero; i++) {
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

namespace detail {

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

} // namespace detail

template <typename P, typename Traits, typename T, typename>
inline void
drawRandomParticles(ParticleSoA<P, Traits>&          a_particles,
                    const T                          a_numParticles,
                    const std::function<RealVect()>& a_distribution)
{
  a_particles.clear();

  // Each rank draws its own share (unit weight) into the buffer. Routing to owners happens later when
  // the buffer is added to a ParticleContainer.
  const T numParticles = detail::partitionParticles(a_numParticles);

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
