/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file    CD_LoadBalancingImplem.H
  @details Implementation of CD_LoadBalancing.H
  @author  Robert Marskar
*/

#ifndef CD_LoadBalanceImplem_H
#define CD_LoadBalanceImplem_H

// Std includes
#include <algorithm>
#include <random>
#include <chrono>

// Chombo includes
#include <LoadBalance.H>
#include <BoxLayout.H>
#include <CH_Timer.H>

// Our includes
#include <CD_LoadBalancing.H>
#include <CD_NamespaceHeader.H>

template <class T>
void
LoadBalancing::makeBalance(Vector<int>& a_ranks, const Vector<T>& a_loads, const Vector<Box>& a_boxes)
{
  CH_TIME("LoadBalancing::makeBalance");

  //  LoadBalance(a_ranks, a_loads, a_boxes);
  Loads rankLoads;
  LoadBalancing::makeBalance<T>(a_ranks, rankLoads, a_loads, a_boxes);
}

template <class T>
void
LoadBalancing::makeBalance(Vector<int>&       a_ranks,
                           Loads&             a_rankLoads,
                           const Vector<T>&   a_boxLoads,
                           const Vector<Box>& a_boxes)
{
  CH_TIME("LoadBalancing::makeBalance");

  // Convert everything to floating points
  Vector<Real> boxLoads;
  for (int i = 0; i < a_boxLoads.size(); i++) {
    boxLoads.push_back(1.0 * a_boxLoads[i]);
  }

  // Minimum number of grid subsets is the number of boxe, and the maximum number of grid subsets
  // is the number of ranks.
  const int numBoxes   = a_boxes.size();
  const int numRanks   = numProc();
  const int numSubsets = std::min(numBoxes, numRanks);

  a_ranks.resize(numBoxes);

  if (numSubsets > 0) {

    // Figure out the total and target load (load per subset) on this level.
    Real totalLoad = 0.0;
    for (int ibox = 0; ibox < numBoxes; ibox++) {
      totalLoad += boxLoads[ibox];
    }

    Real staticTargetLoad = totalLoad / numSubsets;

    // Build the grid subsets. When we do this we iterate through the boxes and try to ensure that we partition
    // the subsets such that the subsetLoad is as close to the dynamic targetLoad as possible.
    //
    // The pair contains the starting index for the subset and the computational load for the subset.
    using Span   = std::pair<int, int>;
    using Subset = std::pair<Span, Real>;

    std::vector<Subset> subsets(numSubsets);

    int firstSubsetBox = 0;

    Real remainingLoad = totalLoad;

    for (int curSubset = 0; curSubset < numSubsets; curSubset++) {

      // The firstSubsetBox is the index for the first box in this subset (always assigned).
      Real subsetLoad = boxLoads[firstSubsetBox];

      int lastSubsetBox = firstSubsetBox;

      const int subsetsLeft = numSubsets - (curSubset + 1);
      const int boxesLeft   = numBoxes - (firstSubsetBox + 1);

      if (boxesLeft > subsetsLeft) {
        for (int ibox = firstSubsetBox + 1; ibox < numBoxes; ibox++) {

          // Hook for catching case when we add too many boxes to this subset. Each remaining subset must have at least one box.
          if (numBoxes - lastSubsetBox - 1 <= subsetsLeft) {
            break;
          }

          // Check if we should add this box - we do this by making sure that the dynamically moving target load stays as close
          // to the static load as possible.
          //
          // In the below, '1' is the load without ibox, and '2' is the load with ibox
          const Real load1 = subsetLoad;
          const Real load2 = subsetLoad + boxLoads[ibox];

          // Check if we should add this box.
          bool addBoxToSubset = false;

          if (boxLoads[ibox] <= std::numeric_limits<Real>::epsilon()) {
            addBoxToSubset = true;
          }
          else if (load1 > staticTargetLoad) {
            addBoxToSubset = false;
          }
          else if (load2 <= staticTargetLoad) {
            addBoxToSubset = true;
          }
          else if (load1 <= staticTargetLoad && load2 > staticTargetLoad) {
            // Compute the new average load if we add or don't add this box to the current subset. Accept the answer
            // that leads to a smallest deviation from the static target load.
            const Real loadErrWithoutBox = std::abs(load1 - staticTargetLoad);
            const Real loadErrWithBox    = std::abs(load2 - staticTargetLoad);

            if (loadErrWithBox <= loadErrWithoutBox) {
              addBoxToSubset = true;
            }
          }

          // Add box or break out of box iteration.
          if (addBoxToSubset) {
            subsetLoad    = subsetLoad + boxLoads[ibox];
            lastSubsetBox = ibox;

            continue;
          }
          else {
            lastSubsetBox = ibox - 1;

            break;
          }
        }
      }

      // Create the subset
      subsets[curSubset] = std::make_pair(std::make_pair(firstSubsetBox, lastSubsetBox), subsetLoad);

      // Update the remaining load.
      remainingLoad    = remainingLoad - subsetLoad;
      staticTargetLoad = remainingLoad / subsetsLeft;

      // Update start box for next iteration.
      firstSubsetBox = lastSubsetBox + 1;
    }

    // Sort the subsets from largest to smallest computational load.
    std::sort(subsets.begin(), subsets.end(), [](const Subset& A, const Subset& B) -> bool {
      return A.second > B.second;
    });

    // Get the accumulated loads per rank, sorted by lowest-to-highest load.
    const std::vector<std::pair<int, Real>> sortedRankLoads = a_rankLoads.getSortedLoads();

    // Assign the most expensive grid subset to the rank with the lowest accumulated load.
    for (int i = 0; i < subsets.size(); i++) {
      const int  startIndex = subsets[i].first.first;
      const int  endIndex   = subsets[i].first.second;
      const Real subsetLoad = subsets[i].second;
      const int  rank       = sortedRankLoads[i].first;

      // Assign to rank
      for (int ibox = startIndex; ibox <= endIndex; ibox++) {
        a_ranks[ibox] = rank;
      }

      // Update load on this rank.
      a_rankLoads.incrementLoad(rank, subsetLoad);
    }
  }
  else {
    a_ranks.resize(0);
  }
}

template <class T>
std::vector<std::pair<Box, T>>
LoadBalancing::packPairs(const Vector<Box>& a_boxes, const Vector<T>& a_loads)
{
  CH_TIME("LoadBalancing::packPairs");

  std::vector<std::pair<Box, T>> vec;
  for (int i = 0; i < a_boxes.size(); i++) {
    vec.emplace_back(a_boxes[i], a_loads[i]);
  }

  return vec;
}

template <class T>
void
LoadBalancing::unpackPairs(Vector<Box>& a_boxes, Vector<T>& a_loads, const std::vector<std::pair<Box, T>>& a_pairs)
{
  CH_TIME("LoadBalancing::unpackPairs");

  // Reconstruct boxes and loads
  a_boxes.resize(0);
  a_loads.resize(0);

  for (const auto& v : a_pairs) {
    a_boxes.push_back(v.first);
    a_loads.push_back(v.second);
  }
}

template <typename T>
void
LoadBalancing::sort(Vector<Vector<Box>>& a_boxes, Vector<Vector<T>>& a_loads, const BoxSorting a_which)
{
  CH_TIME("LoadBalancing::sort");

  for (int lvl = 0; lvl < a_boxes.size(); lvl++) {
    LoadBalancing::sort(a_boxes[lvl], a_loads[lvl], a_which);
  }
}

template <typename T>
void
LoadBalancing::sort(Vector<Box>& a_boxes, Vector<T>& a_loads, const BoxSorting a_which)
{
  CH_TIME("LoadBalancing::sort");

  switch (a_which) {
  case BoxSorting::None: {
    break;
  }
  case BoxSorting::Std: {
    LoadBalancing::standardSort(a_boxes, a_loads);

    break;
  }
  case BoxSorting::Shuffle: {
    LoadBalancing::shuffleSort(a_boxes, a_loads);

    break;
  }
  case BoxSorting::Morton: {
    LoadBalancing::mortonSort(a_boxes, a_loads);

    break;
  }
  default: {
    MayDay::Abort("LoadBalancing::sort_boxes - unknown algorithm requested");

    break;
  }
  }
}

template <class T>
void
LoadBalancing::standardSort(Vector<Box>& a_boxes, Vector<T>& a_loads)
{
  CH_TIME("LoadBalancing::standardSort");

  std::vector<std::pair<Box, T>> vec = packPairs(a_boxes, a_loads);

  // Call std::sort, using box1 < box2 lambda as sorting criterion.
  std::sort(std::begin(vec), std::end(vec), [](const std::pair<Box, T>& v1, const std::pair<Box, T>& v2) {
    return v1.first < v2.first;
  });

  unpackPairs(a_boxes, a_loads, vec);
}

template <class T>
void
LoadBalancing::shuffleSort(Vector<Box>& a_boxes, Vector<T>& a_loads)
{
  CH_TIME("LoadBalancing::shuffleSort");

  auto vec = packPairs(a_boxes, a_loads);

  // Set up RNG
  int seed = std::chrono::system_clock::now().time_since_epoch().count();
#ifdef CH_MPI // Broadcast
  MPI_Bcast(&seed, 1, MPI_INT, 0, Chombo_MPI::comm);
#endif

  // Shuffle vector
  std::default_random_engine e(seed);
  std::shuffle(vec.begin(), vec.end(), e);

  // Split boxes and loads

  unpackPairs(a_boxes, a_loads, vec);
}

template <class T>
void
LoadBalancing::mortonSort(Vector<Box>& a_boxes, Vector<T>& a_loads)
{
  CH_TIME("LoadBalancing::mortonSort");

  auto vec = packPairs(a_boxes, a_loads);

  // Get max bits
  std::vector<Box>& b    = a_boxes.stdVector();
  int               bits = maxBits(b.begin(), b.end());

  // Morton sort.
  std::sort(std::begin(vec), std::end(vec), [bits](const std::pair<Box, T>& v1, const std::pair<Box, T>& v2) -> bool {
    return mortonComparator(bits, v1, v2);
  });

  // Put back in normal form
  unpackPairs(a_boxes, a_loads, vec);
}

template <class T>
bool
LoadBalancing::mortonComparator(const int a_maxBits, const std::pair<Box, T>& a_lhs, const std::pair<Box, T>& a_rhs)
{
  const Box& lbox = a_lhs.first;
  const Box& rbox = a_rhs.first;

  const IntVect l = lbox.smallEnd();
  const IntVect r = rbox.smallEnd();

  for (int i = a_maxBits; i > 0; i--) {

    // March from most significant bit to least.
    const int N = (1 << i);

    for (int dir = CH_SPACEDIM - 1; dir >= 0; dir--) {
      if ((l[dir] / N) < (r[dir] / N)) {
        return true;
      }
      else if ((l[dir] / N) > (r[dir] / N)) {
        return false;
      }
    }
  }

  return false;
}

#include <CD_NamespaceFooter.H>

#endif
