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

  // Minimum number of grid subsets is the number of boxe, and the maximum number of grid subsets
  // is the number of ranks.
  const int numBoxes   = a_boxes.size();
  const int numRanks   = numProc();
  const int numSubsets = std::min(numBoxes, numRanks);

  if (numSubsets > 0) {

    // Figure out the total and target load (load per subset) on this level.
    Real totalLoad  = 0;
    Real targetLoad = 0;

    for (int ibox = 0; ibox < numBoxes; ibox++) {
      totalLoad += 1.0 * a_boxLoads[ibox];
      targetLoad += 1.0 * a_boxLoads[ibox] / numSubsets;
    }

    // Build the grid subsets. When we do this we iterate through the boxes and try to ensure that we partition
    // the subsets such that the subsetLoad is as close to the targetLoad as possible, and that the total load
    // for the subsets we've calculated so far is as close as possible to to the targetLoad * number of subsets so far.
    //
    // The pair contains the starting index for the subset and the computational load for the subset.
    using Span   = std::pair<int, int>;
    using Subset = std::pair<Span, Real>;

    std::vector<Subset> subsets;

    // If we have fewer boxes than ranks, each subset is just one of the box and one load.
    if (numBoxes <= numRanks) {
      for (int ibox = 0; ibox < numBoxes; ibox++) {
        subsets.emplace_back(std::make_pair(std::make_pair(ibox, ibox), 1.0 * a_boxLoads[ibox]));
      }
    }
    else {
      // If we have more boxes than ranks, this gets a bit more complicated. But still doable.
      int startIndex = 0;
      int endIndex   = -1;

      Real curLoadSoFar = 0.0;
      Real avgLoadSoFar = 0.0;
      int  curSubset    = 1;

      for (int ibox = 1; ibox < numBoxes; ibox++) {
        const Real curLoadWithoutThisBox = curLoadSoFar;
        const Real curLoadWithThisBox    = curLoadSoFar + a_boxLoads[ibox];

        // If we exceed the target load, we must include either ibox-1 or ibox
        if (curLoadWithThisBox > targetLoad * curSubset) {
        }
      }
    }

    // Sort the subsets from largest to smallest computational load.
    std::sort(subsets.begin(), subsets.end(), [](const Subset& A, const Subset& B) -> bool {
      return A.second >= B.second;
    });

    // Get the accumulated loads per rank, sorted by lowest-to-highest load.
    const std::vector<std::pair<int, Real>>& sortedRankLoads = a_rankLoads.getSortedLoads();

    // Assign the most expensive grid subset to the rank with the lowest accumulated load - then update the loads.
  }
  else {
    a_ranks.resize(0);
  }

#if 1 // Chombo fallback code - to be removed later.
  LoadBalance(a_ranks, a_boxLoads, a_boxes);
#endif
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
