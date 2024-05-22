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

  LoadBalance(a_ranks, a_loads, a_boxes);
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
