/*!
  @file    load_balanceI.H
  @details Implementation of load_balance.H
  @author  Robert Marskar
  @date    Jan. 2018
*/

#ifndef _LOADBALANCEI_H_
#define _LOADBALANCEI_H_

#include "load_balance.H"

#include <algorithm>
#include <random>
#include <chrono>

template <class T>
void load_balance::make_balance(Vector<int>& a_levelRanks, const Vector<T>& a_loads, const Vector<Box>& a_levelBoxes) {
  LoadBalance(a_levelRanks, a_loads, a_levelBoxes);
}

template <class T>
void load_balance::level_by_level(Vector<Vector<int> >&       a_procs,
				  const Vector<Vector<T> >&   a_loads,
				  const Vector<Vector<Box> >& a_boxes){

  const int nlevels = a_boxes.size();

  for (int lvl = 0; lvl < nlevels; lvl++){
    load_balance::make_balance(a_procs[lvl], a_loads[lvl], a_boxes[lvl]);
  }
}

template <class T>
void load_balance::hierarchy(Vector<Vector<int> >&       a_procs,
			     const Vector<Vector<T> >&   a_loads,
			     const Vector<Vector<Box> >& a_boxes){

  const int nlevels = a_procs.size();
  
  std::map<std::pair<int, int>, int> boxmap;
  
  Vector<Box> allBoxes;
  Vector<T>   allLoads;
  Vector<int> allProcs;

  int gbox = 0;

  // Make the hierarchy to level mapping
  for (int lvl = 0; lvl < nlevels; lvl++){
    for (int ibox = 0; ibox < a_boxes[lvl].size(); ibox++){

      allBoxes.push_back(a_boxes[lvl][ibox]);
      allLoads.push_back(a_loads[lvl][ibox]);

      boxmap.emplace(std::pair<int, int>(lvl, ibox), gbox);

      gbox++;
    }
  }

  // Do global load balancing
  load_balance::make_balance(allProcs, allLoads, allBoxes);

  // Resize the a_procs vector to correct size.
  a_procs.resize(nlevels);
  for (int lvl = 0; lvl < nlevels; lvl++){
    a_procs[lvl].resize(a_boxes[lvl].size(), -1);
  }

  // Do the global to local mapping
  for (const auto& bmap : boxmap){
    const int gid = bmap.second;        // Global index
    const int lvl = bmap.first.first;   // Current level
    const int lid = bmap.first.second;  // Local box index on level lvl
    const int pid = allProcs[gid];      // Processor ownership

    a_procs[lvl][lid] = pid;
  }
  
}

template <class T>
std::vector<std::pair<Box, T> > load_balance::makeVectorPair(Vector<Box>& a_boxes, Vector<T>& a_loads){
  
  std::vector<std::pair<Box, T> > vec;
  for (int i = 0; i < a_boxes.size(); i++){
    vec.emplace_back(a_boxes[i], a_loads[i]);
  }

  return vec;
}

template <class T>
void load_balance::splitVectorPair(Vector<Box>& a_boxes, Vector<T>& a_loads, const std::vector<std::pair<Box, T> >& a_pairs){

  // Reconstruct boxes and loads
  a_boxes.resize(0);
  a_loads.resize(0);
  for (const auto& v : a_pairs){
    a_boxes.push_back(v.first);
    a_loads.push_back(v.second);
  }
}

template <typename T>
void load_balance::sort(Vector<Box>& a_boxes, Vector<T>& a_loads, const box_sorting a_which){

  switch(a_which){
  case box_sorting::none:
    break;
  case box_sorting::std:
    load_balance::std_sort(a_boxes, a_loads);
    break;
  case box_sorting::shuffle:
    load_balance::shuffle_sort(a_boxes, a_loads);
    break;
  case box_sorting::morton:
    load_balance::morton_sort(a_boxes, a_loads);
    break;
  default:
    MayDay::Abort("load_balance::sort_boxes - unknown algorithm requested");
    break;
  }
}

template <class T>
void load_balance::std_sort(Vector<Box>& a_boxes, Vector<T>& a_loads){
  std::vector<std::pair<Box, T> > vec = makeVectorPair(a_boxes, a_loads);

  // Call std::sort, using box1 < box2 lambda as sorting criterion. 
  std::sort(std::begin(vec),
	    std::end(vec),
	    [](const std::pair<Box, T>& v1,
	       const std::pair<Box, T>& v2){

	      return v1.first < v2.first;
	    });

  splitVectorPair(a_boxes, a_loads, vec);
}

template <class T>
void load_balance::shuffle_sort(Vector<Box>& a_boxes, Vector<T>& a_loads){
  auto vec = makeVectorPair(a_boxes, a_loads);

  // Set up RNG
  int seed = std::chrono::system_clock::now().time_since_epoch().count();
#ifdef CH_MPI // Broadcast
  MPI_Bcast(&seed, 1, MPI_INT, 0, Chombo_MPI::comm);
#endif

  // Shuffle vector
  std::default_random_engine e(seed);
  std::shuffle(vec.begin(), vec.end(), e);

  // Split boxes and loads

  splitVectorPair(a_boxes, a_loads, vec);
}

template <class T>
void load_balance::morton_sort(Vector<Box>& a_boxes, Vector<T>& a_loads){
  MayDay::Abort("load_balance::morton_sort(boxes, loads) - not implemented");
}

#endif
