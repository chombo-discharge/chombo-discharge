/*!
  @file    load_balanceI.H
  @details Implementation of load_balance.H
  @author  Robert Marskar
  @date    Jan. 2018
*/

#ifndef _LOADBALANCEI_H_
#define _LOADBALANCEI_H_

#include "load_balance.H"

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
  
  Vector<Box> allBoxes(0);
  Vector<T>   allLoads(0);
  Vector<int> allProcs(0);

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

#endif
