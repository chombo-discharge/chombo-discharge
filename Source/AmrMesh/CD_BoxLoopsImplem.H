/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file  CD_BoxLoopsImplem.H
  @brief Implementation of CD_BoxLoops.H
*/

#ifndef CD_BoxLoopsImplem_H
#define CD_BoxLoopsImplem_H

// Our includes
#include <CD_NamespaceHeader.H>

namespace BoxLoops {
  
  template <typename Functor>
  ALWAYS_INLINE
  void loop(const Box& computeBox, Functor&& func){

    // TLDR: This runs through all cells in the computeBox and calls the kernel function.
  
    const int* lo = computeBox.loVect();
    const int* hi = computeBox.hiVect();

#if CH_SPACEDIM==3
    for (int k = lo[2]; k <= hi[2]; k++){
#endif
      for (int j = lo[1]; j <= hi[1]; j++){
	CD_PRAGMA_SIMD
	  for (int i = lo[0]; i <= hi[0]; i++){
	    func(IntVect(D_DECL(i,j,k)));
	  }
      }
#if CH_SPACEDIM==3      
    }
#endif
  }

  template <typename Iterator, typename Functor>
  ALWAYS_INLINE
  void loop(Iterator& iter, Functor&& func){

    // TLDR: This runs through all cells in the vof-iterator and calls the kernel. 
    for (iter.reset(); iter.ok(); ++iter){
      func(iter());
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
