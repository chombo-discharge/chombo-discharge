/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file  CD_DataOpsImplem.H
  @brief Implementation of CD_DataOps.H
*/

#ifndef CD_DataOpsImplem_H
#define CD_DataOpsImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_NamespaceHeader.H>

template <typename Functor>
ALWAYS_INLINE
void DataOps::loop(const Box& computeBox, Functor&& func){

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

template <typename Functor>
ALWAYS_INLINE
void DataOps::loop(VoFIterator& vofit, Functor&& func){

  // TLDR: This runs through all cells in the vof-iterator and calls the kernel. 
  for (vofit.reset(); vofit.ok(); ++vofit){
    func(vofit());
  }
}    

template <typename T> void DataOps::setValue(LevelData<MFInterfaceFAB<T> >& a_lhs, const T& a_value){
  CH_TIME("DataOps::setValue");
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()].setVal(a_value);
  }
}

template <typename T> int DataOps::sgn(const T a_value){
  return ((a_value > 0.0) - (a_value < 0.0));
}

#include <CD_NamespaceFooter.H>

#endif
