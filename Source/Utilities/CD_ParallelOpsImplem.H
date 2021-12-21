/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ParallelOpsImplem.H
  @brief  Implementation of CD_ParallelOps.H
  @author Robert Marskar
*/

#ifndef CD_ParallelOpsImplem_H
#define CD_ParallelOpsImplem_H

// Chombo includes
#include <SPMD.H>

// Our includes
#include <CD_NamespaceHeader.H>

inline  
Real ParallelOps::Min(const Real& a_input) {
#ifdef CH_MPI
  Real ret = 1.0;
  const int result = MPI_Allreduce(&a_input, &ret, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("In file ParallelOps::Min -- MPI communication error");
  }
  return ret;
#else
  return a_input;
#endif
}

inline
Real ParallelOps::Max(const Real& a_input) {

#ifdef CH_MPI
  Real ret = 1.0;
  const int result = MPI_Allreduce(&a_input, &ret, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("In file ParallelOps::Max -- MPI communication error");
  }
  return ret;
#else
  return a_input;
#endif
}

inline  
int ParallelOps::Min(const int& a_input) {
#ifdef CH_MPI
  int ret = 1.0;
  const int result = MPI_Allreduce(&a_input, &ret, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("In file ParallelOps::Min -- MPI communication error");
  }
  return ret;
#else
  return a_input;
#endif
}

inline
int ParallelOps::Max(const int& a_input) {

#ifdef CH_MPI
  int ret = 1.0;
  const int result = MPI_Allreduce(&a_input, &ret, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("In file ParallelOps::Max -- MPI communication error");
  }
  return ret;
#else
  return a_input;
#endif
}

inline
void ParallelOps::VectorSum(Vector<long int>& a_data) {
#ifdef CH_MPI
  const Vector<long int> tmp = a_data;
  const int result = MPI_Allreduce(&(tmp[0]),&(a_data[0]), a_data.size(), MPI_LONG, MPI_SUM, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("In file ParallelOps::VectorSum -- MPI communication error");
  }
#endif
}

#include <CD_NamespaceFooter.H>

#endif
