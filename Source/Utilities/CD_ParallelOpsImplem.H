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

inline Real
ParallelOps::min(const Real& a_input)
{
#ifdef CH_MPI
  Real      ret    = 1.0;
  const int result = MPI_Allreduce(&a_input, &ret, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("In file ParallelOps::min -- MPI communication error");
  }
  return ret;
#else
  return a_input;
#endif
}

inline Real
ParallelOps::max(const Real& a_input)
{

#ifdef CH_MPI
  Real      ret    = 1.0;
  const int result = MPI_Allreduce(&a_input, &ret, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("In file ParallelOps::max -- MPI communication error");
  }
  return ret;
#else
  return a_input;
#endif
}

inline int
ParallelOps::min(const int& a_input)
{
#ifdef CH_MPI
  int       ret    = 1;
  const int result = MPI_Allreduce(&a_input, &ret, 1, MPI_INT, MPI_MIN, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("In file ParallelOps::min -- MPI communication error");
  }
  return ret;
#else
  return a_input;
#endif
}

inline int
ParallelOps::max(const int& a_input)
{

#ifdef CH_MPI
  int       ret    = 1;
  const int result = MPI_Allreduce(&a_input, &ret, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("In file ParallelOps::max -- MPI communication error");
  }
  return ret;
#else
  return a_input;
#endif
}

inline Real
ParallelOps::sum(const Real& a_value)
{
#ifdef CH_MPI
  Real ret = 0.0;

  const int result = MPI_Allreduce(&a_value, &ret, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("In file ParallelOps::sum -- MPI communication error");
  }

  return ret;
#else
  return a_value;
#endif
}

inline void
ParallelOps::vectorSum(Vector<long int>& a_data)
{
#ifdef CH_MPI
  const Vector<long int> tmp = a_data;
  const int result = MPI_Allreduce(&(tmp[0]), &(a_data[0]), a_data.size(), MPI_LONG, MPI_SUM, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("In file ParallelOps::vectorSum -- MPI communication error");
  }
#endif
}

#include <CD_NamespaceFooter.H>

#endif
