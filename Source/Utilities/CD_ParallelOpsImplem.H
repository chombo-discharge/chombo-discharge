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

inline void
ParallelOps::barrier()
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
}

inline std::pair<int, int>
ParallelOps::partition(const int a_N) noexcept
{
#ifdef CH_MPI
  const int equalChunk = a_N / numProc();
  const int remainder  = a_N % numProc();

  int begin = procID() * equalChunk;
  int end   = begin + equalChunk - 1;

  // Remainder goes to last rank.
  if (procID() == numProc() - 1) {
    end += remainder;
  }

  return std::make_pair(begin, end);
#else
  return std::make_pair(0, N - 1);
#endif
}

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

inline std::pair<Real, RealVect>
ParallelOps::min(const Real& a_val, const RealVect& a_pos) noexcept
{
#ifdef CH_MPI
  const int nProc = numProc();

  Real  d = a_val;
  Real* D = new Real[nProc];

  Real  x = a_pos[0];
  Real* X = new Real[nProc];

  Real  y = a_pos[0];
  Real* Y = new Real[nProc];

#if CH_SPACEDIM == 3
  Real  z = a_pos[2];
  Real* Z = new Real[nProc];
#endif

  MPI_Allgather(&d, 1, MPI_CH_REAL, D, 1, MPI_CH_REAL, Chombo_MPI::comm);
  MPI_Allgather(&x, 1, MPI_CH_REAL, X, 1, MPI_CH_REAL, Chombo_MPI::comm);
  MPI_Allgather(&y, 1, MPI_CH_REAL, Y, 1, MPI_CH_REAL, Chombo_MPI::comm);
#if CH_SPACEDIM == 3
  MPI_Allgather(&z, 1, MPI_CH_REAL, Z, 1, MPI_CH_REAL, Chombo_MPI::comm);
#endif

  Real dRet = D[0];
  Real xRet = X[0];
  Real yRet = Y[0];
#if CH_SPACEDIM == 3
  Real zRet = Z[0];
#endif

  for (int i = 1; i < nProc; i++) {
    if (D[i] < dRet) {
      dRet = D[i];
      xRet = X[i];
      yRet = Y[i];
#if CH_SPACEDIM == 3
      zRet = Z[i];
#endif
    }
  }

  delete[] D;
  delete[] X;
  delete[] Y;
#if CH_SPACEDIM == 3
  delete[] Z;
#endif

  return std::make_pair(dRet, RealVect(D_DECL(xRet, yRet, zRet)));
#else
  return std::make_pair(a_val, a_pos);
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

inline std::pair<Real, RealVect>
ParallelOps::max(const Real& a_val, const RealVect& a_pos) noexcept
{
#ifdef CH_MPI
  const int nProc = numProc();

  Real  d = a_val;
  Real* D = new Real[nProc];

  Real  x = a_pos[0];
  Real* X = new Real[nProc];

  Real  y = a_pos[0];
  Real* Y = new Real[nProc];

#if CH_SPACEDIM == 3
  Real  z = a_pos[2];
  Real* Z = new Real[nProc];
#endif

  MPI_Allgather(&d, 1, MPI_CH_REAL, D, 1, MPI_CH_REAL, Chombo_MPI::comm);
  MPI_Allgather(&x, 1, MPI_CH_REAL, X, 1, MPI_CH_REAL, Chombo_MPI::comm);
  MPI_Allgather(&y, 1, MPI_CH_REAL, Y, 1, MPI_CH_REAL, Chombo_MPI::comm);
#if CH_SPACEDIM == 3
  MPI_Allgather(&z, 1, MPI_CH_REAL, Z, 1, MPI_CH_REAL, Chombo_MPI::comm);
#endif

  Real dRet = D[0];
  Real xRet = X[0];
  Real yRet = Y[0];
#if CH_SPACEDIM == 3
  Real zRet = Z[0];
#endif

  for (int i = 1; i < nProc; i++) {
    if (D[i] > dRet) {
      dRet = D[i];
      xRet = X[i];
      yRet = Y[i];
#if CH_SPACEDIM == 3
      zRet = Z[i];
#endif
    }
  }

  delete[] D;
  delete[] X;
  delete[] Y;
#if CH_SPACEDIM == 3
  delete[] Z;
#endif

  return std::make_pair(dRet, RealVect(D_DECL(xRet, yRet, zRet)));
#else
  return std::make_pair(a_val, a_pos);
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
ParallelOps::vectorSum(Vector<Real>& a_data)
{
#ifdef CH_MPI
  const Vector<Real> tmp = a_data;
  const int result = MPI_Allreduce(&(tmp[0]), &(a_data[0]), a_data.size(), MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("In file ParallelOps::vectorSum -- MPI communication error");
  }
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
