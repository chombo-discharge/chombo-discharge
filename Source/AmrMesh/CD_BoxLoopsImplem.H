/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BoxLoopsImplem.H
  @brief  Implementation of CD_BoxLoops.H
  @author Robert Marskar
*/

#ifndef CD_BoxLoopsImplem_H
#define CD_BoxLoopsImplem_H

// Our includes
#include <CD_NamespaceHeader.H>

template <typename Functor>
ALWAYS_INLINE void
BoxLoops::loop(const Box& computeBox, Functor&& func, const IntVect& a_stride)
{

  CH_assert(a_stride > IntVect::Zero);

  const int* lo = computeBox.loVect();
  const int* hi = computeBox.hiVect();

  // TLDR: This runs through all cells in the computeBox and calls the kernel function.
#if CH_SPACEDIM == 3
  for (int k = lo[2]; k <= hi[2]; k += a_stride[2]) {
#endif
    for (int j = lo[1]; j <= hi[1]; j += a_stride[1]) {
      CD_PRAGMA_SIMD
      for (int i = lo[0]; i <= hi[0]; i += a_stride[0]) {
        func(IntVect(D_DECL(i, j, k)));
      }
    }
#if CH_SPACEDIM == 3
  }
#endif
}

template <typename Functor>
ALWAYS_INLINE void
BoxLoops::loop(const IntVectSet& a_ivs, Functor&& a_kernel)
{
  for (IVSIterator iter(a_ivs); iter.ok(); ++iter) {
    a_kernel(iter());
  }
}

template <typename Functor>
ALWAYS_INLINE void
BoxLoops::loop(const DenseIntVectSet& a_ivs, Functor&& a_kernel)
{
  for (DenseIntVectSetIterator iter(a_ivs); iter.ok(); ++iter) {
    a_kernel(iter());
  }
}

template <typename Functor>
ALWAYS_INLINE void
BoxLoops::loop(VoFIterator& iter, Functor&& a_kernel)
{

  // TLDR: This runs through all cells in the vof-iterator and calls the kernel.
  for (iter.reset(); iter.ok(); ++iter) {
    a_kernel(iter());
  }
}

template <typename Functor>
ALWAYS_INLINE void
BoxLoops::loop(FaceIterator& iter, Functor&& a_kernel)
{

  // TLDR: This runs through all cells in the vof-iterator and calls the kernel.
  for (iter.reset(); iter.ok(); ++iter) {
    a_kernel(iter());
  }
}

template <typename T, typename Functor>
ALWAYS_INLINE void
BoxLoops::loop(const Vector<T>& a_subset, Functor&& a_kernel)
{
  const std::vector<T>& stdVec = ((Vector<T>&)a_subset).stdVector();

  for (const auto& v : stdVec) {
    a_kernel(v);
  }
}

#include <CD_NamespaceFooter.H>

#endif
