/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BoxLoops.H
  @brief  Declaration of a namespace for proto-typing grid and EB loops.
  @author Robert Marskar
*/

#ifndef CD_BoxLoops_H
#define CD_BoxLoops_H

// Chombo includes
#include <Box.H>
#include <IntVectSet.H>
#include <DenseIntVectSet.H>
#include <VoFIterator.H>
#include <FaceIterator.H>
#include <Vector.H>

// Our includes
#include <CD_Decorations.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Namespace for encapsulating loop types.
  @details BoxLoops provides several overloads of @c loop for iterating over structured and
  unstructured grid data. The general pattern is:

  @code{.cpp}
  BoxLoops::loop(domain, kernel);
  @endcode

  where @c kernel is a callable (lambda or functor) taking the element type for that overload.

  <b>Box loop — compile-time stride:</b>
  Strides are non-type template parameters, making their value visible to the compiler so the
  inner loop can be auto-vectorized. All strides must be >= 1 (enforced by @c static_assert).
  Use @c D_DECL for dimension-independent call sites. Unit-stride iteration is @c D_DECL(1,1,1):

  @code{.cpp}
  BoxLoops::loop<D_DECL(1, 1, 1)>(box, [&](const IntVect& iv) {
      phi(iv, 0) = 0.0;
  });

  BoxLoops::loop<D_DECL(2, 2, 2)>(colorBox, [&](const IntVect& iv) {
      phi(iv, 0) = 0.0;
  });
  @endcode

  <b>VoF, face, and set iterators:</b>

  @code{.cpp}
  BoxLoops::loop(vofIt,  [&](const VolIndex&  vof)  { ... });
  BoxLoops::loop(faceIt, [&](const FaceIndex& face) { ... });
  BoxLoops::loop(ivs,    [&](const IntVect&   iv)   { ... });
  @endcode
*/
namespace BoxLoops {

  /*!
    @brief Launch a C++ kernel over a regular grid with compile-time per-dimension strides.
    @details Strides are encoded as non-type template parameters, making them visible to the
    compiler as constants and enabling auto-vectorization of the innermost loop. All strides must
    be >= 1 (enforced by @c static_assert).

    Note: @c IntVect cannot be a non-type template parameter in C++14, so per-dimension @c int
    parameters are used instead. Call sites can use @c D_DECL for dimension independence, e.g.:

    @code{.cpp}
    BoxLoops::loop<D_DECL(2, 2, 2)>(box, kernel);
    @endcode

    @param[in]    a_computeBox Computation box
    @param[inout] a_kernel     Kernel to launch.
  */
#if CH_SPACEDIM == 2
  template <int Si, int Sj, typename Functor>
  ALWAYS_INLINE void
  loop(const Box& a_computeBox, Functor&& kernel);
#else
  template <int Si, int Sj, int Sk, typename Functor>
  ALWAYS_INLINE void
  loop(const Box& a_computeBox, Functor&& kernel);
#endif

  /*!
    @brief Launch a C++ kernel over a subset of cells.
    @param[inout] a_cells   Grid cells where we launch the kernel.
    @param[inout] a_kernel  Kernel to launch
  */
  template <typename Functor>
  ALWAYS_INLINE void
  loop(const IntVectSet& a_ivs, Functor&& a_kernel);

  /*!
    @brief Launch a C++ kernel over a subset of cells.
    @param[inout] a_cells   Grid cells where we launch the kernel.
    @param[inout] a_kernel  Kernel to launch
  */
  template <typename Functor>
  ALWAYS_INLINE void
  loop(const DenseIntVectSet& a_ivs, Functor&& a_kernel);

  /*!
    @brief Launch a C++ kernel over cut-cells.
    @param[inout] a_vofit   Grid cells where we launch the kernel.
    @param[inout] a_kernel  Kernel to launch
  */
  template <typename Functor>
  ALWAYS_INLINE void
  loop(VoFIterator& a_iter, Functor&& a_kernel);

  /*!
    @brief Launch a C++ kernel over faces
    @param[inout] a_faceit  Grid faces where we launch the kernel.
    @param[inout] a_kernel  Kernel to launch
  */
  template <typename Functor>
  ALWAYS_INLINE void
  loop(FaceIterator& a_iter, Functor&& a_kernel);

  /*!
    @brief Launch a C++ kernel over a specified subset of something.
    @param[inout] a_subset Subset. E.g. a std::vector<VolIndex>
    @param[inout] a_kernel Kernel to launch
  */
  template <typename T, typename Functor>
  ALWAYS_INLINE void
  loop(const Vector<T>& a_subset, Functor&& a_kernel);
} // namespace BoxLoops

#include <CD_NamespaceFooter.H>

#include <CD_BoxLoopsImplem.H>

#endif
