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
*/
namespace BoxLoops {

  /*!
    @brief Launch a C++ kernel over a regular grid. 
    @details The stride argument lets the user use different strides along the different coordinate axes. The default 
    is to use a stride of one along each coordinate. Note that the stride must be positive for each coordinate axis. 
    @param[in]    a_computeBox Computation box
    @param[inout] a_kernel     Kernel to launch. 
    @param[in]    a_stride     Stride for each loop. 
  */
  template <typename Functor>
  ALWAYS_INLINE void
  loop(const Box& a_computeBox, Functor&& kernel, const IntVect& a_stride = IntVect::Unit);

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
