/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BoxLoops.H
  @brief  Agglomeration of useful data operations
  @author Robert Marskar
*/

#ifndef CD_BoxLoops_H
#define CD_BoxLoops_H

// Chombo includes
#include <Box.H>
#include <VoFIterator.H>

// Our includes
#include <CD_Decorations.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Namespace for encapsulating loop types
  @note All methods in this class assumes that data are allocated over the same realm. If you are trying to do something across realms, you need to copy
  the data onto realms that lie on the same data holders first!.
*/
namespace BoxLoops {
  
  /*!
    @brief Launch a C++ kernel over a regular grid. This version uses simd decorations to ensure that the innermost loop can be vectorized. 
    @param[in] a_computeBox Computation box
    @param[in] a_kernel     Kernel to launch. 
  */
  template <typename Functor>
  ALWAYS_INLINE
  static void loop(const Box& a_computeBox, Functor&& kernel);

  /*!
    @brief Launch a C++ kernel over a cut-cells. 
    @param[in] a_vofit   Grid cells where we launch the kernel. 
    @param[in] a_kernel  Kernel to launch
  */  
  template <typename Functor>
  ALWAYS_INLINE
  static void loop(VoFIterator& a_vofit, Functor&& a_kernel);

}

#include <CD_NamespaceFooter.H>

#include <CD_BoxLoopsImplem.H>

#endif
