/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */
/*!
  @file   CD_EBAMRTags.H
  @brief  Contains an alias declaration of EBAMRTags.
  @author Robert Marskar
*/

#ifndef CD_EBAMRTags_H
#define CD_EBAMRTags_H

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of cell tags
  @details EBAMRTags is an AMR alias for storing a DenseIntVectSet in each grid patch. It is used to flag cells for refinement, but only within grid patches. 
  @note DenseIntVectSet is sued for performance reasons. 
*/
using EBAMRTags = Vector<RefCountedPtr<LayoutData<DenseIntVectSet> > >;

#include <CD_NamespaceFooter.H>

#endif
