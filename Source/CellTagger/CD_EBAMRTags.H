/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */
/*!
  @file   CD_EBAMRTags.H
  @brief  Simple declaration of tags
  @author Robert Marskar
*/

#ifndef CD_EBAMRTags_H
#define CD_EBAMRTags_H

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of cell tags
*/
using EBAMRTags = Vector<RefCountedPtr<LayoutData<DenseIntVectSet> > >;

#include <CD_NamespaceFooter.H>

#endif
