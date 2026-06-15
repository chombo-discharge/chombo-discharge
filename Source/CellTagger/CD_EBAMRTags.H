/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_EBAMRTags.H
  @brief  Contains an alias declaration of EBAMRTags.
  @author Robert Marskar
*/

#ifndef CD_EBAMRTAGS_H
#define CD_EBAMRTAGS_H

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Declaration of cell tags
  @details EBAMRTags is an AMR alias for storing a DenseIntVectSet in each grid patch. It is used to flag cells for refinement, but only within grid patches.
  @note DenseIntVectSet is used for performance reasons.
*/
using EBAMRTags = Vector<RefCountedPtr<LayoutData<DenseIntVectSet>>>;

#include <CD_NamespaceFooter.H>

#endif
