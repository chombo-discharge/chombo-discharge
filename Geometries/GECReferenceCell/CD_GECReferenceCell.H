/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_GECReferenceCell.H
  @brief  Simplified GEC reference cell geometry.
  @author Robert Marskar
*/

#ifndef CD_GECREFERENCECELL_H
#define CD_GECREFERENCECELL_H

// our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Simplified GEC reference cell geometry.
*/
class GECReferenceCell : public ComputationalGeometry
{
public:
  /**
    @brief Constructor
  */
  GECReferenceCell();

  /**
    @brief Destructor
  */
  ~GECReferenceCell();
};

#include <CD_NamespaceFooter.H>

#endif
