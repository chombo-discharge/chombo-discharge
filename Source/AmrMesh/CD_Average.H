/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Average.H
  @brief  Declaration of averaging methods.
  @author Robert Marskar
*/

#ifndef CD_AVERAGE_H
#define CD_AVERAGE_H

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Various averaging methods. 
  @details 'Arithmetic'/'Harmonic' mean what you think, while 'Conservative' also incorporate volume, face, 
  and EB face fractions such that the total quantity is conserved. 
*/
enum class Average
{
  Arithmetic,
  Harmonic,
  Conservative,
  Geometric
};

#include <CD_NamespaceFooter.H>

#endif
