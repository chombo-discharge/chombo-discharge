/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Tessellation.H
  @brief  Declaration of a geometry class for surface tessellation.
  @author Robert Marskar
*/

#ifndef CD_TESSELLATION_H
#define CD_TESSELLATION_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief A geometry consisting of a single polygon mesh read from an STL file.
*/
class Tessellation : public ComputationalGeometry
{
public:
  /**
    @brief Constructor. Reads mesh file and input options via ParmParse.
  */
  Tessellation();

  /**
    @brief Destructor.
  */
  ~Tessellation();
};

#include <CD_NamespaceFooter.H>

#endif
