/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_CoarseFineDeposition.H
  @brief  Declaration of deposition algorithms in an AMR context. 
  @author Robert Marskar
*/

#ifndef CD_COARSEFINEDEPOSITION_H
#define CD_COARSEFINEDEPOSITION_H

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Coarse-fine deposition types (see CD_EBAMRParticleMesh for how these are handled).
  @details Interp means that the mass on the fine level is put there via interpolation of the coarse level.
  @details Halo means that we deposit the particles on the coarse side of the refinement boundary directly on the fine level, using
           special deposition routines where the particles have the same physical width as on the coarse level. 
  @details HaloNGP means that particles on the coarse side of the refinement boundary are deposited using an NGP scheme. 
*/
enum class CoarseFineDeposition
{
  Interp,
  Halo,
  HaloNGP,
  Transition
};

#include <CD_NamespaceFooter.H>

#endif
