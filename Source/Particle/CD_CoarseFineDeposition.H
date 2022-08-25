/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CoarseFineDeposition.H
  @brief  Declaration of deposition algorithms in an AMR context. 
  @author Robert Marskar
*/

#ifndef CD_CoarseFineDeposition_H
#define CD_CoarseFineDeposition_H

// Our includes
#include <CD_NamespaceHeader.H>

/*!
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
  HaloNGP
};

#include <CD_NamespaceFooter.H>

#endif
