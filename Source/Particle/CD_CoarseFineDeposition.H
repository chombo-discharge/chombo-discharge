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
*/
enum class CoarseFineDeposition {
  PVR,
  Halo
};

#include <CD_NamespaceFooter.H>

#endif
