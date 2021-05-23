/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TimeCode.H
  @brief  Declaration of various identifiers for "time codes"
  @author Robert Marskar
*/

#ifndef CD_TimeCode_H
#define CD_TimeCode_H

// Our includes
#include <CD_NamespaceHeader.H>
  
enum class TimeCode {
  Advection,
  Diffusion,
  AdvectionDiffusion,  
  Source,
  RelaxationTime,
  Restricted,
  Hardcap,
  Error,
  Physics
};

#include <CD_NamespaceFooter.H>

#endif
