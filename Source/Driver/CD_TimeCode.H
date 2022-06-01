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

/*!
  @brief Enum class for identifying how a time step was restricted. 
  @details If the time code you are looking for is not in here, put it in and submit a pull request. 
*/
enum class TimeCode
{
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
