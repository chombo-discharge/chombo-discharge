/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_TimeCode.H
  @brief  Declaration of various identifiers for "time codes"
  @author Robert Marskar
*/

#ifndef CD_TIMECODE_H
#define CD_TIMECODE_H

// Our includes
#include <CD_NamespaceHeader.H>

/**
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
