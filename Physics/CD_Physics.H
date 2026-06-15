/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_Physics.H
   @brief  Declaration of the Physics namespace
   @author Robert Marskar
*/

#ifndef CD_PHYSICS_H
#define CD_PHYSICS_H

/**
   @brief Namespace containing physics models for use with chombo-discharge.

   @details This namespace collects self-contained physics modules built on top of
   the core chombo-discharge solvers and infrastructure. Each sub-module defines a complete
   problem specification including species, time steppers, and cell taggers.
*/
namespace Physics {};

#endif
