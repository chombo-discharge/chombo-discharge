/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_AdvectionDiffusion.H
   @brief  Declaration of the AdvectionDiffusion namespace
   @author Robert Marskar
*/

#ifndef CD_ADVECTIONDIFFUSION_H
#define CD_ADVECTIONDIFFUSION_H

namespace Physics {

  /**
     @brief Namespace for encapsulating advection-diffusion physics.

     @details This namespace contains a simple model for advection-diffusion of a passive scalar on
     an adaptive mesh. It includes a time stepper, a species definition, and a cell tagger.
  */
  namespace AdvectionDiffusion {};

} // namespace Physics

#endif
