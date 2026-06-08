/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_BrownianWalker.H
   @brief  Declaration of the BrownianWalker namespace
   @author Robert Marskar
*/

#ifndef CD_BROWNIANWALKER_H
#define CD_BROWNIANWALKER_H

namespace Physics {

  /**
     @brief Namespace for encapsulating Brownian walker physics.
     @details This namespace contains a particle-based physics module in which particles undergo
     stochastic Brownian motion and deterministic advection on an adaptive mesh. It includes a
     species definition, a time stepper, and a cell tagger.
  */
  namespace BrownianWalker {}

} // namespace Physics

#endif
