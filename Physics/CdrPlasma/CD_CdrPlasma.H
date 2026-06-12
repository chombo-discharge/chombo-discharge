/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_CdrPlasma.H
   @brief  Declaration of the CdrPlasma namespace
   @author Robert Marskar
*/

#ifndef CD_CDRPLASMA_H
#define CD_CDRPLASMA_H

namespace Physics {

  /**
     @brief Namespace for encapsulating CDR plasma physics.
     @details This namespace contains a physics module for simulating weakly ionized plasmas
     in the drift-diffusion (CDR) approximation, self-consistently coupled to Poisson's equation
     for the electric field and to radiative transfer equations for photoionization. It includes
     abstract and concrete time steppers, cell taggers, plasma kinetics interfaces, and JSON-driven
     plasma models.
  */
  namespace CdrPlasma {};

} // namespace Physics

#endif
