/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_BrownianWalkerSpecies.H
   @brief  Declaration of the Physics::BrownianWalker::BrownianWalkerSpecies class
   @author Robert Marskar
*/

#ifndef CD_BROWNIANWALKERSPECIES_H
#define CD_BROWNIANWALKERSPECIES_H

// Std includes
#include <random>

// Our includes
#include <CD_ItoSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace BrownianWalker {

    /**
       @brief Concrete ItoSpecies for the Brownian walker physics module.
       @details Initial particles are drawn from a Gaussian distribution centered at a
       user-specified position with a user-specified radius. The number of particles, blob
       center, and blob radius are read from ParmParse at construction time. Each particle
       is assigned unit weight.
    */
    class BrownianWalkerSpecies : public ItoSpecies
    {
    public:
      /**
         @brief Constructor. Reads input parameters and draws the initial particle distribution.
      */
      BrownianWalkerSpecies();

      /**
         @brief Destructor.
      */
      ~BrownianWalkerSpecies();

    protected:
      /**
         @brief Number of initial computational particles to draw.
      */
      int m_numParticles;

      /**
         @brief Standard deviation of the initial Gaussian particle distribution.
      */
      Real m_blobRadius;

      /**
         @brief Center of the initial Gaussian particle distribution.
      */
      RealVect m_blobCenter;

      /**
         @brief Draw the initial particle positions from a Gaussian distribution and assign unit weights.
      */
      void
      drawInitParticles();
    };
  } // namespace BrownianWalker
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
