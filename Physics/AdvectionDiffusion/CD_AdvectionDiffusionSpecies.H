/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_AdvectionDiffusionSpecies.H
  @brief  Declaration of the Physics::AdvectionDiffusion::AdvectionDiffusionSpecies class
  @author Robert Marskar
*/

#ifndef CD_ADVECTIONDIFFUSIONSPECIES_H
#define CD_ADVECTIONDIFFUSIONSPECIES_H

// Our includes
#include <CD_CdrSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace AdvectionDiffusion {

    /**
      @brief Concrete CdrSpecies implementation for the advection-diffusion test problem.
      @details This class defines a single passive scalar species for use with
      AdvectionDiffusionStepper. The initial data, mobility, and diffusivity are set at
      construction time. The species carries no net charge.
    */
    class AdvectionDiffusionSpecies : public CdrSpecies
    {
    public:
      /**
        @brief Full constructor.
        @param[in] a_initialData Functor returning the initial scalar concentration at a given position.
        @param[in] a_mobile      If true, the species is advected by the velocity field.
        @param[in] a_diffusive   If true, the species undergoes diffusion.
      */
      AdvectionDiffusionSpecies(const std::function<Real(const RealVect& a_pos)>& a_initialData,
                                const bool                                        a_mobile,
                                const bool                                        a_diffusive)
      {
        m_initialData  = a_initialData;
        m_isMobile     = a_mobile;
        m_isDiffusive  = a_diffusive;
        m_chargeNumber = 0;
        m_name         = "AdvectionDiffusionSpecies";
      }

      /**
        @brief Destructor.
      */
      virtual ~AdvectionDiffusionSpecies()
      {}

      /**
        @brief Return the initial scalar concentration at a given position and time.
        @param[in] a_pos  Physical position.
        @param[in] a_time Current time (unused; initial data is time-independent).
        @return Initial concentration evaluated at a_pos.
      */
      Real
      initialData(const RealVect a_pos, const Real a_time) const override
      {
        return m_initialData(a_pos);
      }

    protected:
      /**
        @brief Functor encoding the initial concentration profile.
      */
      std::function<Real(const RealVect& a_pos)> m_initialData;
    };
  } // namespace AdvectionDiffusion
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
