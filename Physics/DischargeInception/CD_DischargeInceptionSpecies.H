/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_DischargeInceptionSpecies.H
   @brief  Declaration of the Physics::DischargeInception::DischargeInceptionSpecies CDR species
   @author Robert Marskar
*/

#ifndef CD_DISCHARGEINCEPTIONSPECIES_H
#define CD_DISCHARGEINCEPTIONSPECIES_H

// Our includes
#include <CD_CdrSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace DischargeInception {

    /**
       @brief CdrSpecies subclass for use with DischargeInceptionStepper.
       @details Represents a charged or neutral species whose initial density, mobility, and
       diffusivity are all set at construction time via function and flag arguments. The initial
       density is provided as a callable that maps a position to a density value.
    */
    class DischargeInceptionSpecies : public CdrSpecies
    {
    public:
      /**
         @brief Full constructor.
         @param[in] a_initialData Callable returning the initial density at a given position.
         @param[in] a_mobile      True if the species has a non-zero drift velocity.
         @param[in] a_diffusive   True if the species has a non-zero diffusion coefficient.
      */
      DischargeInceptionSpecies(const std::function<Real(const RealVect& a_pos)>& a_initialData,
                                const bool                                        a_mobile,
                                const bool                                        a_diffusive)
      {
        m_initialData  = a_initialData;
        m_isMobile     = a_mobile;
        m_isDiffusive  = a_diffusive;
        m_chargeNumber = 0;
        m_name         = "Ion distribution";
      }

      /**
         @brief Destructor.
      */
      virtual ~DischargeInceptionSpecies()
      {}

      /**
         @brief Return the initial species density at the given position.
         @param[in] a_pos  Physical position.
         @param[in] a_time Current time (not used; initial data is time-independent).
         @return Initial density value from the stored function.
      */
      Real
      initialData(const RealVect a_pos, const Real a_time) const override
      {
        return m_initialData(a_pos);
      }

    protected:
      /**
         @brief Callable that maps a position to the initial species density.
      */
      std::function<Real(const RealVect& a_pos)> m_initialData;
    };
  } // namespace DischargeInception
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
