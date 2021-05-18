/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_AdvectionDiffusionSpecies.H
  @brief  Simple species for AdvectionDiffusion test problem
  @author Robert Marskar
*/

#ifndef CD_AdvectionDiffusionSpecies_H
#define CD_AdvectionDiffusionSpecies_H

// Our includes
#include <CD_CdrSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace AdvectionDiffusion {

    class AdvectionDiffusionSpecies : public CdrSpecies {
    public:
      AdvectionDiffusionSpecies();
      ~AdvectionDiffusionSpecies();

      Real initialData(const RealVect a_pos, const Real a_time) const override;
    protected:
      Real     m_blob_radius;
      Real     m_blob_amplitude;
      RealVect m_blob_center;
    
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
