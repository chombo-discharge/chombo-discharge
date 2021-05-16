/*!
  @file   advection_diffusion_species.H
  @brief  Simple species for advection_diffusion test problem
  @author Robert Marskar
  @date   March 2020
*/

#ifndef _ADVECTION_DIFFUSION_SPECIES_
#define _ADVECTION_DIFFUSION_SPECIES_

#include "cdr_species.H"

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace advection_diffusion {

    class advection_diffusion_species : public cdr_species {
    public:
      advection_diffusion_species();
      ~advection_diffusion_species();

      Real initialData(const RealVect a_pos, const Real a_time) const;
    protected:
      Real     m_blob_radius;
      Real     m_blob_amplitude;
      RealVect m_blob_center;
    
    };
  }
}
#include "CD_NamespaceFooter.H"

#endif
