/*!
  @file   ito_plasma_air2.H
  @brief  Simple class for 
  @author Robert Marskar
  @date   June 2020
*/

#ifndef _ITO_PLASMA_AIR2_
#define _ITO_PLASMA_AIR2_

#include "ito_plasma_physics.H"

namespace physics {
  namespace ito_plasma{

    /*!
      @brief Abstract clas for specifying plasma kinetics
    */
    class ito_plasma_air2 : public ito_plasma_physics {
    public:
      ito_plasma_air2();
      ~ito_plasma_air2();

      /*!
	@brief Compute Townsend ionization length
      */
      Real compute_alpha(const RealVect a_E) const override;

      Vector<RealVect> compute_ito_velocities(const Real         a_time,
					      const RealVect     a_pos,
					      const RealVect     a_E,
					      const Vector<Real> a_cdr_densities) const override;
    };
  }
}
#endif
