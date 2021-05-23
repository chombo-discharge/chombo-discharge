/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PerlinSlabSdf.H
  @brief  Declaration of a noisy slab, i.e. a box with surface roughness
  @author Robert Marskar
*/

#ifndef CD_PerlinSlabSdf_H
#define CD_PerlinSlabSdf_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Implicit function class describing a hexagonal rod
*/
class PerlinSlabSdf: public BaseIF {
public:

  /*!
    @brief Constructor
    @param[in] a_ccPoint     Center point on the slab side with surface noise
    @param[in] a_normal      Normal vector of slab side with surface noise
    @param[in] a_xyz         Width, depth, and height of slab
    @param[in] a_curv        Corner curvatures
    @param[in] a_fluidInside Fluid on inside or not
  */
  PerlinSlabSdf(const RealVect a_ccPoint,
		 const RealVect a_normal,
		 const RealVect a_xyz,
		 const RealVect a_noiseFreq,
		 const int      a_octaves,
		 const Real     a_noiseAmp,
		 const Real     a_persistence,
		 const Real     a_cornerCurv,
		 const bool     a_reseed,
		 const bool     a_fluidInside);

  PerlinSlabSdf(const PerlinSlabSdf& a_inputIF);
  PerlinSlabSdf() = delete;

  virtual ~PerlinSlabSdf();
  virtual Real value(const RealVect& a_pos) const override;
  virtual BaseIF* newImplicitFunction() const override;

protected:

  RefCountedPtr<BaseIF> m_baseif;
  bool m_fluidInside;

};

#include <CD_NamespaceFooter.H>

#endif
