/*!
  @file   perlin_slab.H
  @brief  Declaration of a noisy slab, i.e. a box with surface roughness
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _PERLIN_SLAB_IF_
#define _PERLIN_SLAB_IF_

#include <BaseIF.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Implicit function class describing a hexagonal rod
*/
class perlin_slab_if: public BaseIF {
public:

  /*!
    @brief Constructor
    @param[in] a_ccPoint     Center point on the slab side with surface noise
    @param[in] a_normal      Normal vector of slab side with surface noise
    @param[in] a_xyz         Width, depth, and height of slab
    @param[in] a_curv        Corner curvatures
    @param[in] a_fluidInside Fluid on inside or not
  */
  perlin_slab_if(const RealVect a_ccPoint,
		 const RealVect a_normal,
		 const RealVect a_xyz,
		 const RealVect a_noiseFreq,
		 const int      a_octaves,
		 const Real     a_noiseAmp,
		 const Real     a_persistence,
		 const Real     a_cornerCurv,
		 const bool     a_reseed,
		 const bool     a_fluidInside);

  perlin_slab_if(const perlin_slab_if& a_inputIF);
  perlin_slab_if() = delete;

  virtual ~perlin_slab_if();
  virtual Real value(const RealVect& a_pos) const override;
  virtual BaseIF* newImplicitFunction() const override;

protected:

  RefCountedPtr<BaseIF> m_baseif;
  bool m_fluidInside;

};
#include "CD_NamespaceFooter.H"
#endif
