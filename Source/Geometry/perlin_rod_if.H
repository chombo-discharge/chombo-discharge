/*!
  @file perlin_rod_if.H
  @brief Declaration of a rod with a rounded, noisy tip
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _PERLIN_ROD_IF_
#define _PERLIN_ROD_IF_

#include <BaseIF.H>

#include <CD_CylinderSdf.H>
#include <CD_GradedPerlinSphereSdf.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief perlin_rod_if function
*/
class perlin_rod_if : public PerlinSphereSdf {
public:  
  /*!
    @brief Constructor
  */
  perlin_rod_if(const Real&     a_rad,
		const RealVect& a_center1,
		const RealVect& a_center2,
		const bool&     a_inside,
		const Real&     a_noiseAmp,
		const RealVect& a_noiseFreq,
		const Real&     a_persistence,
		const int&      a_octaves,
		const bool&     a_reseed = false);

  /*!
    @brief Copy constructor
  */

  perlin_rod_if(const perlin_rod_if& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~perlin_rod_if();

  /*!
    @brief Value function
  */
  virtual Real value(const RealVect& a_pos) const;

  /*!
    @brief Factory function
  */
  virtual BaseIF* newImplicitFunction() const;
  
protected:

  /*!
    @brief Basis rod, used for copying
  */
  RefCountedPtr<BaseIF> m_baseif;
};
#include "CD_NamespaceFooter.H"
#endif
