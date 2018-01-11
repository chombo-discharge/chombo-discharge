/*!
  @file perlin_sphere_if.H
  @brief Declaration of a perlin_sphere_if function
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _PERLIN_SPHERE_IF_
#define _PERLIN_SPHERE_IF_

#include <BaseIF.H>
#include <SphereIF.H>

#include "perlin_if.H"

/*!
  @brief perlin_sphere_if function
*/
class perlin_sphere_if : public BaseIF {
public:

  /*!
   */
  perlin_sphere_if(){}
  /*!
    @brief Constructor
  */
  perlin_sphere_if(const Real&     a_rad,
		   const RealVect& a_center,
		   const bool&     a_inside,
		   const Real&     a_noiseAmp,
		   const RealVect& a_noiseFreq,
		   const Real&     a_persistence,
		   const int&      a_octaves,
		   const bool&     a_reseed = false);

  /*!
    @brief Copy constructor
  */
  perlin_sphere_if(const perlin_sphere_if& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~perlin_sphere_if();

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
    @brief Radius
  */
  Real m_rad;

  /*!
    @brief Center
  */
  RealVect m_center;

  /*!
    @brief Inside/outside
  */
  bool m_inside;
  
  /*!
    @brief Noise function
  */
  RefCountedPtr<BaseIF> m_perlinIF;
};
#endif
