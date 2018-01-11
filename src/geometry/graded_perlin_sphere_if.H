/*!
  @file graded_perlin_sphere_if.H
  @brief Declaration of a graded_perlin_sphere_if function
  @author Robert Marskar
  @date Oct. 2017
*/

#ifndef _GRADED_PERLIN_SPHERE_IF_
#define _GRADED_PERLIN_SPHERE_IF_

#include "perlin_sphere_if.H"

/*!
  @brief graded_perlin_sphere_if function
*/
class graded_perlin_sphere_if : public perlin_sphere_if {
public:  
  /*!
    @brief Constructor
  */
  graded_perlin_sphere_if(const Real&     a_rad,
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

  graded_perlin_sphere_if(const graded_perlin_sphere_if& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~graded_perlin_sphere_if();

  /*!
    @brief Value function
  */
  virtual Real value(const RealVect& a_pos) const;

  /*!
    @brief Factory function
  */
  virtual BaseIF* newImplicitFunction() const;
  
protected:
};
#endif
