/*!
  @file perlin_plane_if.H
  @brief Declaration of a noise plane level-set function
  @details This class is essentially a PlaneIF with PerlinSdf noise
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _PERLIN_PLANE_IF_
#define _PERLIN_PLANE_IF_

#include <BaseIF.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief perlin_plane_if function
*/
class perlin_plane_if : public BaseIF {
public:  
  /*!
    @brief Constructor
  */
  perlin_plane_if(const RealVect a_normal,
		  const RealVect a_point,
		  const bool     a_inside,
		  const Real     a_noiseAmp,
		  const RealVect a_noiseFreq,
		  const Real     a_persistence,
		  const int      a_octaves,
		  const bool     a_reseed = false);
  /*!
    @brief Copy constructor
  */
  perlin_plane_if(const perlin_plane_if& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~perlin_plane_if();

  /*!
    @brief Level-set function
  */
  virtual Real value(const RealVect& a_pos) const;

  /*!
    @brief Factory method
  */
  virtual BaseIF* newImplicitFunction() const;
  
protected:

  /*!
    @brief Plane normal
  */
  RealVect m_normal;

  /*!
    @brief Plane point
  */
  RealVect m_point;

  /*!
    @brief Rod function
  */
  RefCountedPtr<BaseIF> m_plane;

  /*!
    @brief Noise function
  */
  RefCountedPtr<BaseIF> m_perlin;
};
#include "CD_NamespaceFooter.H"
#endif
