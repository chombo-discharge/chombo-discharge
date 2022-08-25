/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PointParticle.H
  @brief  Declaration of a computational point particle 
  @author Robert Marskar
*/

#ifndef CD_PointParticle_H
#define CD_PointParticle_H

// Std includes
#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

// Our includes
#include <CD_GenericParticle.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief A particle class that only has a position and a weight.
*/
class PointParticle : public GenericParticle<1, 0>
{
public:
  /*!
    @brief Default constructor. Sets weight and position to zero.
  */
  inline PointParticle();

  /*!
    @brief Constructor. Subsequently calls define
    @param[in] a_position Position (in physical space)
    @param[in] a_weight     Particle weight
  */
  inline PointParticle(const RealVect a_position, const Real a_weight);

  /*!
    @brief Dtor - does nothing.
  */
  inline virtual ~PointParticle();

  /*!
    @brief Full define function
    @param[in] a_position Position (in physical space)
    @param[in] a_weight     Particle weight
  */
  inline void
  define(const RealVect a_position, const Real a_weight);

  /*!
    @brief Get weight
  */
  inline Real&
  weight();

  /*! 
    @brief Const retrieve weight
  */
  inline const Real&
  weight() const;
};

#include <CD_NamespaceFooter.H>

#include <CD_PointParticleImplem.H>

#endif
