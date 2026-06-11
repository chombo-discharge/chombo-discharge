/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_PointParticle.H
  @brief  Declaration of a computational point particle 
  @author Robert Marskar
*/

#ifndef CD_POINTPARTICLE_H
#define CD_POINTPARTICLE_H

// Std includes
#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

// Our includes
#include <CD_GenericParticle.H>
#include <CD_NamespaceHeader.H>

/**
  @brief A particle class that only has a position and a weight.
*/
class PointParticle : public GenericParticle<1, 0>
{
public:
  /**
    @brief Default constructor. Sets weight and position to zero.
  */
  inline PointParticle();

  /**
    @brief Constructor. Subsequently calls define
    @param[in] a_position Position (in physical space)
    @param[in] a_weight     Particle weight
  */
  inline PointParticle(const RealVect& a_position, Real a_weight);

  /**
    @brief Dtor - does nothing.
  */
  inline ~PointParticle() override;

  /**
    @brief Copy constructor.
    @param[in] a_other Other particle.
  */
  PointParticle(const PointParticle& a_other) = default;

  /**
    @brief Copy assignment operator.
    @param[in] a_other Other particle.
    @return Reference to this.
  */
  PointParticle&
  operator=(const PointParticle& a_other) = default;

  /**
    @brief Move constructor.
    @param[in] a_other Other particle.
  */
  PointParticle(PointParticle&& a_other) = default;

  /**
    @brief Move assignment operator.
    @param[in] a_other Other particle.
    @return Reference to this.
  */
  PointParticle&
  operator=(PointParticle&& a_other) = default;

  /**
    @brief Full define function
    @param[in] a_position Position (in physical space)
    @param[in] a_weight     Particle weight
  */
  inline void
  define(const RealVect& a_position, Real a_weight);

  /**
    @brief Get weight
    @return Return value
  */
  inline Real&
  weight();

  /** 
    @brief Const retrieve weight
    @return Return value
  */
  inline const Real&
  weight() const;
};

#include <CD_NamespaceFooter.H>

#include <CD_PointParticleImplem.H>

#endif
