/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file    CD_MushroomIF.H
  @brief   Declaration of a "mushroom" implicit function, i.e. a cylinder with a circular flat plate on top. 
  @author  Robert Marskar
*/

#ifndef CD_MUSHROOMIF_H
#define CD_MUSHROOMIF_H

// Chombo includes
#include <BaseIF.H>
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief "Mushroom" implicit function. Consists of a cylinder with a circular plate on top. 
  @details Oriented along +z in 3D and +y in 2D. 
*/
class MushroomIF : public BaseIF
{
public:
  /**
    @brief Full constructor
    @param[in] a_center      Center point of the circular platew
    @param[in] a_R           Radius of the flat plate. 
    @param[in] a_r           Radius of the cylinder foot
    @param[in] a_L           Length of the cylinder foot
    @param[in] a_d           Plate thickness
    @param[in] a_curv        Corner curvature
    @param[in] a_fluidInside Fluid on inside or not
  */
  MushroomIF(const RealVect& a_center,
             const Real      a_R,
             const Real      a_r,
             const Real      a_L,
             const Real      a_d,
             const Real      a_curv,
             const bool      a_fluidInside);

  /**
    @brief Copy constructor
    @param[in] a_inputIF Other input
  */
  MushroomIF(const MushroomIF& a_inputIF);

  /**
    @brief Destructor (does nothing)
  */
  virtual ~MushroomIF();

  /**
    @brief Value function
    @param[in] a_pos Position
    @return Function value at the given point
  */
  virtual Real
  value(const RealVect& a_pos) const override;

  /**
    @brief Factory function.
    @return Returns new MushroomIF object, casted to BaseIF*
  */
  virtual BaseIF*
  newImplicitFunction() const override;

protected:
  /**
    @brief Implicit function -- this is a smooth union between the foot and the circular plate
  */
  RefCountedPtr<BaseIF> m_baseIF;

  /**
    @brief Fluid on the inside or not
  */
  bool m_fluidInside;
};

#include <CD_NamespaceFooter.H>

#endif
