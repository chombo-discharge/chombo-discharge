/* chombo-discharge
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NeedleIF.H
  @brief  Declaration of needle electrode (cylinder with one pointy end and one flat end)
  @author Fanny Skirbekk
*/

#ifndef CD_NeedleIF.H
#define CD_NeedleIF.H

// Chombo includes
#include <BaseIF.H>
#include <IntersectionIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Cylinder with one pointed end and one flat end
*/
class NeedleIF : public BaseIF{
public:
  /*!
    @brief Constructor
    @param[in] a_centerTipSide
  */
  NeedleIFconst RealVect& a_centerTipSide, const Realvect& a_centerBack, const Real& a_radius, const bool& a_fluidInside, const Real& a_tipRadius, RealVect& a_tip, Real& a_angle);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input Needle
  */
  NeedleIF(const NeedleIF& a_inputIF);

  /*!
    @brief Value function. Returns distance to object.
    @param[in] a_point Physical position.
  */
  virtual Real value(const RealVect& a_point) const override;

  /*!
    @brief IF factory method
  */
  virtual BaseIF* newImpicitFunction() const override;

protected:
  RefCountedPtr<BaseIF> m_baseif;
};

#include <CD_NamespaceFooter.H>

#endif
