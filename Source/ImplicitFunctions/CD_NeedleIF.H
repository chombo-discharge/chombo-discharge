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

#ifndef CD_NeedleIF_H
#define CD_NeedleIF_H

// Chombo includes
#include <BaseIF.H>
#include <IntersectionIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

class NeedleIF : public BaseIF{
public:
  NeedleIF(const RealVect& a_centerTip, const RealVect& a_centerBack, const Real& a_radius, const bool& a_fluidInside, const Real& a_tipLength, const Real& a_tipRadius);

  NeedleIF(const NeedleIF& a_inputIF);

protected:
  RefCountedPtr<BaseIF> m_baseif;
};

#include <CD_NamespaceFooter.H>

#endif
