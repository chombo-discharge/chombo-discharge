/* chombo-discharge
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NeedleIF.cpp
  @brief  Implementation of CD_NeedleIF.H
  @author Fanny Skirbekk
*/

// Our includes
#include <CD_NeedleIF.H>
#include <CD_CylinderSdf.H>
#include <CD_NamespaceHeader.H>

NeedleIF::NeedleIF(const RealVect& a_centerTip, const RealVect& a_centerBack, const Real& a_radius, const bool& a_fluidInside, const Real& a_tipLength, const Real& a_tipRadius){

  // in rodif they find new centers --why?

  // Build the needle-parts
  Vector<BaseIF*> isects;
  isects.push_back(static_cast<BaseIF*> (new CylinderSdf(a_centerTip, a_centerBack, a_radius, a_fluidInside)));

  // Build the needle
  m_baseif = RefCountedPtr<BaseIF>(new IntersectionIF(isects));

  // Delete everything we have allocated so far
  for(int i = 0; i < isects.size(); ++i){
    delete isects[i];
  }
}

NeedleIF::NeedleIF(const NeedleIF& a_inputIF){
  this->m_baseif = a_inputIF.m_baseif;
}

#include <CD_NamespaceFooter.H>
