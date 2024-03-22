/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   RegularGeometry.cpp
  @brief  Implementation of CD_RegularGeometry.H
  @author Robert Marskar
*/

// Our includes
#include <CD_RegularGeometry.H>
#include <CD_NamespaceHeader.H>

RegularGeometry::RegularGeometry()
{
  m_dielectrics.resize(0);
  m_electrodes.resize(0);
}

RegularGeometry::~RegularGeometry()
{}

#include <CD_NamespaceFooter.H>
