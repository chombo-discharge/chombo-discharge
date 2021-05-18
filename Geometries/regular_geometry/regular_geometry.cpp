/*!
  @file   regular_geometry.H
  @brief  Regular geometry (no EBs)
  @author Robert Marskar
  @date   July 2018
*/

#include "regular_geometry.H"

#include "CD_NamespaceHeader.H"

regular_geometry::regular_geometry(){
  m_dielectrics.resize(0);
  m_electrodes.resize(0);
}

regular_geometry::~regular_geometry(){
  
}
#include "CD_NamespaceFooter.H"
