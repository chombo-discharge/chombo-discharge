/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzDirichletEBBC.cpp
  @brief  Implementation of CD_MFHelmholtzDirichletEBBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzDirichletEBBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzDirichletEBBC::MFHelmholtzDirichletEBBC(const int a_phase, const RefCountedPtr<JumpBC>& a_jumpBC){
  m_phase  = a_phase;
  m_jumpBC = a_jumpBC;

  m_order       = -1;
  m_weight      = -1;
  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzDirichletEBBC::~MFHelmholtzDirichletEBBC(){

}

#include <CD_NamespaceFooter.H>
