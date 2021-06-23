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

MFHelmholtzDirichletEBBC::MFHelmholtzDirichletEBBC() : EBHelmholtzDirichletEBBC() {

}

MFHelmholtzDirichletEBBC::~MFHelmholtzDirichletEBBC(){

}

void MFHelmholtzDirichletEBBC::define() {
  EBHelmholtzDirichletEBBC::define();
}

#include <CD_NamespaceFooter.H>
