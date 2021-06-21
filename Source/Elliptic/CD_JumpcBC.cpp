/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_JumpBC.cpp
  @brief  Implementatino of CD_JumpBC.H
  @author Robert Marskar
*/


// Our includes
#include <CD_JumpBC.H>
#include <CD_NamespaceHeader.H>

JumpBC::JumpBC(const MFLevelGrid& a_mflg, const LevelData<MFBaseIVFAB>& a_Bcoef, const Real a_dx, const int a_order, const int a_weight, const int a_radius){
  MayDay::Warning("JumpBC::JumpBC -- not implemented");
}

JumpBC::~JumpBC(){

}

#include <CD_NamespaceFooter.H>
