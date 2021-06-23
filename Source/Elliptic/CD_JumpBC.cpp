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
#include <CD_VofUtils.H>
#include <CD_LeastSquares.H>
#include <CD_NamespaceHeader.H>

JumpBC::JumpBC(const MFLevelGrid& a_mflg, const LevelData<MFBaseIVFAB>& a_Bcoef, const Real a_dx, const int a_order, const int a_weight, const int a_radius){
  CH_TIME("JumpBC::JumpBC");

  this->define(a_mflg, a_Bcoef, a_dx, a_order, a_weight, a_radius);
}

JumpBC::~JumpBC(){

}

void JumpBC::define(const MFLevelGrid& a_mflg, const LevelData<MFBaseIVFAB>& a_Bcoef, const Real a_dx, const int a_order, const int a_weight, const int a_radius){
  //  MayDay::Warning("JumpBC::define - not implemented");
}

const BaseIVFAB<VoFStencil>& JumpBC::getAvgStencils(const int a_phase, const DataIndex& a_dit) const {
  return m_avgStencils[a_dit].getIVFAB(a_phase);
}

const BaseIVFAB<Real>& JumpBC::getInhomogeneousContribution(const int a_phase, const DataIndex& a_dit) const {
  return m_inhomogeneousContrib[a_dit].getIVFAB(a_phase);
}

const BaseIVFAB<Real>& JumpBC::getHomogeneousContribution(const int a_phase, const DataIndex& a_dit) const {
  return m_homogeneousContrib[a_dit].getIVFAB(a_phase);
}

#include <CD_NamespaceFooter.H>
