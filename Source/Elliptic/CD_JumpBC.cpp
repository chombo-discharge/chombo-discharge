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

const BaseIVFAB<Real>& JumpBC::getBndryPhi(const int a_phase, const DataIndex& a_dit) const {
  return m_boundaryPhi[a_dit].getIVFAB(a_phase);
}

bool JumpBC::getLeastSquaresBoundaryGradStencil(std::pair<Real, VoFStencil>& a_stencil,
						const VolIndex&              a_vof,
						const VofUtils::Neighborhood a_neighborhood,
						const DataIndex&             a_dit,
						const int                    a_order,
						const int                    a_phase) const {
  bool foundStencil = false;

  const EBLevelGrid& eblg = m_mflg.getEBLevelGrid(a_phase);
  const EBISBox& ebisbox  = eblg.getEBISL()[a_dit];
  const RealVect normal   = ebisbox.normal(a_vof);  
    
  const VoFStencil gradientStencil = LeastSquares::getBndryGradSten(a_vof, a_neighborhood, ebisbox, m_dx, a_order, m_weight, a_order);

  if(gradientStencil.size() > 0 && normal != RealVect::Zero){
    
    const VoFStencil DphiDnStencil =  LeastSquares::projectGradSten(gradientStencil, -normal);
    const Real boundaryWeight      = -LeastSquares::sumAllWeights(DphiDnStencil);

    a_stencil = std::make_pair(boundaryWeight, DphiDnStencil);

    foundStencil = true;
  }

  return foundStencil;
}

#include <CD_NamespaceFooter.H>
