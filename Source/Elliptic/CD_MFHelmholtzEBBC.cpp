/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzEBBC.cpp
  @brief  Implementation of CD_MFHelmholtzEBBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzEBBC::MFHelmholtzEBBC(const int a_phase, const RefCountedPtr<JumpBC>& a_jumpBC){
  m_phase  = a_phase;
  m_jumpBC = a_jumpBC;
}

MFHelmholtzEBBC::~MFHelmholtzEBBC(){

}
  
void MFHelmholtzEBBC::applyEBFlux(VoFIterator&       a_vofit,
				  EBCellFAB&         a_Lphi,
				  const EBCellFAB&   a_phi,
				  const DataIndex&   a_dit,
				  const Real&        a_beta,
				  const bool&        a_homogeneousPhysBC) const {
  
  VoFIterator& singlePhaseVofs = m_jumpBC->getSinglePhaseVofs(m_phase, a_dit);
  VoFIterator& multiPhaseVofs  = m_jumpBC->getMultiPhaseVofs (m_phase, a_dit);

  this->applyEBFlux(singlePhaseVofs, multiPhaseVofs, a_Lphi, a_phi, a_dit, a_beta, a_homogeneousPhysBC);
}

#include <CD_NamespaceFooter.H>
