/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_JumpBCImplem.H
  @brief  Implementation of CD_JumpBC.H
  @author Robert Marskar
  @todo   The denominator should be precomputed for efficiency reasons. 
  @todo   Inhomogeneous contribution is missing. 
*/

#ifndef CD_JumpBCImplem_H
#define CD_JumpBCImplem_H

// Chombo includes
#include <NeighborIterator.H>

// Our includes
#include <CD_JumpBC.H>
#include <CD_NamespaceHeader.H>

inline
bool JumpBC::isStencilValidCF(const VoFStencil& a_stencil, const DataIndex& a_dit) const {
  const DisjointBoxLayout& dbl = m_mflg.getGrids();
  const ProblemDomain& domain  = m_mflg.getDomain();

  // Construct boxes that contain all the valid cells for this stencil. 
  std::vector<Box> validBoxes;
  
  Box curBox = dbl[a_dit];
  curBox.grow(m_ghostCF);
  curBox &= domain;
  validBoxes.emplace_back(curBox);

  NeighborIterator nit(dbl);
  for (nit.begin(a_dit); nit.ok(); ++nit){
    Box neighBox = dbl[nit()];
    neighBox.grow(m_ghostCF);
    neighBox &= domain;

    validBoxes.emplace_back(neighBox);
  }

  // Now check that the stencil if the stencil. We set valid = false
  // if any of the stencil points reaches out of the ghosted boxes. 
  bool valid = true;
  
  for (int i = 0; i < a_stencil.size(); i++){
    const VolIndex& vof = a_stencil.vof(i);
    const IntVect& iv   = vof.gridIndex();
    
    bool insideOneBox = false;
    for (const auto& b : validBoxes){
      if(b.contains(iv)) insideOneBox = true;
    }

    if(!insideOneBox) {
      valid = false;
      break;
    }
  }

  return valid;
}

inline
Real JumpBC::applyStencil(const VoFStencil& a_stencil, const EBCellFAB& a_phi) const {
  Real contrib = 0.0;
  for (int i = 0; i < a_stencil.size(); i++){
    contrib += a_stencil.weight(i)*a_phi(a_stencil.vof(i), m_comp);
  }

  return contrib;
}

inline
void JumpBC::matchBC(const MFCellFAB&       a_phi,
		     const BaseIVFAB<Real>& a_jump,
		     const bool             a_homogeneousPhysBC,
		     const DataIndex&       a_dit) const {
  CH_assert(m_multiPhase);
  
  constexpr int vofComp     = 0;
  constexpr int firstPhase  = 0;
  constexpr int secondPhase = 1;

  const EBCellFAB& phiPhase0 = a_phi.getPhase(firstPhase );
  const EBCellFAB& phiPhase1 = a_phi.getPhase(secondPhase);

  const EBISBox& ebisBoxPhase0 = phiPhase0.getEBISBox();
  const EBISBox& ebisBoxPhase1 = phiPhase1.getEBISBox();

  BaseIVFAB<Real>& bndryPhiPhase0 = m_boundaryPhi[a_dit].getIVFAB(firstPhase );
  BaseIVFAB<Real>& bndryPhiPhase1 = m_boundaryPhi[a_dit].getIVFAB(secondPhase);
  
  for (IVSIterator ivsIt(m_ivs[a_dit]); ivsIt.ok(); ++ivsIt){
    const IntVect iv    = ivsIt();
    const VolIndex vof0 = VolIndex(iv, vofComp);

    const VoFStencil& derivStenPhase0 = m_avgStencils[a_dit].getIVFAB(firstPhase) (vof0, vofComp);
    const VoFStencil& derivStenPhase1 = m_avgStencils[a_dit].getIVFAB(secondPhase)(vof0, vofComp);

    const Real& weightPhase0 = m_avgWeights[a_dit].getIVFAB(firstPhase) (vof0, vofComp);
    const Real& weightPhase1 = m_avgWeights[a_dit].getIVFAB(secondPhase)(vof0, vofComp);

    // Compute the average jump.
    Real jump = 0.0; // Jump not computed (yet)
    // Vector<VolIndex> vofs = ebisBoxPhase0.getVofs(blablabla)

    // Do the matching 
    const Real denom         = 1./(weightPhase0 + weightPhase1); 
    const Real contribPhase0 = this->applyStencil(derivStenPhase0, phiPhase0);
    const Real contribPhase1 = this->applyStencil(derivStenPhase1, phiPhase1);
    const Real phiBndry      = denom*(jump - contribPhase0 - contribPhase1);

    // Copy the result to the individual phase data holders
    Vector<VolIndex> vofsPhase0 = ebisBoxPhase0.getVoFs(iv);
    Vector<VolIndex> vofsPhase1 = ebisBoxPhase1.getVoFs(iv);

    for (const auto& v : vofsPhase0.stdVector()) bndryPhiPhase0(v, m_comp) = phiBndry;
    for (const auto& v : vofsPhase1.stdVector()) bndryPhiPhase1(v, m_comp) = phiBndry;
  }
}

#include <CD_NamespaceFooter.H>

#endif
