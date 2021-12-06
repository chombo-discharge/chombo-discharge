/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_AdvectionDiffusionSpecies.cpp
  @brief  Implementation of CD_AdvectionDiffusionSpecies.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <PolyGeom.H>

// Our includes
#include <CD_AdvectionDiffusionSpecies.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::AdvectionDiffusion;

AdvectionDiffusionSpecies::AdvectionDiffusionSpecies(){

  ParmParse pp("AdvectionDiffusion");

  m_chargeNumber = 0;
  m_name   = "scalar species";

  // Parse options
  Vector<Real> v;
  pp.get   ("diffusion",      m_isDiffusive);
  pp.get   ("advection",      m_isMobile);
  pp.get   ("blob_amplitude", m_blobAmplitude);
  pp.get   ("blob_radius",    m_blobRadius);
  pp.getarr("blob_center",    v, 0, SpaceDim); m_blobCenter = RealVect(D_DECL(v[0], v[1], v[2]));
}

AdvectionDiffusionSpecies::~AdvectionDiffusionSpecies(){

}

Real AdvectionDiffusionSpecies::initialData(const RealVect a_pos, const Real a_time) const{

  // TLDR: Make a super-Gaussian blob centered on m_blobCenter and radius/amplitude
  //       m_blobRadius/m_blobAmplitude
  
  const RealVect d  = a_pos - m_blobCenter;
  const Real     d2 = PolyGeom::dot(d,d);
  const Real     r2 = m_blobRadius*m_blobRadius;

  return m_blobAmplitude*exp(-0.5*d2*d2/(r2*r2));
}

#include <CD_NamespaceFooter.H>
