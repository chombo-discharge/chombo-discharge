/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzNeumannDomainBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>
#include <CH_Timer.H>

// Our includes
#include <CD_BoxLoops.H>
#include <CD_EBHelmholtzNeumannDomainBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzNeumannDomainBC::EBHelmholtzNeumannDomainBC()
{
  CH_TIME("EBHelmholtzNeumannDomainBC::EBHelmholtzNeumannDomainBC()");

  m_multByBco = true;

  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzNeumannDomainBC::EBHelmholtzNeumannDomainBC(const Real a_DphiDn)
{
  CH_TIME("EBHelmholtzNeumannDomainBC::EBHelmholtzNeumannDomainBC(Real)");

  this->setDphiDn(a_DphiDn);
}

EBHelmholtzNeumannDomainBC::EBHelmholtzNeumannDomainBC(const std::function<Real(const RealVect& a_pos)>& a_DphiDn)
{
  CH_TIME("EBHelmholtzNeumannDomainBC::EBHelmholtzNeumannDomainBC(std::function<Real(RealVect)>)");

  this->setDphiDn(a_DphiDn);
}

EBHelmholtzNeumannDomainBC::~EBHelmholtzNeumannDomainBC()
{
  CH_TIME("EBHelmholtzNeumannDomainBC::~EBHelmholtzNeumannDomainBC()");
}

void
EBHelmholtzNeumannDomainBC::setDphiDn(const Real a_DphiDn)
{
  CH_TIME("EBHelmholtzNeumannDomainBC::setDphiDn(Real)");

  m_multByBco = true;

  m_useConstant = true;
  m_useFunction = false;

  m_constantDphiDn = a_DphiDn;
}

void
EBHelmholtzNeumannDomainBC::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn)
{
  CH_TIME("EBHelmholtzNeumannDomainBC::setDphiDn(std::function<Real(RealVect)>)");

  m_multByBco = true;

  m_useConstant = false;
  m_useFunction = true;

  m_functionDphiDn = a_DphiDn;
}

void
EBHelmholtzNeumannDomainBC::setBxDphiDn(const Real a_BxDphiDn)
{
  CH_TIME("EBHelmholtzNeumannDomainBC::setBxDphiDn(Real)");

  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

void
EBHelmholtzNeumannDomainBC::setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn)
{
  CH_TIME("EBHelmholtzNeumannDomainBC::setBxDphiDn(std::function<Real(RealVect)>)");

  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

void
EBHelmholtzNeumannDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                        const BaseFab<Real>&  a_phi,
                                        const BaseFab<Real>&  a_Bcoef,
                                        const int&            a_dir,
                                        const Side::LoHiSide& a_side,
                                        const DataIndex&      a_dit,
                                        const bool            a_useHomogeneous) const
{
  CH_TIME("EBHelmholtzNeumannDomainBC::getFaceFlux(BaseFab<Real>, BaseFab<Real>, int, DataIndex, bool)");

  CH_assert(a_faceFlux.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);

  const int isign = (a_side == Side::Lo) ? -1 : 1;

  if (a_useHomogeneous) {
    a_faceFlux.setVal(0.0);
  }
  else {
    // Note the conspicuous minus-sign. It's because we assume that a positive DphiDn puts a flux INTO the domain.
    if (m_useConstant) {
      a_faceFlux.setVal(-isign * m_constantDphiDn);
    }
    else if (m_useFunction) {

      // Kernel -- get dphi/dn at spatial position.
      auto kernel = [&](const IntVect& iv) -> void {
        const RealVect pos    = this->getBoundaryPosition(iv, a_dir, a_side);
        const Real     DphiDn = m_functionDphiDn(pos);

        a_faceFlux(iv, m_comp) = -isign * DphiDn;
      };

      // Execute the kernel.
      BoxLoops::loop(a_faceFlux.box(), kernel);
    }

    // Multiply by B-coefficient. We always do this unless the user specifically called setBxDphiDn in which case the input value
    // is already multiplied by the B-coefficient.
    if (m_multByBco) {
      this->multiplyByBcoef(a_faceFlux, a_Bcoef, a_dir, a_side);
    }
  }
}

Real
EBHelmholtzNeumannDomainBC::getFaceFlux(const VolIndex&       a_vof,
                                        const EBCellFAB&      a_phi,
                                        const EBFaceFAB&      a_Bcoef,
                                        const int&            a_dir,
                                        const Side::LoHiSide& a_side,
                                        const DataIndex&      a_dit,
                                        const bool            a_useHomogeneous) const
{
  CH_TIME("EBHelmholtzNeumannDomainBC::getFaceFlux(VolIndex, EBCellFAB, int, Side::LoHiSide, DataIndex, bool)");

  Real centroidFlux = 0.0;

  if (a_useHomogeneous) {
    centroidFlux = 0.0;
  }
  else {
    const int            isign   = (a_side == Side::Lo) ? -1 : 1;
    const IntVect        iv      = a_vof.gridIndex();
    const EBISBox&       ebisbox = m_eblg.getEBISL()[a_dit];
    const ProblemDomain& domain  = m_eblg.getDomain();

    const Vector<FaceIndex> faces = ebisbox.getFaces(a_vof, a_dir, a_side);
    if (faces.size() > 0) {
      if (faces.size() == 1) { // Get an interpolation stencil, using centered differences
        IntVectSet        cfivs;
        const FaceStencil faceSten = EBArith::getInterpStencil(faces[0], cfivs, ebisbox, domain);

        for (int i = 0; i < faceSten.size(); i++) {
          const Real&      weight = faceSten.weight(i);
          const FaceIndex& face   = faceSten.face(i);
          const VolIndex&  curVof = face.getVoF(flip(a_side));

          // Get dphi/dx on the boundary
          Real centeredDphiDn;
          if (m_useConstant) {
            centeredDphiDn = m_constantDphiDn;
          }
          else if (m_useFunction) {
            centeredDphiDn = m_functionDphiDn(this->getBoundaryPosition(curVof.gridIndex(), a_dir, a_side));
          }
          else {
            centeredDphiDn = 0.0;
            MayDay::Error("EBHelmholtzNeumannDomainBC::getFaceFlux - logic bust");
          }

          centroidFlux += weight * centeredDphiDn;
        }

        // Multiply by b-coefficient and aperture.
        const FaceIndex& bndryFace = faces[0];
        const Real       Bco       = m_multByBco ? a_Bcoef(bndryFace, m_comp) : 1.0;
        const Real       area      = ebisbox.areaFrac(bndryFace);

        centroidFlux *= -isign * area * Bco;
      }
      else {
        // If this triggers then the cell edge has multiple phases. It should be possible to handle this somehow.
        MayDay::Error(
          "EBHelmholtzNeumannDomainBC -- boundary face is multivalued and EBHelmholtzNeumannDomainBC does not support that (yet)");
      }
    }
  }

  return centroidFlux;
}

#include <CD_NamespaceFooter.H>
