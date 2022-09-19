/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzDirichletDomainBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>
#include <CH_Timer.H>

// Our includes
#include <CD_BoxLoops.H>
#include <CD_EBHelmholtzDirichletDomainBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC()
{
  CH_TIME("EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC()");

  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC(const Real a_value)
{
  CH_TIME("EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC(Real)");

  this->setValue(a_value);
}

EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC(const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("EBHelmholtzDirichletDomainBC::EBHelmholtzDirichletDomainBC(std::function<(const RealVect a_pos)>)");

  this->setValue(a_value);
}

EBHelmholtzDirichletDomainBC::~EBHelmholtzDirichletDomainBC()
{
  CH_TIME("EBHelmholtzDirichletDomainBC::~EBHelmholtzDirichletDomainBC()");
}

void
EBHelmholtzDirichletDomainBC::setValue(const Real a_value)
{
  CH_TIME("EBHelmholtzDirichletDomainBC::setValue(Real)");

  m_useConstant = true;
  m_useFunction = false;

  m_constantValue = a_value;
}

void
EBHelmholtzDirichletDomainBC::setValue(const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("EBHelmholtzDirichletDomainBC::setValue(std::function<Real(const RealVect)>)");

  m_useConstant = false;
  m_useFunction = true;

  m_functionValue = a_value;
}

void
EBHelmholtzDirichletDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                          const BaseFab<Real>&  a_phi,
                                          const int&            a_dir,
                                          const Side::LoHiSide& a_side,
                                          const DataIndex&      a_dit,
                                          const bool            a_useHomogeneous) const
{
  CH_TIME(
    "EBHelmholtzDirichletDomainBC::getFaceFlux(BaseFab<Real>, BaseFab<Real>, int, Side::LoHiSide, DataIndex, bool)");

  CH_assert(m_useConstant || m_useFunction);
  CH_assert(a_faceFlux.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);

  // TLDR: This fill the regular flux on the domain edge/face.

  const Real ihdx =
    2.0 /
    m_dx; // Spatial step is dx/2 because the finite difference goes between cell center and domain edge (so half a cell).
  const Real sign = (a_side == Side::Lo) ? -1 : 1; // For getting the direction of the derivative correctly.

  std::function<void(const IntVect&)> kernel;

  // Need to figure which kernel we should compute. If we have homogeneous BCs then the boundary value is zero. Likewise, with a non-zero value
  // we have the dphi/dn = (phi-bc_value)/(dx/2) on the low side and (bc_value - phi)/(dx/2) on the high side. So make a switch between homogeneous/inhomogeneous
  // and constant/non-constant values.
  if (a_useHomogeneous) {
    kernel = [&](const IntVect& iv) -> void {
      a_faceFlux(iv, m_comp) = -sign * ihdx * a_phi(iv, m_comp);
    };
  }
  else {                 // Physical BCs, select whether or not we use a constant value of spatially varying value.
    if (m_useConstant) { // Constant value.
      kernel = [&](const IntVect& iv) -> void {
        a_faceFlux(iv, m_comp) = sign * ihdx * (m_constantValue - a_phi(iv, m_comp));
      };
    }
    else if (m_useFunction) { // Spatially varying.
      kernel = [&](const IntVect& iv) -> void {
        const RealVect pos   = this->getBoundaryPosition(iv, a_dir, a_side);
        const Real     value = m_functionValue(pos);

        a_faceFlux(iv, m_comp) = sign * ihdx * (value - a_phi(iv, m_comp));
      };
    }
    else {
      MayDay::Error("EBHelmholtzDirichletDomainBC::getFaceFlux -- logic bust");
    }
  }

  // Execute the kernel.
  BoxLoops::loop(a_faceFlux.box(), kernel);

  // Multiplies by B-coefficient so that a_faceFlux = B*dphi/dn.
  const BaseFab<Real>& Bco = (*m_Bcoef)[a_dit][a_dir].getSingleValuedFAB();
  this->multiplyByBcoef(a_faceFlux, Bco, a_dir, a_side);
}

Real
EBHelmholtzDirichletDomainBC::getFaceFlux(const VolIndex&       a_vof,
                                          const EBCellFAB&      a_phi,
                                          const int&            a_dir,
                                          const Side::LoHiSide& a_side,
                                          const DataIndex&      a_dit,
                                          const bool            a_useHomogeneous) const
{
  CH_TIME("EBHelmholtzDirichletDomainBC::getFaceFlux(VolIndex, EBCellFAB, int, Side::LoHiSide, DataIndex, bool)");

  CH_assert(m_useConstant || m_useFunction);

  const int            isign   = (a_side == Side::Lo) ? -1 : 1;
  const Real           ihdx    = 2.0 / m_dx;
  const IntVect        iv      = a_vof.gridIndex();
  const EBISBox&       ebisbox = m_eblg.getEBISL()[a_dit];
  const ProblemDomain& domain  = m_eblg.getDomain();

  Real centroidFlux = 0.0;

  const Vector<FaceIndex> faces = ebisbox.getFaces(a_vof, a_dir, a_side);

  if (faces.size() > 0) {
    if (faces.size() == 1) { // Get an interpolation stencil, using centered differences
      IntVectSet        cfivs;
      const FaceStencil faceSten = EBArith::getInterpStencil(faces[0], cfivs, ebisbox, domain);

      // Get dphi/dn on the boundary face. Use interpolation from face centers to the face centroid.
      for (int i = 0; i < faceSten.size(); i++) {
        const Real&      weight = faceSten.weight(i);
        const FaceIndex& face   = faceSten.face(i);
        const VolIndex&  curVof = face.getVoF(flip(a_side));
        const IntVect&   curIV  = curVof.gridIndex();

        Real value;
        if (a_useHomogeneous) {
          value = 0.0;
        }
        else {
          if (m_useConstant) {
            value = m_constantValue;
          }
          else if (m_useFunction) {
            value = m_functionValue(this->getBoundaryPosition(curIV, a_dir, a_side));
          }
          else {
            MayDay::Error("EBHelmholtzDirichletDomainBC::getFaceFlux - logic bust");
          }
        }

        const Real centeredFaceFlux = isign * ihdx * (value - a_phi(curVof, m_comp));

        centroidFlux += centeredFaceFlux * weight;
      }

      // Multiply by b-coefficient and aperture.
      const FaceIndex& bndryFace = faces[0];
      const Real       Bco       = (*m_Bcoef)[a_dit][a_dir](bndryFace, m_comp);
      const Real       area      = ebisbox.areaFrac(bndryFace);

      centroidFlux *= Bco * area;
    }
    else {
      MayDay::Error(
        "EBHelmholtzDirichletDomainBC -- boundary face is multivalued and EBHelmholtzDirichletDomainBC does not support this (yet)");
    }
  }

  return centroidFlux;
}

#include <CD_NamespaceFooter.H>
