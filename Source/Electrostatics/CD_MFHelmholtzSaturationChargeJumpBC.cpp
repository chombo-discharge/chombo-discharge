/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_MFHelmholtzSaturationChargeJumpBC.cpp
 * @brief  Implementation of CD_MFHelmholtzSaturationChargeJumpBC.H
 * @author Robert Marskar
 */

// Our includes
#include <CD_MFHelmholtzSaturationChargeJumpBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzSaturationChargeJumpBC::MFHelmholtzSaturationChargeJumpBC(const phase::which_phase     a_phase,
                                                                     const Location::Cell         a_dataLocation,
                                                                     const MFLevelGrid&           a_mflg,
                                                                     const BcoefPtr&              a_Bcoef,
                                                                     const AmrMask&               a_validCells,
                                                                     const Real                   a_dx,
                                                                     const int                    a_order,
                                                                     const int                    a_weight,
                                                                     const int                    a_radius,
                                                                     const int                    a_ghostCF,
                                                                     const IntVect                a_ghostPhi,
                                                                     const RealVect&              a_probLo,
                                                                     const RefCountedPtr<BaseIF>& a_interface)
  : MFHelmholtzJumpBC(a_dataLocation,
                      a_mflg,
                      a_Bcoef,
                      a_validCells,
                      a_dx,
                      a_order,
                      a_weight,
                      a_radius,
                      a_ghostCF,
                      a_ghostPhi,
                      a_probLo,
                      a_interface),
    m_phase(a_phase)
{
  CH_TIME("MFHelmholtzSaturationChargeJumpBC::MFHelmholtzSaturationChargeJumpBC");
}

MFHelmholtzSaturationChargeJumpBC::~MFHelmholtzSaturationChargeJumpBC()
{
  CH_TIME("MFHelmholtzSaturationChargeJumpBC::~MFHelmholtzSaturationChargeJumpBC");
}

void
MFHelmholtzSaturationChargeJumpBC::matchBC(BaseIVFAB<Real>& a_jump,
                                           const MFCellFAB& a_phi,
                                           const bool /*a_homogeneousPhysBC*/,
                                           const DataIndex& a_dit) const
{
  CH_assert(m_multiPhase);

  // TLDR: Normally, the boundary potential would be computed from dphi/dn1 + dphi/dn2 = sigma and since dphi/dn can be
  // represented as
  //
  //          dphi/dn = wb*phiB + sum[wi * phi(i)]
  //
  //       we would have an equation which could be solved for phiB. That value would be used by MFHelmholtzEBBC for
  //       putting in the finite volume fluxes.
  //
  //       This function overrides that functionality and computes phiB by assuming that dphi/dn=0 on one of the phases
  //       (m_phase), and rather that a_jump is the free parameter. In this form we have
  //
  //           phiB = -sum[wi * phi(i)]/wb.
  //
  //       We happen to know that the stencils were stored in "flux-form" by multiplying by the b-coefficient AND the
  //       denominator 1/(wq*bq + wp*bp). So, what we have stored in the weights are actually b*wB and for the stencils
  //       we have stored b*wi/(bq*wq + bp*wp). So, both the weight and the stencil were multiplied by the b-coefficient
  //       so we don't have to worry about the coefficient when we divide that out. But, the term (bq*wq + bp*wp) has
  //       been multiplied into the stencil weights so we need to divide that back out.

  constexpr int vofComp     = 0;
  constexpr int firstPhase  = 0;
  constexpr int secondPhase = 1;

  const EBCellFAB& phiPhase0 = a_phi.getPhase(firstPhase);
  const EBCellFAB& phiPhase1 = a_phi.getPhase(secondPhase);

  const EBISBox& ebisBoxPhase0 = phiPhase0.getEBISBox();
  const EBISBox& ebisBoxPhase1 = phiPhase1.getEBISBox();

  BaseIVFAB<Real>& bndryPhiPhase0 = m_boundaryPhi[a_dit].getIVFAB(firstPhase);
  BaseIVFAB<Real>& bndryPhiPhase1 = m_boundaryPhi[a_dit].getIVFAB(secondPhase);

  for (IVSIterator ivsIt(m_ivs[a_dit]); ivsIt.ok(); ++ivsIt) {
    const IntVect iv = ivsIt();

    const VolIndex vof0 = VolIndex(iv, vofComp);

    const VoFStencil& derivStenPhase0 = m_avgStencils[a_dit].getIVFAB(firstPhase)(vof0, vofComp);
    const VoFStencil& derivStenPhase1 = m_avgStencils[a_dit].getIVFAB(secondPhase)(vof0, vofComp);

    const Real& denomPhase0 = m_denom[a_dit].getIVFAB(firstPhase)(vof0, vofComp);

    const Real& weightPhase0 = m_avgWeights[a_dit].getIVFAB(firstPhase)(vof0, vofComp);
    const Real& weightPhase1 = m_avgWeights[a_dit].getIVFAB(secondPhase)(vof0, vofComp);

    const Real contribPhase0 = this->applyStencil(derivStenPhase0, phiPhase0);
    const Real contribPhase1 = this->applyStencil(derivStenPhase1, phiPhase1);

    // The side with zero normal derivative (z) and the other side (o). The stencils were multiplied by the denominator,
    // so divide it back out. The denominator is zero only if neither side has a stencil, and then there is no flux.
    const bool zIsPhase0 = (m_phase == phase::gas);

    const Real weightZ  = zIsPhase0 ? weightPhase0 : weightPhase1;
    const Real weightO  = zIsPhase0 ? weightPhase1 : weightPhase0;
    const Real contribZ = zIsPhase0 ? contribPhase0 : contribPhase1;
    const Real contribO = zIsPhase0 ? contribPhase1 : contribPhase0;

    Real phiBndry = 0.0;
    Real jump     = 0.0;

    if (denomPhase0 != 0.0) {
      if (weightZ > 0.0) {
        // Zero normal derivative on side z: phiB = -sum(wi*phi(i))/wb, with the weights divided by bq*wq. BC is
        // b*dphi/dn = sigma on the other side, so compute the jump from that.
        phiBndry = contribZ / denomPhase0; // Divide because we multiplied the stencil by bq/(wq*bq + wp*bp)
        phiBndry *= -1. / weightZ;         // Divide by boundary weight, which is actually bq*wq.

        jump += weightO * phiBndry;
        jump += contribO / denomPhase0;
      }
      else {
        // Side z has no stencil. Impose zero normal derivative on the other side instead, so that there is no flux
        // through the interface and no surface charge.
        phiBndry = contribO / denomPhase0;
        phiBndry *= -1. / weightO;
      }
    }

    // Copy the result to the individual phase data holders
    Vector<VolIndex> vofsPhase0 = ebisBoxPhase0.getVoFs(iv);
    Vector<VolIndex> vofsPhase1 = ebisBoxPhase1.getVoFs(iv);

    for (const auto& v : vofsPhase0.stdVector()) {
      bndryPhiPhase0(v, m_comp) = phiBndry;
    }
    for (const auto& v : vofsPhase1.stdVector()) {
      bndryPhiPhase1(v, m_comp) = phiBndry;
    }

    for (const auto& v : vofsPhase0.stdVector()) {
      a_jump(v, m_comp) = jump;
    }
  }
}

#include <CD_NamespaceFooter.H>
