/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBParticleMeshImplem.H
  @brief  Implementatoin of CD_EBParticleMesh.H
  @author Robert Marskar
*/

#ifndef CD_EBParticleMeshImplem_H
#define CD_EBParticleMeshImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBParticleMesh.H>
#include <CD_BoxLoops.H>
#include <CD_ParticleOps.H>
#include <CD_NamespaceHeader.H>

template <class P, class Ret, Ret (P::*MemberFunc)() const>
void
EBParticleMesh::deposit(EBCellFAB&           a_meshData,
                        const List<P>&       a_particles,
                        const DepositionType a_depositionType,
                        const Real           a_widthScale,
                        const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit");
  if (m_verbose) {
    pout() << "EBParticleMesh::deposit" << endl;
  }

  constexpr int numComp = EBParticleMesh::sanitize<Ret>();

  Real invVol = 1.0;
  for (int dir = 0; dir < SpaceDim; dir++) {
    invVol /= m_dx[dir];
  }

  // These are the default _total_ particle widths.
  const RealVect cicWidth = 1 * a_widthScale * RealVect::Unit;
  const RealVect tscWidth = 2 * a_widthScale * RealVect::Unit;

  const Box cicBox(-1 * IntVect::Unit, 1 * IntVect::Unit);
  const Box tscBox(-2 * IntVect::Unit, 2 * IntVect::Unit);

  Real ngpFactor = invVol;
  Real cicFactor = invVol;
  Real tscFactor = invVol;

  // Various scaling factors since particle deposition functions require different normalizations.
  for (int dir = 0; dir < SpaceDim; dir++) {
    cicFactor *= 1.0 / cicWidth[dir];
    tscFactor *= 2.0 / tscWidth[dir];
  }

  // Switch between deposition types.
  switch (a_depositionType) {
  case DepositionType::NGP: {
    for (ListIterator<P> lit(a_particles); lit.ok(); ++lit) {
      const Ret w = (lit().*MemberFunc)();

      this->depositParticleNGP(a_meshData, lit().position(), ngpFactor, &w, numComp);
    }

    break;
  }
  case DepositionType::CIC: {
    for (ListIterator<P> lit(a_particles); lit.ok(); ++lit) {
      const Ret w = (lit().*MemberFunc)();

      this->depositParticleCIC(a_meshData, lit().position(), cicWidth, cicBox, cicFactor, &w, numComp, a_forceIrregNGP);
    }

    break;
  }
  case DepositionType::TSC: {
    for (ListIterator<P> lit(a_particles); lit.ok(); ++lit) {
      const Ret w = (lit().*MemberFunc)();

      this->depositParticleTSC(a_meshData, lit().position(), tscWidth, tscBox, tscFactor, &w, numComp, a_forceIrregNGP);
    }

    break;
  }
  default: {
    MayDay::Abort("EBParticleMesh::deposit -- logic bust");

    break;
  }
  }
}

template <class P, class Ret, Ret (P::*MemberFunc)()>
void
EBParticleMesh::interpolate(List<P>&             a_particleList,
                            const EBCellFAB&     a_meshData,
                            const DepositionType a_interpType,
                            const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::interpolate(Real)");

  static_assert(std::is_same<Ret, Real&>::value || std::is_same<Ret, RealVect&>::value,
                "Ret should be Real& or RealVect&");

  constexpr int numComp = std::is_same<Ret, Real&>::value ? 1 : SpaceDim;

  CH_assert(a_meshData.nComp() == numComp);

  const Interval variables = Interval(0, numComp - 1);

  // validBox is the domain in which we can use the specified interpolation scheme. E.g., TSC fails nears the domain boundaries
  // since we don't have enough ghost cells to interpolate from.
  //
  // gatherBox is the box from which we gather the field.
  //
  Box validBox;
  Box gatherBox;

  switch (a_interpType) {
  case DepositionType::NGP: {
    validBox  = m_domain.domainBox();
    gatherBox = Box(IntVect::Zero, IntVect::Zero);

    break;
  }
  case DepositionType::CIC: {
    validBox  = grow(m_domain.domainBox(), -1);
    gatherBox = Box(IntVect::Zero, IntVect::Unit);

    break;
  }
  case DepositionType::TSC: {
    validBox  = grow(m_domain.domainBox(), -2);
    gatherBox = Box(IntVect::Zero, 2 * IntVect::Unit);

    break;
  }
  default: {
    MayDay::Error("EBParticleMesh::interpolate - logic bust");
  }
  }

  /*!
    @brief Local helper class for turning a Real or a RealVect into a pointer to the beginning of the object. I'm doing this so that I can
    use a single version of interpolateParticle that does not require template specialization. 
  */
  struct GetPointer
  {
    /*!
      @brief Return as pointer
    */
    Real*
    operator()(Real& r) const
    {
      return &r;
    }

    /*!
      @brief Return pointer to start of RealVect object
    */
    auto
    operator()(RealVect& r) const
    {
      return r.dataPtr();
    }
  };

  switch (a_interpType) {
  case DepositionType::NGP: {
    for (ListIterator<P> lit(a_particleList); lit; ++lit) {
      auto&       p = lit();
      auto&&      w = (p.*MemberFunc)();
      const auto& x = p.position();

      this->interpolateParticleNGP(GetPointer{}(w), a_meshData, x, numComp);
    }

    break;
  }
  case DepositionType::CIC: {
    for (ListIterator<P> lit(a_particleList); lit; ++lit) {
      auto&       p = lit();
      auto&&      w = (p.*MemberFunc)();
      const auto& x = p.position();

      this->interpolateParticleCIC(GetPointer{}(w), a_meshData, validBox, gatherBox, x, numComp, a_forceIrregNGP);
    }

    break;
  }
  case DepositionType::TSC: {
    for (ListIterator<P> lit(a_particleList); lit; ++lit) {
      auto&       p = lit();
      auto&&      w = (p.*MemberFunc)();
      const auto& x = p.position();

      this->interpolateParticleTSC(GetPointer{}(w), a_meshData, validBox, gatherBox, x, numComp, a_forceIrregNGP);
    }

    break;
  }
  default: {
    MayDay::Abort("EBParticleMesh::interpolate -- logic bust");

    break;
  }
  }
}

inline void
EBParticleMesh::depositParticleNGP(EBCellFAB&      a_rho,
                                   const RealVect& a_position,
                                   const Real&     a_volumeFactor,
                                   const Real*     a_strength,
                                   const int&      a_numComp) const noexcept
{
  CH_TIME("EBParticleMesh::depositParticleNGP");

  const IntVect particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);

  CH_assert(m_region.contains(particleIV));

  FArrayBox& rho = a_rho.getFArrayBox();

  for (int comp = 0; comp < a_numComp; comp++) {
    rho(particleIV, comp) += a_strength[comp] * a_volumeFactor;
  }
}

inline void
EBParticleMesh::depositParticleCIC(EBCellFAB&      a_rho,
                                   const RealVect& a_position,
                                   const RealVect& a_particleWidth,
                                   const Box&      a_particleBox,
                                   const Real&     a_volumeFactor,
                                   const Real*     a_strength,
                                   const int&      a_numComp,
                                   const bool      a_forceIrregNGP) const noexcept
{
  CH_TIME("EBParticleMesh::depositParticleCIC");

  const IntVect particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);

  CH_assert(m_region.contains(particleIV));

  FArrayBox& rho = a_rho.getFArrayBox();

  if (m_ebisbox.isIrregular(particleIV) && a_forceIrregNGP) {
    for (int comp = 0; comp < a_numComp; comp++) {
      rho(particleIV, comp) += a_strength[comp] * a_volumeFactor;
    }
  }
  else {
    auto cicKernel = [&](const IntVect& iv) -> void {
      Real weight = a_volumeFactor;

      // Below, we've transformed the coordinates so that the cell center of cell iv is at the origin.
      for (int dir = 0; dir < SpaceDim; dir++) {
        const Real a = (m_probLo[dir] - a_position[dir]) / m_dx[dir] + iv[dir];
        const Real b = a + 1.0;
        const Real L = 0.5 * a_particleWidth[dir];

        weight *= std::max(0.0, std::min(b, L) - std::max(a, -L));
      }

      for (int comp = 0; comp < a_numComp; comp++) {
        rho(iv, comp) += weight * a_strength[comp];
      }
    };

    const Box particleBox = a_particleBox + particleIV;

    BoxLoops::loop(particleBox, cicKernel);
  }
}

inline void
EBParticleMesh::depositParticleTSC(EBCellFAB&      a_rho,
                                   const RealVect& a_position,
                                   const RealVect& a_particleWidth,
                                   const Box&      a_particleBox,
                                   const Real&     a_volumeFactor,
                                   const Real*     a_strength,
                                   const int&      a_numComp,
                                   const bool      a_forceIrregNGP) const noexcept
{
  CH_TIME("EBParticleMesh::depositParticleTSC");

  const IntVect particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);

  CH_assert(m_region.contains(particleIV));

  FArrayBox& rho = a_rho.getFArrayBox();

  auto F = [](const Real x) -> Real {
    const Real a = (x <= 0.0) ? 2 * x * (1 + x) : 2 * x * (1 - x);

    return (std::abs(x) <= 0.5) ? a : 0.0;
  };

  if (m_ebisbox.isIrregular(particleIV) && a_forceIrregNGP) {
    for (int comp = 0; comp < a_numComp; comp++) {
      rho(particleIV, comp) += a_strength[comp] * a_volumeFactor;
    }
  }
  else {
    auto tscKernel = [&](const IntVect& iv) -> void {
      Real weight = a_volumeFactor;

      for (int dir = 0; dir < SpaceDim; dir++) {
        const Real a = (m_probLo[dir] - a_position[dir]) / m_dx[dir] + iv[dir];
        const Real b = a + 1.0;
        const Real L = a_particleWidth[dir];

        const Real alpha  = std::max(a, -0.5 * L);
        const Real beta   = std::min(b, +0.5 * L);
        const Real factor = (alpha < beta) ? 1.0 : 0.0;

        weight *= factor * (beta - alpha) - (beta * std::abs(beta) - alpha * std::abs(alpha)) / L;
      }

      for (int comp = 0; comp < a_numComp; comp++) {
        rho(iv, comp) += weight * a_strength[comp];
      }
    };

    const Box particleBox = a_particleBox + particleIV;

    BoxLoops::loop(particleBox, tscKernel);
  }
}

inline void
EBParticleMesh::interpolateParticleNGP(Real*            a_particleField,
                                       const EBCellFAB& a_meshData,
                                       const RealVect&  a_position,
                                       const int&       a_numComp) const noexcept
{
  CH_TIME("EBParticleMesh::interpolateParticleNGP");

  const IntVect    particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);
  const FArrayBox& meshData   = a_meshData.getFArrayBox();
  const Real       factor     = m_ebisbox.isCovered(particleIV) ? 0.0 : 1.0;

  for (int comp = 0; comp < a_numComp; comp++) {
    a_particleField[comp] = factor * meshData(particleIV, comp);
  }
}

inline void
EBParticleMesh::interpolateParticleCIC(Real*            a_particleField,
                                       const EBCellFAB& a_meshData,
                                       const Box&       a_validBox,
                                       const Box&       a_gatherBox,
                                       const RealVect&  a_position,
                                       const int&       a_numComp,
                                       const bool       a_forceIrregNGP) const noexcept
{
  CH_TIME("EBParticleMesh::interpolateParticleCIC");

  // particleIV is the cell index where the particle center lies.
  // loIndex is the cell index of the lower-left corner of the particle.
  const IntVect particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);
  const IntVect loIndex    = ParticleOps::getParticleCellIndex(a_position - 0.5 * m_dx, m_probLo, m_dx);

  CH_assert(m_region.contains(particleIV));

  // Get regular data.
  const FArrayBox& meshData = a_meshData.getFArrayBox();

  // Reset the particle field.
  for (int comp = 0; comp < a_numComp; comp++) {
    a_particleField[comp] = 0.0;
  }

  // NGP interpolation if called for it, or particle abuts the domain boundary.
  if ((m_ebisbox.isIrregular(particleIV) && a_forceIrregNGP) || !(a_validBox.contains(particleIV))) {
    for (int comp = 0; comp < a_numComp; comp++) {
      a_particleField[comp] = meshData(particleIV, comp);
    }
  }
  else if (!(m_ebisbox.isCovered(particleIV))) {
    auto cicKernel = [&](const IntVect& iv) -> void {
      Real weight = 1.0;

      const RealVect L = (m_probLo - a_position) / m_dx + (RealVect(iv) + 0.5 * RealVect::Unit);

      for (int dir = 0; dir < SpaceDim; dir++) {
        weight *= (1. - std::abs(L[dir]));
      }

      for (int comp = 0; comp < a_numComp; comp++) {
        a_particleField[comp] += weight * meshData(iv, comp);
      }
    };

    const Box gatherBox = a_gatherBox + loIndex;

    BoxLoops::loop(gatherBox, cicKernel);
  }
}

inline void
EBParticleMesh::interpolateParticleTSC(Real*            a_particleField,
                                       const EBCellFAB& a_meshData,
                                       const Box&       a_validBox,
                                       const Box&       a_gatherBox,
                                       const RealVect&  a_position,
                                       const int&       a_numComp,
                                       const bool       a_forceIrregNGP) const noexcept
{
  CH_TIME("EBParticleMesh::interpolateParticleTSC");

  // particleIV is the cell index where the particle center lies.
  // loIndex is the cell index of the lower-left corner of the particle.
  const IntVect particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);
  const IntVect loIndex    = ParticleOps::getParticleCellIndex(a_position - m_dx, m_probLo, m_dx);

  CH_assert(m_region.contains(particleIV));

  // Get regular data.
  const FArrayBox& meshData = a_meshData.getFArrayBox();

  // Reset the particle field.
  for (int comp = 0; comp < a_numComp; comp++) {
    a_particleField[comp] = 0.0;
  }

  // NGP interpolation if called for it, or particle abuts the domain boundary.
  if ((m_ebisbox.isIrregular(particleIV) && a_forceIrregNGP) || !(a_validBox.contains(particleIV))) {
    for (int comp = 0; comp < a_numComp; comp++) {
      a_particleField[comp] = meshData(particleIV, comp);
    }
  }
  else if (!(m_ebisbox.isCovered(particleIV))) {
    auto tscKernel = [&](const IntVect& iv) -> void {
      Real weight = 1.0;

      const RealVect L = (m_probLo - a_position) / m_dx + (RealVect(iv) + 0.5 * RealVect::Unit);

      for (int dir = 0; dir < SpaceDim; dir++) {
        const Real& l = std::abs(L[dir]);

        if (l < 0.5) {
          weight *= 0.75 - l * l;
        }
        else {
          weight *= 0.5 * (1.5 - l) * (1.5 - l);
        }
      }

      for (int comp = 0; comp < a_numComp; comp++) {
        a_particleField[comp] += weight * meshData(iv, comp);
      }
    };

    const Box gatherBox = a_gatherBox + loIndex;

    BoxLoops::loop(gatherBox, tscKernel);
  }
}

#include <CD_NamespaceFooter.H>

#endif
