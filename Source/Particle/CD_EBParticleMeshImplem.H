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

// Std incldues
#include <type_traits>

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBParticleMesh.H>
#include <CD_BoxLoops.H>
#include <CD_ParticleOps.H>
#include <CD_NamespaceHeader.H>

inline Real*
to_bar_arg(Real& r)
{
  return &r;
}
inline auto
to_bar_arg(RealVect& v) -> decltype(v.dataPtr())
{
  return v.dataPtr();
}

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

  Box validBox = m_domain.domainBox();

  switch (a_interpType) {
  case DepositionType::NGP: {
    validBox = m_domain.domainBox();

    break;
  }

  case DepositionType::CIC: {
    validBox = grow(validBox, -1);

    break;
  }
  case DepositionType::TSC: {
    validBox = grow(validBox, -2);

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

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    P&              p = lit();
    const RealVect& x = p.position();

    auto&& w = (p.*MemberFunc)();

    this->interpolateParticle(GetPointer{}(w),
                              a_meshData,
                              validBox,
                              m_probLo,
                              m_dx,
                              x,
                              variables,
                              a_interpType,
                              a_forceIrregNGP);
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

  CH_assert(m_region.contains(iv));

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

  CH_assert(m_region.contains(iv));

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

  CH_assert(m_region.contains(iv));

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
EBParticleMesh::interpolateParticle(Real*                a_particleField,
                                    const EBCellFAB&     a_meshField,
                                    const Box&           a_validBox,
                                    const RealVect&      a_probLo,
                                    const RealVect&      a_dx,
                                    const RealVect&      a_position,
                                    const Interval&      a_interval,
                                    const DepositionType a_interpType,
                                    const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::interpolateParticle");

  const int startComp = a_interval.begin();
  const int endComp   = a_interval.end();

  CH_assert(a_meshField.nComp() >= endComp - 1);

  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.

  // Nifty lambda for converting RealVect position to IntVect (lower-left corner). Note that the input vector
  // is the displacement from a_probLo and not the physical position.
  const auto& dx = a_dx;

  auto getCellIndex = [&dx](const RealVect a_rv) -> IntVect {
    return IntVect(D_DECL(std::floor(a_rv[0] / dx[0]), std::floor(a_rv[1] / dx[1]), std::floor(a_rv[2] / dx[2])));
  };

  // Get grid cell corresponding to this particle.
  const IntVect particleIndex = getCellIndex(a_position - a_probLo);

  // Assertion -- particle must live on this patch.
  CH_assert(m_region.contains(particleIndex));

  // Nifty little lambda for computing the displacement between the particle and a grid cell.
  auto particleDisplacement = [&](const IntVect& iv) -> RealVect {
    return (a_probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * a_dx - a_position) / a_dx;
  };

  // Get regular data.
  const FArrayBox& meshField = a_meshField.getFArrayBox();

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side.
  if ((m_ebisbox.isIrregular(particleIndex) && a_forceIrregNGP) || !(a_validBox.contains(particleIndex))) {
    for (int comp = startComp; comp <= endComp; comp++) {
      a_particleField[comp] = meshField(particleIndex, comp);
    }
  }
  else if (m_ebisbox.isCovered(particleIndex)) { // Need to set to something.
    for (int comp = startComp; comp <= endComp; comp++) {
      a_particleField[comp] = 0.0;
    }
  }
  else {
    switch (a_interpType) {
    case DepositionType::NGP: // Hook for nearest-grid-point interpolation.
    {
      for (int comp = startComp; comp <= endComp; comp++) {
        a_particleField[comp] = meshField(particleIndex, comp);
      }

      break;
    }
    case DepositionType::CIC: // Hook for cloud-in-cell interpolation.
    {
      // Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
      // the box corresponding to the cloud. This region contains the cells where we gather the field.
      const IntVect loIndex = getCellIndex(a_position - a_probLo - 0.5 * a_dx);
      const Box     cicBox  = Box(loIndex, loIndex + IntVect::Unit);

      // This is the cloud-in-cell interpolation kernel. It computes the weight from the current cell and
      // gathers the force.
      auto cicKernel = [&](const IntVect& iv) -> void {
        Real weight = 1.0;

        // L is the distance between the grid cell and the particle.
        const RealVect L = particleDisplacement(iv);
        for (int dir = 0; dir < SpaceDim; dir++) {
          weight *= (1. - std::abs(L[dir]));
        }

        for (int comp = startComp; comp <= endComp; comp++) {
          a_particleField[comp] += weight * meshField(iv, comp);
        }
      };

      // Run kernel and gather field from mesh.
      for (int comp = startComp; comp <= endComp; comp++) {
        a_particleField[comp] = 0.0;
      }

      BoxLoops::loop(cicBox, cicKernel);

      break;
    }
    case DepositionType::TSC: // Hook for triangle-shaped cloud interpolation.
    {
      // Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
      // the box corresponding to the cloud. This region contains the cells where we gather the field.
      const IntVect loIndex = getCellIndex(a_position - a_probLo - 1.0 * a_dx);
      const Box     tscBox  = Box(loIndex, loIndex + 2 * IntVect::Unit);

      // This is the triangle-shaped-cloud interpolation kernel. It computes the weight from the current cell and
      // gather the force.
      auto tscKernel = [&](const IntVect& iv) -> void {
        Real weight = 1.0;

        // L is the distance between the grid cell and the particle.
        const RealVect L = particleDisplacement(iv);
        for (int dir = 0; dir < SpaceDim; dir++) {
          const Real& l = std::abs(L[dir]);

          if (l < 0.5) {
            weight *= 0.75 - l * l;
          }
          else {
            weight *= 0.5 * (1.5 - l) * (1.5 - l);
          }
        }

        for (int comp = startComp; comp <= endComp; comp++) {
          a_particleField[comp] += weight * meshField(iv, comp);
        }
      };

      // Run kernel and gather field from mesh.
      for (int comp = startComp; comp <= endComp; comp++) {
        a_particleField[comp] = 0.0;
      }

      BoxLoops::loop(tscBox, tscKernel);

      break;
    }
    default:
      MayDay::Error("EBParticleMesh::interpolateParticle(RealVect) - Invalid interpolation type requested.");
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
