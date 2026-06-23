/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_EBParticleMeshSoA.H
  @brief  Dev staging: ParticleSoA port of EBParticleMesh (single-patch deposit/interpolate).
  @author Robert Marskar
  @details
  This is the Dev-folder port of Source/Particle/CD_EBParticleMesh.H onto the ParticleSoA
  leaf. The per-particle EB/cut-cell kernels (the depositParticleXXX / interpolateParticleXXX
  helpers) are adapted from production -- they are storage-agnostic (position + Real* field +
  numComp) and gain only a destination-component offset. Only the public loops change:

    - storage: List<P> -> ParticleSoA<P, Traits>;
    - field selector: pointer-to-member-FUNCTION (&P::weight) -> a variadic pack of
      ParticleSoA payload data-member pointers (deposit<&Payload::vx, ...>), plus a dedicated
      depositWeight() for the container-owned weight column;
    - vector fields are per-component (non-contiguous) columns: the loop marshals each
      particle's components to/from a stack Real[] scratch the workers consume, promoting
      (deposit) / demoting (interpolate) between the column type and Real.

  <h3>Destination components</h3>
  Every public deposit/interpolate comes in two forms: a CORE form taking an explicit starting
  mesh component a_comp (the selected columns map to consecutive components
  [a_comp, a_comp + nColumns)), and a CONVENIENCE overload that omits a_comp and uses the first
  components (a_comp == 0).

  <h3>Why depositWeight is separate</h3>
  weight is a container-OWNED mandatory column, not a payload member, so it cannot be named by
  the deposit<&Payload::field> selector. depositWeight() is its entry point. Mechanically it is
  the same as deposit() -- both are thin wrappers that hand depositCore a strength gatherer --
  it just gathers from weightColumn() instead of payload columns.

  See Dev/PORTING_EBParticleMesh.md. Named with the "SoA" suffix so it coexists with the
  production EBParticleMesh (linked via libSource) until design-freeze.

  @note Per project rule, NO default arguments -- callers specify every argument (the "first
        components" behaviour is provided by an overload, not a default).
  @note depositParticleTSC fixes a partition-of-unity bug present in production (see the kernel).
*/

#ifndef CD_EBPARTICLEMESHSOA_H
#define CD_EBPARTICLEMESHSOA_H

// Std includes
#include <cstddef>
#include <type_traits>
#include <utility>

// Chombo includes
#include <CH_Timer.H>
#include <ProblemDomain.H>
#include <EBISBox.H>
#include <EBCellFAB.H>
#include <RealVect.H>
#include <Box.H>
#include <IntVect.H>

// Our includes
#include <CD_DepositionType.H>
#include <CD_ParticleOps.H>
#include <CD_BoxLoops.H>
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Dev port of EBParticleMesh that deposits/interpolates ParticleSoA leaves on one patch.
*/
class EBParticleMeshSoA
{
public:
  /**
    @brief Default constructor. Must subsequently call define.
  */
  EBParticleMeshSoA() : m_verbose(false)
  {}

  /**
    @brief Full constructor; calls define.
    @param[in] a_domain  Computational domain.
    @param[in] a_region  Cell-centered valid box.
    @param[in] a_ebisbox EBIS box.
    @param[in] a_dx      Grid resolution.
    @param[in] a_probLo  Lower-left corner of the computational domain.
  */
  EBParticleMeshSoA(const ProblemDomain& a_domain,
                    const Box&           a_region,
                    const EBISBox&       a_ebisbox,
                    const RealVect&      a_dx,
                    const RealVect&      a_probLo)
  {
    this->define(a_domain, a_region, a_ebisbox, a_dx, a_probLo);
  }

  /**
    @brief Define function.
    @param[in] a_domain  Computational domain.
    @param[in] a_region  Cell-centered valid box.
    @param[in] a_ebisbox EBIS box.
    @param[in] a_dx      Grid resolution.
    @param[in] a_probLo  Lower-left corner of the computational domain.
  */
  void
  define(const ProblemDomain& a_domain,
         const Box&           a_region,
         const EBISBox&       a_ebisbox,
         const RealVect&      a_dx,
         const RealVect&      a_probLo)
  {
    CH_TIME("EBParticleMeshSoA::define");

    m_verbose = false;
    m_domain  = a_domain;
    m_region  = a_region;
    m_ebisbox = a_ebisbox;
    m_probLo  = a_probLo;
    m_dx      = a_dx;
  }

  /**
    @brief Deposit the container-owned weight column onto mesh component a_comp.
    @details weight is a mandatory column, not a payload member, so it has its own entry point
    (see the class documentation). Equivalent to deposit() but gathers from weightColumn().
    @tparam P      Payload type of the ParticleSoA.
    @tparam Traits Column descriptor of P.
    @param[in,out] a_meshData       Mesh data (INCREMENTED). Must have at least a_comp+1 components.
    @param[in]     a_comp           Destination mesh component.
    @param[in]     a_particles      Particles to deposit.
    @param[in]     a_depositionType Deposition kernel (NGP/CIC/TSC).
    @param[in]     a_widthScale     Particle-width multiplier.
    @param[in]     a_forceIrregNGP  Force NGP in cut cells.
  */
  template <typename P, typename Traits>
  void
  depositWeight(EBCellFAB&                    a_meshData,
                const int                     a_comp,
                const ParticleSoA<P, Traits>& a_particles,
                const DepositionType          a_depositionType,
                const Real                    a_widthScale,
                const bool                    a_forceIrregNGP) const
  {
    CH_TIME("EBParticleMeshSoA::depositWeight");

    this->depositCore(a_meshData,
                      a_comp,
                      a_particles,
                      1,
                      a_depositionType,
                      a_widthScale,
                      a_forceIrregNGP,
                      [&](const std::size_t a_i, Real* a_out) {
                        a_out[0] = a_particles.weight(a_i);
                      });
  }

  /**
    @brief Deposit the weight column onto the first mesh component (a_comp == 0).
    @tparam P      Payload type of the ParticleSoA.
    @tparam Traits Column descriptor of P.
    @param[in,out] a_meshData       Mesh data (INCREMENTED).
    @param[in]     a_particles      Particles to deposit.
    @param[in]     a_depositionType Deposition kernel (NGP/CIC/TSC).
    @param[in]     a_widthScale     Particle-width multiplier.
    @param[in]     a_forceIrregNGP  Force NGP in cut cells.
  */
  template <typename P, typename Traits>
  void
  depositWeight(EBCellFAB&                    a_meshData,
                const ParticleSoA<P, Traits>& a_particles,
                const DepositionType          a_depositionType,
                const Real                    a_widthScale,
                const bool                    a_forceIrregNGP) const
  {
    CH_TIME("EBParticleMeshSoA::depositWeight");

    this->depositWeight(a_meshData, 0, a_particles, a_depositionType, a_widthScale, a_forceIrregNGP);
  }

  /**
    @brief Deposit one or more payload columns onto mesh components [a_comp, a_comp + nColumns).
    @details nColumns = sizeof...(Members) (1 = scalar, SpaceDim = vector). The selected columns
    map to consecutive mesh components starting at a_comp.
    @tparam Members Pack of payload data-member pointers selecting the column(s) to deposit.
    @tparam P       Payload type (deduced).
    @tparam Traits  Column descriptor (deduced).
    @param[in,out] a_meshData       Mesh data (INCREMENTED). Needs at least a_comp+sizeof...(Members) components.
    @param[in]     a_comp           First destination mesh component.
    @param[in]     a_particles      Particles to deposit.
    @param[in]     a_depositionType Deposition kernel (NGP/CIC/TSC).
    @param[in]     a_widthScale     Particle-width multiplier.
    @param[in]     a_forceIrregNGP  Force NGP in cut cells.
  */
  template <auto... Members, typename P, typename Traits>
  void
  deposit(EBCellFAB&                    a_meshData,
          const int                     a_comp,
          const ParticleSoA<P, Traits>& a_particles,
          const DepositionType          a_depositionType,
          const Real                    a_widthScale,
          const bool                    a_forceIrregNGP) const
  {
    CH_TIME("EBParticleMeshSoA::deposit");

    static_assert(sizeof...(Members) == 1 || sizeof...(Members) == SpaceDim,
                  "deposit selects either one scalar column or SpaceDim vector-component columns");

    this->depositCore(a_meshData,
                      a_comp,
                      a_particles,
                      static_cast<int>(sizeof...(Members)),
                      a_depositionType,
                      a_widthScale,
                      a_forceIrregNGP,
                      [&](const std::size_t a_i, Real* a_out) {
                        std::size_t k = 0;
                        ((a_out[k++] = static_cast<Real>(a_particles.template get<Members>(a_i))), ...);
                      });
  }

  /**
    @brief Deposit payload columns onto the first mesh components (a_comp == 0).
    @tparam Members Pack of payload data-member pointers selecting the column(s) to deposit.
    @tparam P       Payload type (deduced).
    @tparam Traits  Column descriptor (deduced).
    @param[in,out] a_meshData       Mesh data (INCREMENTED).
    @param[in]     a_particles      Particles to deposit.
    @param[in]     a_depositionType Deposition kernel (NGP/CIC/TSC).
    @param[in]     a_widthScale     Particle-width multiplier.
    @param[in]     a_forceIrregNGP  Force NGP in cut cells.
  */
  template <auto... Members, typename P, typename Traits>
  void
  deposit(EBCellFAB&                    a_meshData,
          const ParticleSoA<P, Traits>& a_particles,
          const DepositionType          a_depositionType,
          const Real                    a_widthScale,
          const bool                    a_forceIrregNGP) const
  {
    CH_TIME("EBParticleMeshSoA::deposit");

    this->deposit<Members...>(a_meshData, 0, a_particles, a_depositionType, a_widthScale, a_forceIrregNGP);
  }

  /**
    @brief Interpolate mesh components [a_comp, a_comp + nColumns) onto one or more payload columns.
    @details nColumns = sizeof...(Members) (1 = scalar, SpaceDim = vector). Consecutive mesh
    components starting at a_comp are read into the selected columns.
    @tparam Members Pack of payload data-member pointers selecting the target column(s).
    @tparam P       Payload type (deduced).
    @tparam Traits  Column descriptor (deduced).
    @param[in,out] a_particles     Particles whose target columns are overwritten.
    @param[in]     a_meshData      Mesh data. Needs at least a_comp+sizeof...(Members) components.
    @param[in]     a_comp          First source mesh component.
    @param[in]     a_interpType    Interpolation kernel (NGP/CIC/TSC).
    @param[in]     a_forceIrregNGP Force NGP in cut cells.
  */
  template <auto... Members, typename P, typename Traits>
  void
  interpolate(ParticleSoA<P, Traits>& a_particles,
              const EBCellFAB&        a_meshData,
              const int               a_comp,
              const DepositionType    a_interpType,
              const bool              a_forceIrregNGP) const
  {
    CH_TIME("EBParticleMeshSoA::interpolate");

    static_assert(sizeof...(Members) == 1 || sizeof...(Members) == SpaceDim,
                  "interpolate targets either one scalar column or SpaceDim vector-component columns");

    this->interpolateCore(
      a_particles,
      a_meshData,
      a_comp,
      static_cast<int>(sizeof...(Members)),
      a_interpType,
      a_forceIrregNGP,
      [&](const std::size_t a_i, const Real* a_in) {
        std::size_t k = 0;
        ((a_particles.template get<Members>(
            a_i) = static_cast<std::decay_t<decltype(a_particles.template get<Members>(a_i))>>(a_in[k++])),
         ...);
      });
  }

  /**
    @brief Interpolate the first mesh components (a_comp == 0) onto payload columns.
    @tparam Members Pack of payload data-member pointers selecting the target column(s).
    @tparam P       Payload type (deduced).
    @tparam Traits  Column descriptor (deduced).
    @param[in,out] a_particles     Particles whose target columns are overwritten.
    @param[in]     a_meshData      Mesh data.
    @param[in]     a_interpType    Interpolation kernel (NGP/CIC/TSC).
    @param[in]     a_forceIrregNGP Force NGP in cut cells.
  */
  template <auto... Members, typename P, typename Traits>
  void
  interpolate(ParticleSoA<P, Traits>& a_particles,
              const EBCellFAB&        a_meshData,
              const DepositionType    a_interpType,
              const bool              a_forceIrregNGP) const
  {
    CH_TIME("EBParticleMeshSoA::interpolate");

    this->interpolate<Members...>(a_particles, a_meshData, 0, a_interpType, a_forceIrregNGP);
  }

  /**
    @brief Interpolate mesh component a_comp onto the container-owned weight column.
    @details weight is a mandatory column, not a payload member, so it has its own entry point (the
    interpolation mirror of depositWeight). Equivalent to interpolate() but scatters into weight().
    @tparam P      Payload type of the ParticleSoA.
    @tparam Traits Column descriptor of P.
    @param[in,out] a_particles     Particles whose weight column is overwritten.
    @param[in]     a_meshData      Mesh data. Must have at least a_comp+1 components.
    @param[in]     a_comp          Source mesh component.
    @param[in]     a_interpType    Interpolation kernel (NGP/CIC/TSC).
    @param[in]     a_forceIrregNGP Force NGP in cut cells.
  */
  template <typename P, typename Traits>
  void
  interpolateWeight(ParticleSoA<P, Traits>& a_particles,
                    const EBCellFAB&        a_meshData,
                    const int               a_comp,
                    const DepositionType    a_interpType,
                    const bool              a_forceIrregNGP) const
  {
    CH_TIME("EBParticleMeshSoA::interpolateWeight");

    this->interpolateCore(a_particles,
                          a_meshData,
                          a_comp,
                          1,
                          a_interpType,
                          a_forceIrregNGP,
                          [&](const std::size_t a_i, const Real* a_in) {
                            a_particles.weight(a_i) = a_in[0];
                          });
  }

  /**
    @brief Interpolate the first mesh component (a_comp == 0) onto the weight column.
    @tparam P      Payload type of the ParticleSoA.
    @tparam Traits Column descriptor of P.
    @param[in,out] a_particles     Particles whose weight column is overwritten.
    @param[in]     a_meshData      Mesh data.
    @param[in]     a_interpType    Interpolation kernel (NGP/CIC/TSC).
    @param[in]     a_forceIrregNGP Force NGP in cut cells.
  */
  template <typename P, typename Traits>
  void
  interpolateWeight(ParticleSoA<P, Traits>& a_particles,
                    const EBCellFAB&        a_meshData,
                    const DepositionType    a_interpType,
                    const bool              a_forceIrregNGP) const
  {
    CH_TIME("EBParticleMeshSoA::interpolateWeight");

    this->interpolateWeight(a_particles, a_meshData, 0, a_interpType, a_forceIrregNGP);
  }

protected:
  /**
    @brief Computational domain.
  */
  ProblemDomain m_domain;

  /**
    @brief Verbose flag.
  */
  bool m_verbose;

  /**
    @brief Cell-centered valid box.
  */
  Box m_region;

  /**
    @brief EBIS box.
  */
  EBISBox m_ebisbox;

  /**
    @brief Grid resolution.
  */
  RealVect m_dx;

  /**
    @brief Lower-left corner of the computational domain.
  */
  RealVect m_probLo;

  /**
    @brief Shared deposition loop: gathers each particle's strength via a_gather, then calls the
           storage-agnostic per-particle kernel for the selected deposition type.
    @tparam P          Payload type.
    @tparam Traits     Column descriptor.
    @tparam StrengthFn Callable (std::size_t i, Real* out) filling out[0..numComp).
    @param[in,out] a_meshData       Mesh data (INCREMENTED).
    @param[in]     a_comp           First destination mesh component.
    @param[in]     a_particles      Particles to deposit.
    @param[in]     a_numComp        Number of components to deposit (1 or SpaceDim).
    @param[in]     a_depositionType Deposition kernel.
    @param[in]     a_widthScale     Particle-width multiplier.
    @param[in]     a_forceIrregNGP  Force NGP in cut cells.
    @param[in]     a_gather         Strength gatherer (promotes column values to Real).
  */
  template <typename P, typename Traits, typename StrengthFn>
  void
  depositCore(EBCellFAB&                    a_meshData,
              const int                     a_comp,
              const ParticleSoA<P, Traits>& a_particles,
              const int                     a_numComp,
              const DepositionType          a_depositionType,
              const Real                    a_widthScale,
              const bool                    a_forceIrregNGP,
              StrengthFn                    a_gather) const
  {
    CH_TIME("EBParticleMeshSoA::depositCore");

    CH_assert(a_meshData.nComp() >= a_comp + a_numComp);

    Real invVol = 1.0;
    for (int dir = 0; dir < SpaceDim; dir++) {
      invVol /= m_dx[dir];
    }

    const RealVect cicWidth = 1 * a_widthScale * RealVect::Unit;
    const RealVect tscWidth = 2 * a_widthScale * RealVect::Unit;

    Real ngpFactor = invVol;
    Real cicFactor = invVol;
    Real tscFactor = invVol;
    for (int dir = 0; dir < SpaceDim; dir++) {
      cicFactor *= 1.0 / cicWidth[dir];
      tscFactor *= 2.0 / tscWidth[dir];
    }

    const std::size_t numParticles = a_particles.size();

    Real strength[SpaceDim];

    switch (a_depositionType) {
    case DepositionType::NGP: {
      for (std::size_t i = 0; i < numParticles; i++) {
        a_gather(i, strength);
        this->depositParticleNGP(a_meshData, a_comp, a_particles.position(i), ngpFactor, strength, a_numComp);
      }
      break;
    }
    case DepositionType::CIC: {
      for (std::size_t i = 0; i < numParticles; i++) {
        a_gather(i, strength);
        this->depositParticleCIC(a_meshData,
                                 a_comp,
                                 a_particles.position(i),
                                 cicWidth,
                                 cicFactor,
                                 strength,
                                 a_numComp,
                                 a_forceIrregNGP);
      }
      break;
    }
    case DepositionType::TSC: {
      for (std::size_t i = 0; i < numParticles; i++) {
        a_gather(i, strength);
        this->depositParticleTSC(a_meshData,
                                 a_comp,
                                 a_particles.position(i),
                                 tscWidth,
                                 tscFactor,
                                 strength,
                                 a_numComp,
                                 a_forceIrregNGP);
      }
      break;
    }
    default: {
      MayDay::Abort("EBParticleMeshSoA::depositCore -- logic bust");
      break;
    }
    }
  }

  /**
    @brief Shared interpolation loop: calls the storage-agnostic per-particle kernel into a stack
           scratch, then scatters the result into the selected column(s) via a_scatter.
    @tparam P         Payload type.
    @tparam Traits    Column descriptor.
    @tparam ScatterFn Callable (std::size_t i, const Real* in) writing in[0..numComp) to columns.
    @param[in,out] a_particles     Particles whose target columns are overwritten.
    @param[in]     a_meshData      Mesh data.
    @param[in]     a_comp          First source mesh component.
    @param[in]     a_numComp       Number of components (1 or SpaceDim).
    @param[in]     a_interpType    Interpolation kernel.
    @param[in]     a_forceIrregNGP Force NGP in cut cells.
    @param[in]     a_scatter       Result scatterer (demotes Real to the column type).
  */
  template <typename P, typename Traits, typename ScatterFn>
  void
  interpolateCore(ParticleSoA<P, Traits>& a_particles,
                  const EBCellFAB&        a_meshData,
                  const int               a_comp,
                  const int               a_numComp,
                  const DepositionType    a_interpType,
                  const bool              a_forceIrregNGP,
                  ScatterFn               a_scatter) const
  {
    CH_TIME("EBParticleMeshSoA::interpolateCore");

    CH_assert(a_meshData.nComp() >= a_comp + a_numComp);

    // validBox is where the chosen scheme has enough cells; gatherBox is the read stencil.
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
      MayDay::Error("EBParticleMeshSoA::interpolateCore - logic bust");
    }
    }

    const std::size_t numParticles = a_particles.size();

    Real field[SpaceDim];

    switch (a_interpType) {
    case DepositionType::NGP: {
      for (std::size_t i = 0; i < numParticles; i++) {
        this->interpolateParticleNGP(field, a_meshData, a_comp, a_particles.position(i), a_numComp);
        a_scatter(i, field);
      }
      break;
    }
    case DepositionType::CIC: {
      for (std::size_t i = 0; i < numParticles; i++) {
        this->interpolateParticleCIC(field,
                                     a_meshData,
                                     a_comp,
                                     validBox,
                                     gatherBox,
                                     a_particles.position(i),
                                     a_numComp,
                                     a_forceIrregNGP);
        a_scatter(i, field);
      }
      break;
    }
    case DepositionType::TSC: {
      for (std::size_t i = 0; i < numParticles; i++) {
        this->interpolateParticleTSC(field,
                                     a_meshData,
                                     a_comp,
                                     validBox,
                                     gatherBox,
                                     a_particles.position(i),
                                     a_numComp,
                                     a_forceIrregNGP);
        a_scatter(i, field);
      }
      break;
    }
    default: {
      MayDay::Abort("EBParticleMeshSoA::interpolateCore -- logic bust");
      break;
    }
    }
  }

  // ===========================================================================
  // The per-particle EB/cut-cell kernels below are adapted from
  // Source/Particle/CD_EBParticleMeshImplem.H -- storage-agnostic, with an added
  // destination-component offset (a_comp) and the TSC partition-of-unity fix.
  // ===========================================================================

  /**
    @brief NGP deposition of one particle.
    @param[in,out] a_rho          Mesh data.
    @param[in]     a_comp         First destination mesh component.
    @param[in]     a_position     Particle position.
    @param[in]     a_volumeFactor Volume factor.
    @param[in]     a_strength     Per-component strength to deposit.
    @param[in]     a_numComp      Number of components.
  */
  inline void
  depositParticleNGP(EBCellFAB&      a_rho,
                     const int       a_comp,
                     const RealVect& a_position,
                     const Real&     a_volumeFactor,
                     const Real*     a_strength,
                     const int&      a_numComp) const noexcept
  {
    const IntVect particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);

    CH_assert(m_region.contains(particleIV));

    FArrayBox& rho = a_rho.getFArrayBox();

    for (int comp = 0; comp < a_numComp; comp++) {
      rho(particleIV, a_comp + comp) += a_strength[comp] * a_volumeFactor;
    }
  }

  /**
    @brief Cell box that exactly covers the support of a cloud of the given per-direction width.
    @details Anchored at the cloud's lower/upper edge cells, so the iterated box contains ONLY the
    cells the cloud can overlap (cells outside it have zero weight and would be wasted scatter
    updates). For the standard widths this is 2 cells/dir (CIC) or 3 (TSC), versus the
    over-provisioned centered 3^D/5^D box; the weight clamp still zeroes any boundary cell, so the
    deposited field is bit-identical.
    @param[in] a_position      Particle position.
    @param[in] a_particleWidth Cloud width per direction (in cells).
    @return Cell box covering the cloud support.
  */
  inline Box
  cloudBox(const RealVect& a_position, const RealVect& a_particleWidth) const noexcept
  {
    RealVect lo;
    RealVect hi;
    for (int dir = 0; dir < SpaceDim; dir++) {
      const Real halfWidth = 0.5 * a_particleWidth[dir] * m_dx[dir];
      lo[dir]              = a_position[dir] - halfWidth;
      hi[dir]              = a_position[dir] + halfWidth;
    }
    return Box(ParticleOps::getParticleCellIndex(lo, m_probLo, m_dx),
               ParticleOps::getParticleCellIndex(hi, m_probLo, m_dx));
  }

  /**
    @brief CIC deposition of one particle.
    @param[in,out] a_rho           Mesh data.
    @param[in]     a_comp          First destination mesh component.
    @param[in]     a_position      Particle position.
    @param[in]     a_particleWidth Particle width (relative to grid resolution).
    @param[in]     a_volumeFactor  Volume factor.
    @param[in]     a_strength      Per-component strength to deposit.
    @param[in]     a_numComp       Number of components.
    @param[in]     a_forceIrregNGP Force NGP in cut cells.
  */
  inline void
  depositParticleCIC(EBCellFAB&      a_rho,
                     const int       a_comp,
                     const RealVect& a_position,
                     const RealVect& a_particleWidth,
                     const Real&     a_volumeFactor,
                     const Real*     a_strength,
                     const int&      a_numComp,
                     const bool      a_forceIrregNGP) const noexcept
  {
    const IntVect particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);

    CH_assert(m_region.contains(particleIV));

    FArrayBox& rho = a_rho.getFArrayBox();

    // forceIrregNGP first so the common (false) case short-circuits before the EBISBox query.
    if (a_forceIrregNGP && m_ebisbox.isIrregular(particleIV)) {
      for (int comp = 0; comp < a_numComp; comp++) {
        rho(particleIV, a_comp + comp) += a_strength[comp] * a_volumeFactor;
      }
    }
    else {
      auto cicKernel = [&](const IntVect& iv) -> void {
        Real weight = a_volumeFactor;

        for (int dir = 0; dir < SpaceDim; dir++) {
          const Real a = (m_probLo[dir] - a_position[dir]) / m_dx[dir] + iv[dir];
          const Real b = a + 1.0;
          const Real L = 0.5 * a_particleWidth[dir];

          weight *= std::max(0.0, std::min(b, L) - std::max(a, -L));
        }

        for (int comp = 0; comp < a_numComp; comp++) {
          rho(iv, a_comp + comp) += weight * a_strength[comp];
        }
      };

      BoxLoops::loop<D_DECL(1, 1, 1)>(this->cloudBox(a_position, a_particleWidth), cicKernel);
    }
  }

  /**
    @brief TSC deposition of one particle.
    @param[in,out] a_rho           Mesh data.
    @param[in]     a_comp          First destination mesh component.
    @param[in]     a_position      Particle position.
    @param[in]     a_particleWidth Particle width (relative to grid resolution).
    @param[in]     a_volumeFactor  Volume factor.
    @param[in]     a_strength      Per-component strength to deposit.
    @param[in]     a_numComp       Number of components.
    @param[in]     a_forceIrregNGP Force NGP in cut cells.
  */
  inline void
  depositParticleTSC(EBCellFAB&      a_rho,
                     const int       a_comp,
                     const RealVect& a_position,
                     const RealVect& a_particleWidth,
                     const Real&     a_volumeFactor,
                     const Real*     a_strength,
                     const int&      a_numComp,
                     const bool      a_forceIrregNGP) const noexcept
  {
    const IntVect particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);

    CH_assert(m_region.contains(particleIV));

    FArrayBox& rho = a_rho.getFArrayBox();

    // forceIrregNGP first so the common (false) case short-circuits before the EBISBox query.
    if (a_forceIrregNGP && m_ebisbox.isIrregular(particleIV)) {
      for (int comp = 0; comp < a_numComp; comp++) {
        rho(particleIV, a_comp + comp) += a_strength[comp] * a_volumeFactor;
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

          // factor gates the ENTIRE overlap integral so out-of-support cells contribute zero.
          // This FIXES a partition-of-unity bug in production depositParticleTSC
          // (Source/Particle/CD_EBParticleMeshImplem.H), where factor multiplies only the first
          // term and the cloud over-deposits ~2.33x per dimension. Production is intentionally
          // left unpatched for now (TSC deposit is unused there); see PORTING_EBParticleMesh.md.
          weight *= factor * ((beta - alpha) - (beta * std::abs(beta) - alpha * std::abs(alpha)) / L);
        }

        for (int comp = 0; comp < a_numComp; comp++) {
          rho(iv, a_comp + comp) += weight * a_strength[comp];
        }
      };

      BoxLoops::loop<D_DECL(1, 1, 1)>(this->cloudBox(a_position, a_particleWidth), tscKernel);
    }
  }

  /**
    @brief NGP interpolation onto one particle.
    @param[out] a_particleField Per-component field to fill.
    @param[in]  a_meshData      Mesh data.
    @param[in]  a_comp          First source mesh component.
    @param[in]  a_position      Particle position.
    @param[in]  a_numComp       Number of components.
  */
  inline void
  interpolateParticleNGP(Real*            a_particleField,
                         const EBCellFAB& a_meshData,
                         const int        a_comp,
                         const RealVect&  a_position,
                         const int&       a_numComp) const noexcept
  {
    const IntVect    particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);
    const FArrayBox& meshData   = a_meshData.getFArrayBox();
    const Real       factor     = m_ebisbox.isCovered(particleIV) ? 0.0 : 1.0;

    for (int comp = 0; comp < a_numComp; comp++) {
      a_particleField[comp] = factor * meshData(particleIV, a_comp + comp);
    }
  }

  /**
    @brief CIC interpolation onto one particle.
    @param[out] a_particleField Per-component field to fill.
    @param[in]  a_meshData      Mesh data.
    @param[in]  a_comp          First source mesh component.
    @param[in]  a_validBox      Box in which CIC is valid.
    @param[in]  a_gatherBox     Read stencil box.
    @param[in]  a_position      Particle position.
    @param[in]  a_numComp       Number of components.
    @param[in]  a_forceIrregNGP Force NGP in cut cells.
  */
  inline void
  interpolateParticleCIC(Real*            a_particleField,
                         const EBCellFAB& a_meshData,
                         const int        a_comp,
                         const Box&       a_validBox,
                         const Box&       a_gatherBox,
                         const RealVect&  a_position,
                         const int&       a_numComp,
                         const bool       a_forceIrregNGP) const noexcept
  {
    const IntVect particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);
    const IntVect loIndex    = ParticleOps::getParticleCellIndex(a_position - 0.5 * m_dx, m_probLo, m_dx);

    CH_assert(m_region.contains(particleIV));

    const FArrayBox& meshData = a_meshData.getFArrayBox();

    for (int comp = 0; comp < a_numComp; comp++) {
      a_particleField[comp] = 0.0;
    }

    if ((m_ebisbox.isIrregular(particleIV) && a_forceIrregNGP) || !(a_validBox.contains(particleIV))) {
      for (int comp = 0; comp < a_numComp; comp++) {
        a_particleField[comp] = meshData(particleIV, a_comp + comp);
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
          a_particleField[comp] += weight * meshData(iv, a_comp + comp);
        }
      };

      const Box gatherBox = a_gatherBox + loIndex;

      BoxLoops::loop<D_DECL(1, 1, 1)>(gatherBox, cicKernel);
    }
  }

  /**
    @brief TSC interpolation onto one particle.
    @param[out] a_particleField Per-component field to fill.
    @param[in]  a_meshData      Mesh data.
    @param[in]  a_comp          First source mesh component.
    @param[in]  a_validBox      Box in which TSC is valid.
    @param[in]  a_gatherBox     Read stencil box.
    @param[in]  a_position      Particle position.
    @param[in]  a_numComp       Number of components.
    @param[in]  a_forceIrregNGP Force NGP in cut cells.
  */
  inline void
  interpolateParticleTSC(Real*            a_particleField,
                         const EBCellFAB& a_meshData,
                         const int        a_comp,
                         const Box&       a_validBox,
                         const Box&       a_gatherBox,
                         const RealVect&  a_position,
                         const int&       a_numComp,
                         const bool       a_forceIrregNGP) const noexcept
  {
    const IntVect particleIV = ParticleOps::getParticleCellIndex(a_position, m_probLo, m_dx);
    const IntVect loIndex    = ParticleOps::getParticleCellIndex(a_position - m_dx, m_probLo, m_dx);

    CH_assert(m_region.contains(particleIV));

    const FArrayBox& meshData = a_meshData.getFArrayBox();

    for (int comp = 0; comp < a_numComp; comp++) {
      a_particleField[comp] = 0.0;
    }

    if ((m_ebisbox.isIrregular(particleIV) && a_forceIrregNGP) || !(a_validBox.contains(particleIV))) {
      for (int comp = 0; comp < a_numComp; comp++) {
        a_particleField[comp] = meshData(particleIV, a_comp + comp);
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
          a_particleField[comp] += weight * meshData(iv, a_comp + comp);
        }
      };

      const Box gatherBox = a_gatherBox + loIndex;

      BoxLoops::loop<D_DECL(1, 1, 1)>(gatherBox, tscKernel);
    }
  }
};

#include <CD_NamespaceFooter.H>

#endif
