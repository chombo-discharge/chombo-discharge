/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBAMRParticleMesh.H
 * @brief  Declaration of an AMR driver that deposits/interpolates ParticleContainer particles.
 * @author Robert Marskar
 */

#ifndef CD_EBAMRPARTICLEMESH_H
#define CD_EBAMRPARTICLEMESH_H

// Std includes
#include <algorithm>
#include <cstddef>
#include <map>
#include <vector>

// Chombo includes
#include <CH_Timer.H>
#include <RefCountedPtr.H>
#include <Vector.H>
#include <Copier.H>
#include <DataIterator.H>
#include <NeighborIterator.H>
#include <LevelData.H>
#include <BaseFab.H>
#include <FArrayBox.H>
#include <EBCellFAB.H>
#include <EBLevelGrid.H>
#include <EBAlias.H>
#include <IntVectSet.H>
#include <BoxIterator.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_EBAddOp.H>
#include <CD_DataOps.H>
#include <CD_OpenMP.H>
#include <CD_DepositionType.H>
#include <CD_CoarseFineDeposition.H>
#include <CD_IrregularDeposition.H>
#include <CD_EBCoarseFineParticleMesh.H>
#include <CD_EBParticleMesh.H>
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief AMR driver that deposits/interpolates ParticleContainer particles across the hierarchy.
 * @details
 * EBAMRParticleMesh deposits the particles of a ParticleContainer onto AMR mesh data and interpolates
 * mesh data back onto particle columns, by driving the per-patch EBParticleMesh leaf on every box of
 * every level. define() builds the per-level leaves, the level valid+ghost -> valid Copiers, the
 * EBCoarseFineParticleMesh operators, and the outer-halo masks; it must be called before any deposit or
 * interpolate.
 *
 * <h3>Interpolation</h3>
 * interpolate() and interpolateWeight() are purely per-level/per-patch: the mesh is assumed to be
 * valid+ghosted by the caller, so no coarse-fine coupling is needed.
 *
 * <h3>Deposition</h3>
 * deposit(), depositWeight(), and depositGathered() conserve mass across refinement boundaries, and
 * select how the coarse-fine overlap is handled through @c CoarseFineDeposition:
 * - <b>Interp</b>     -- move the over-boundary mass between levels with the coarse-fine mesh operators.
 * - <b>Halo</b>       -- copy the coarse-side halo particles into the finer patches their clouds reach, as
 *                        coarse ghost particles, and deposit them there with the coarse width.
 * - <b>HaloNGP</b>    -- deposit the near-boundary halo particles with NGP so their clouds do not spread
 *                        across the boundary.
 * - <b>Transition</b> -- deposit the coarse-side transition particles with the fine-grid width, on their own
 *                        level and as coarse ghost particles on the finer one.
 * Halo, HaloNGP, and Transition require at least one ghost cell. Halo and Transition need the realm's
 * coarse-to-fine particle ghost masks (see getCoarseToFineMaskWidths).
 *
 * <h3>Cut cells, and what this class does NOT do</h3>
 * The deposition routines take an @c IrregularDeposition, but this class only implements the part of it that the
 * leaf kernels can answer: NGP in a cut cell, yes or no. Everything else is the CALLER's responsibility.
 * - <b>Native</b>, <b>NGP</b> -- handled here. These are the two kernel behaviours.
 * - <b>Mirror</b> -- NOT handled here. It is a second pass that reflects each band particle across the embedded
 *                    boundary and deposits the image, so it needs an image container, a scratch field and a
 *                    collective remap, none of which this class owns. See @c ItoSolver::mirrorPass.
 *                    The one part handled here: Interp's coarse-to-fine step does not rescale cut cells by volume
 *                    fraction, because with the mirror they hold the density rather than kappa times it.
 * - <b>Redistribute</b>, <b>RedistributeBlended</b> -- NOT handled here. They act on the field AFTER the deposit.
 *                    See @c ItoSolver::redistributeAMR and @c McPhoto::depositHybridDivergence.
 *
 * The consequence is worth stating plainly, because it is silent: passing @c Mirror or @c Redistribute to a
 * deposition routine here and then omitting the caller's half produces @c Native behaviour (for @c Mirror with
 * Interp, Native without the volume-fraction rescaling), with no error and no diagnostic. A caller that selects
 * one of those values owns the rest of its meaning.
 *
 * This split is deliberate but not ideal -- the enum is named in one place and interpreted in three, so there is no
 * single source of truth for what a selector means. Centralizing the mirror and the redistribution here is the
 * natural end state; it needs @c ParticleContainer::define's minimum block size, level tiles, valid cells and realm
 * name, which this class is not currently given.
 */
class EBAMRParticleMesh
{
public:
  /**
   * @name Construction and lifetime
   */
  ///@{

  /**
   * @brief Default constructor. Leaves the object undefined; call define().
   */
  EBAMRParticleMesh() : m_ghost(0), m_finestLevel(-1), m_isDefined(false)
  {}

  /**
   * @brief Destructor.
   */
  ~EBAMRParticleMesh() = default;

  /**
   * @brief Copy construction is deleted (shared RefCountedPtr state).
   */
  EBAMRParticleMesh(const EBAMRParticleMesh&) = delete;

  /**
   * @brief Copy assignment is deleted.
   * @return Reference to this.
   */
  EBAMRParticleMesh&
  operator=(const EBAMRParticleMesh&) = delete;

  /**
   * @brief Define the per-level leaves, coarse-fine mesh operators, and halo/transition masks.
   * @param[in] a_eblgs       Per-level EB level grids.
   * @param[in] a_refRat      Per-level refinement ratio (entry l is the ratio between level l and l+1).
   * @param[in] a_dx          Per-level (isotropic) grid spacing.
   * @param[in] a_probLo      Lower-left corner of the physical domain.
   * @param[in] a_ghost       Number of ghost cells in the mesh data holders.
   * @param[in] a_finestLevel Finest AMR level index.
   */
  void
  define(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
         const Vector<int>&                        a_refRat,
         const Vector<Real>&                       a_dx,
         const RealVect&                           a_probLo,
         const int                                 a_ghost,
         const int                                 a_finestLevel)
  {
    CH_TIME("EBAMRParticleMesh::define");

    CH_assert(a_eblgs.size() > a_finestLevel);
    CH_assert(a_refRat.size() > a_finestLevel);
    CH_assert(a_dx.size() > a_finestLevel);

    m_eblgs       = a_eblgs;
    m_refRat      = a_refRat;
    m_dx          = a_dx;
    m_probLo      = a_probLo;
    m_ghost       = a_ghost;
    m_finestLevel = a_finestLevel;

    this->defineCoarseFineMotion();
    this->defineEBParticleMesh();
    this->defineLevelMotion();
    this->defineOuterHaloMasks();

    m_isDefined = true;
  }

  ///@}

  /**
   * @name Interpolation
   */
  ///@{

  /**
   * @brief Interpolate a mesh field onto one or more payload columns, across all levels.
   * @details Per level, per patch: calls EBParticleMesh::interpolate. The mesh is assumed to be
   * valid+ghosted by the caller (no coarse-fine coupling is needed for interpolation).
   * @tparam Members Pack of payload data-member pointers selecting the target column(s).
   * @tparam P       Payload type (deduced).
   * @tparam Traits  Column descriptor (deduced).
   * @param[in,out] a_particles     Particles whose target columns are overwritten.
   * @param[in]     a_meshData      AMR mesh data (>= sizeof...(Members) components).
   * @param[in]     a_interpType    Interpolation kernel (NGP/CIC/TSC).
   * @param[in]     a_forceIrregNGP Force NGP in cut cells.
   */
  template <auto... Members, typename P, typename Traits>
  void
  interpolate(ParticleContainer<P, Traits>& a_particles,
              const EBAMRCellData&          a_meshData,
              const DepositionType          a_interpType,
              const bool                    a_forceIrregNGP) const
  {
    CH_TIME("EBAMRParticleMesh::interpolate");

    CH_assert(m_isDefined);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl = m_eblgs[lvl]->getDBL();
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&      din    = dit[mybox];
        const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];
        const EBCellFAB&      data   = (*a_meshData[lvl])[din];

        interp.template interpolate<Members...>(a_particles[lvl][din], data, a_interpType, a_forceIrregNGP);
      }
    }
  }

  /**
   * @brief Interpolate a mesh field onto the container-owned weight column, across all levels.
   * @details Per level, per patch: calls EBParticleMesh::interpolateWeight. The mesh is assumed to
   * be valid+ghosted by the caller (no coarse-fine coupling is needed for interpolation). This is the
   * interpolation mirror of depositWeight.
   * @tparam P      Payload type (deduced).
   * @tparam Traits Column descriptor (deduced).
   * @param[in,out] a_particles     Particles whose weight column is overwritten.
   * @param[in]     a_meshData      AMR mesh data (>= 1 component).
   * @param[in]     a_interpType    Interpolation kernel (NGP/CIC/TSC).
   * @param[in]     a_forceIrregNGP Force NGP in cut cells.
   */
  template <typename P, typename Traits>
  void
  interpolateWeight(ParticleContainer<P, Traits>& a_particles,
                    const EBAMRCellData&          a_meshData,
                    const DepositionType          a_interpType,
                    const bool                    a_forceIrregNGP) const
  {
    CH_TIME("EBAMRParticleMesh::interpolateWeight");

    CH_assert(m_isDefined);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl = m_eblgs[lvl]->getDBL();
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&      din    = dit[mybox];
        const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];
        const EBCellFAB&      data   = (*a_meshData[lvl])[din];

        interp.interpolateWeight(a_particles[lvl][din], data, a_interpType, a_forceIrregNGP);
      }
    }
  }

  ///@}

  /**
   * @name Deposition (all coarse-fine strategies)
   */
  ///@{

  /**
   * @brief Deposit the container-owned weight onto the AMR mesh.
   * @details Dispatches on the coarse-fine strategy (Interp / Halo / HaloNGP / Transition); see the
   * class documentation. Halo/HaloNGP/Transition need ghost >= 1.
   * @tparam P          Payload type (deduced).
   * @tparam Traits     Column descriptor (deduced).
   * @tparam MaskGetter Callable (int width) -> const AMRParticleGhostMask&.
   * @param[in,out] a_meshData             AMR mesh data (one component; RESET then filled).
   * @param[in]     a_particles            Particles to deposit.
   * @param[in]     a_depositionType       Deposition kernel (NGP/CIC/TSC).
   * @param[in]     a_coarseFineDeposition Coarse-fine handling (Interp / Halo / HaloNGP / Transition).
   * @param[in]     a_irregularDeposition        Cut-cell deposition strategy.
   * @param[in]     a_coarToFineMasks      Callable int -> const AMRParticleGhostMask&, giving the realm's
   * coarse-to-fine particle ghost mask of a width (see getCoarseToFineMaskWidths).
   */
  template <typename P, typename Traits, typename MaskGetter>
  void
  depositWeight(EBAMRCellData&                      a_meshData,
                const ParticleContainer<P, Traits>& a_particles,
                const DepositionType                a_depositionType,
                const CoarseFineDeposition          a_coarseFineDeposition,
                const IrregularDeposition           a_irregularDeposition,
                const MaskGetter&                   a_coarToFineMasks)
  {
    CH_TIME("EBAMRParticleMesh::depositWeight");

    CH_assert(m_isDefined);

    // The leaf kernels answer one question -- NGP in a cut cell, yes or no -- so the selector collapses to a bool
    // here. Mirror and Redistribute leave the kernel alone and are the CALLER's responsibility; see the class
    // documentation for why that is silent when the caller forgets.
    const bool forceIrregNGP = (a_irregularDeposition == IrregularDeposition::NGP);

    // Mirror deposition makes the cut cells hold the density, which Interp's coarse-to-fine step must not rescale.
    const bool cutCellsHoldDensity = (a_irregularDeposition == IrregularDeposition::Mirror);

    const auto patch =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.depositWeight(a_rho, a_soa, a_depositionType, a_widthScale, forceIrregNGP);
      };
    const auto patchNGP =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.depositWeight(a_rho, a_soa, DepositionType::NGP, a_widthScale, forceIrregNGP);
      };
    const auto patchGhosts =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.depositCoarseGhostsWeight(a_rho, 0, a_soa, a_depositionType, a_widthScale);
      };

    switch (a_coarseFineDeposition) {
    case CoarseFineDeposition::Interp: {
      this->depositInterpCore(a_meshData, a_particles, 1, cutCellsHoldDensity, patch);

      break;
    }
    case CoarseFineDeposition::Halo: {
      this->template depositHaloCore<1>(a_meshData,
                                        a_particles,
                                        a_depositionType,
                                        forceIrregNGP,
                                        a_coarToFineMasks,
                                        patch,
                                        patchGhosts);

      break;
    }
    case CoarseFineDeposition::HaloNGP: {
      this->template depositHaloNGPCore<1>(a_meshData, a_particles, patch, patchNGP);

      break;
    }
    case CoarseFineDeposition::Transition: {
      this->template depositTransitionCore<1>(a_meshData,
                                              a_particles,
                                              a_depositionType,
                                              forceIrregNGP,
                                              cutCellsHoldDensity,
                                              a_coarToFineMasks,
                                              patch,
                                              patchGhosts);

      break;
    }
    default: {
      this->errorUnsupported();

      break;
    }
    }
  }

  /**
   * @brief Deposit a custom per-particle scalar (computed by a_cellGather) onto the AMR mesh.
   * @details Generalizes depositWeight() to an arbitrary per-particle value (a product of columns such as
   * weight*mobility), dispatching on the coarse-fine strategy exactly like depositWeight(). Deposits one
   * mesh component.
   * @tparam P          Payload type (deduced).
   * @tparam Traits     Column descriptor (deduced).
   * @tparam CellGather Callable (const ParticleSoA<P, Traits>& a_leaf, std::size_t a_i) -> Real.
   * @tparam MaskGetter Callable (int width) -> const AMRParticleGhostMask&.
   * @param[in,out] a_meshData             AMR mesh data (one component; RESET then filled).
   * @param[in]     a_particles            Particles to deposit.
   * @param[in]     a_depositionType       Deposition kernel (NGP/CIC/TSC).
   * @param[in]     a_coarseFineDeposition Coarse-fine handling (Interp / Halo / HaloNGP / Transition).
   * @param[in]     a_irregularDeposition        Cut-cell deposition strategy.
   * @param[in]     a_cellGather           Per-particle value gatherer (leaf, index) -> Real.
   * @param[in]     a_coarToFineMasks      Callable int -> const AMRParticleGhostMask&, giving the realm's
   * coarse-to-fine particle ghost mask of a width (see getCoarseToFineMaskWidths).
   */
  template <typename P, typename Traits, typename CellGather, typename MaskGetter>
  void
  depositGathered(EBAMRCellData&                      a_meshData,
                  const ParticleContainer<P, Traits>& a_particles,
                  const DepositionType                a_depositionType,
                  const CoarseFineDeposition          a_coarseFineDeposition,
                  const IrregularDeposition           a_irregularDeposition,
                  CellGather                          a_cellGather,
                  const MaskGetter&                   a_coarToFineMasks)
  {
    CH_TIME("EBAMRParticleMesh::depositGathered");

    CH_assert(m_isDefined);

    // See depositWeight for why the selector collapses to a bool here.
    const bool forceIrregNGP = (a_irregularDeposition == IrregularDeposition::NGP);

    // Mirror deposition makes the cut cells hold the density, which Interp's coarse-to-fine step must not rescale.
    const bool cutCellsHoldDensity = (a_irregularDeposition == IrregularDeposition::Mirror);

    const auto patch =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.depositGathered(a_rho,
                                 0,
                                 a_soa,
                                 a_depositionType,
                                 a_widthScale,
                                 forceIrregNGP,
                                 [&](const std::size_t a_i, Real* a_out) {
                                   a_out[0] = a_cellGather(a_soa, a_i);
                                 });
      };
    const auto patchNGP =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.depositGathered(a_rho,
                                 0,
                                 a_soa,
                                 DepositionType::NGP,
                                 a_widthScale,
                                 forceIrregNGP,
                                 [&](const std::size_t a_i, Real* a_out) {
                                   a_out[0] = a_cellGather(a_soa, a_i);
                                 });
      };
    const auto patchGhosts =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.depositCoarseGhostsGathered(a_rho,
                                             0,
                                             a_soa,
                                             a_depositionType,
                                             a_widthScale,
                                             [&](const std::size_t a_i, Real* a_out) {
                                               a_out[0] = a_cellGather(a_soa, a_i);
                                             });
      };

    switch (a_coarseFineDeposition) {
    case CoarseFineDeposition::Interp: {
      this->depositInterpCore(a_meshData, a_particles, 1, cutCellsHoldDensity, patch);

      break;
    }
    case CoarseFineDeposition::Halo: {
      this->template depositHaloCore<1>(a_meshData,
                                        a_particles,
                                        a_depositionType,
                                        forceIrregNGP,
                                        a_coarToFineMasks,
                                        patch,
                                        patchGhosts);

      break;
    }
    case CoarseFineDeposition::HaloNGP: {
      this->template depositHaloNGPCore<1>(a_meshData, a_particles, patch, patchNGP);

      break;
    }
    case CoarseFineDeposition::Transition: {
      this->template depositTransitionCore<1>(a_meshData,
                                              a_particles,
                                              a_depositionType,
                                              forceIrregNGP,
                                              cutCellsHoldDensity,
                                              a_coarToFineMasks,
                                              patch,
                                              patchGhosts);

      break;
    }
    default: {
      this->errorUnsupported();

      break;
    }
    }
  }

  /**
   * @brief Deposit one or more payload columns onto the AMR mesh.
   * @details See depositWeight for the coarse-fine dispatch.
   * @tparam Members    Pack of payload data-member pointers selecting the column(s) to deposit.
   * @tparam P          Payload type (deduced).
   * @tparam Traits     Column descriptor (deduced).
   * @tparam MaskGetter Callable (int width) -> const AMRParticleGhostMask&.
   * @param[in,out] a_meshData             AMR mesh data (sizeof...(Members) components; RESET then filled).
   * @param[in]     a_particles            Particles to deposit.
   * @param[in]     a_depositionType       Deposition kernel (NGP/CIC/TSC).
   * @param[in]     a_coarseFineDeposition Coarse-fine handling (Interp / Halo / HaloNGP / Transition).
   * @param[in]     a_irregularDeposition        Cut-cell deposition strategy.
   * @param[in]     a_coarToFineMasks      Callable int -> const AMRParticleGhostMask&, giving the realm's
   * coarse-to-fine particle ghost mask of a width (see getCoarseToFineMaskWidths).
   */
  template <auto... Members, typename P, typename Traits, typename MaskGetter>
  void
  deposit(EBAMRCellData&                      a_meshData,
          const ParticleContainer<P, Traits>& a_particles,
          const DepositionType                a_depositionType,
          const CoarseFineDeposition          a_coarseFineDeposition,
          const IrregularDeposition           a_irregularDeposition,
          const MaskGetter&                   a_coarToFineMasks)
  {
    CH_TIME("EBAMRParticleMesh::deposit");

    CH_assert(m_isDefined);

    // See depositWeight for why the selector collapses to a bool here.
    const bool forceIrregNGP = (a_irregularDeposition == IrregularDeposition::NGP);

    // Mirror deposition makes the cut cells hold the density, which Interp's coarse-to-fine step must not rescale.
    const bool cutCellsHoldDensity = (a_irregularDeposition == IrregularDeposition::Mirror);

    const auto patch =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.template deposit<Members...>(a_rho, a_soa, a_depositionType, a_widthScale, forceIrregNGP);
      };
    const auto patchNGP =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.template deposit<Members...>(a_rho, a_soa, DepositionType::NGP, a_widthScale, forceIrregNGP);
      };
    const auto patchGhosts =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.template depositCoarseGhosts<Members...>(a_rho, 0, a_soa, a_depositionType, a_widthScale);
      };

    switch (a_coarseFineDeposition) {
    case CoarseFineDeposition::Interp: {
      this->depositInterpCore(a_meshData,
                              a_particles,
                              static_cast<int>(sizeof...(Members)),
                              cutCellsHoldDensity,
                              patch);

      break;
    }
    case CoarseFineDeposition::Halo: {
      this->template depositHaloCore<static_cast<int>(sizeof...(Members))>(a_meshData,
                                                                           a_particles,
                                                                           a_depositionType,
                                                                           forceIrregNGP,
                                                                           a_coarToFineMasks,
                                                                           patch,
                                                                           patchGhosts);

      break;
    }
    case CoarseFineDeposition::HaloNGP: {
      this->template depositHaloNGPCore<static_cast<int>(sizeof...(Members))>(a_meshData, a_particles, patch, patchNGP);

      break;
    }
    case CoarseFineDeposition::Transition: {
      this->template depositTransitionCore<static_cast<int>(sizeof...(Members))>(a_meshData,
                                                                                 a_particles,
                                                                                 a_depositionType,
                                                                                 forceIrregNGP,
                                                                                 cutCellsHoldDensity,
                                                                                 a_coarToFineMasks,
                                                                                 patch,
                                                                                 patchGhosts);
      break;
    }
    default: {
      this->errorUnsupported();

      break;
    }
    }
  }

  /**
   * @brief Width, in finer cells, of the coarse-to-fine particle ghost mask that Halo deposition fills with.
   * @details The reach into the finer level of a coarser-width cloud: refRat/2 for CIC, refRat for TSC, and 0 for NGP,
   * whose cloud never leaves the particle's own cell.
   * @param[in] a_depositionType Deposition kernel.
   * @param[in] a_refRat         Refinement ratio between the level and the finer one.
   * @return The mask width.
   */
  static int
  getHaloGhostWidth(const DepositionType a_depositionType, const int a_refRat) noexcept
  {
    switch (a_depositionType) {
    case DepositionType::CIC: {
      return a_refRat / 2;
    }
    case DepositionType::TSC: {
      return a_refRat;
    }
    default: {
      return 0;
    }
    }
  }

  /**
   * @brief Transition band width, in finer cells: the coarser particles within this distance of the finer grid deposit
   * with the finer-level width.
   * @details refRat/2 for CIC and refRat for TSC. NGP has no band (Transition with NGP deposits as Interp).
   * @param[in] a_depositionType Deposition kernel.
   * @param[in] a_refRat         Refinement ratio between the level and the finer one.
   * @return The band width.
   */
  static int
  getTransitionMaskWidth(const DepositionType a_depositionType, const int a_refRat) noexcept
  {
    return getHaloGhostWidth(a_depositionType, a_refRat);
  }

  /**
   * @brief Width, in finer cells, of the coarse-to-fine particle ghost mask that Transition deposition fills with.
   * @details The band width plus the one cell that a finer-width CIC or TSC cloud reaches.
   * @param[in] a_depositionType Deposition kernel.
   * @param[in] a_refRat         Refinement ratio between the level and the finer one.
   * @return The mask width, or 0 for NGP.
   */
  static int
  getTransitionGhostWidth(const DepositionType a_depositionType, const int a_refRat) noexcept
  {
    const int maskWidth = getTransitionMaskWidth(a_depositionType, a_refRat);

    return (maskWidth > 0) ? maskWidth + 1 : 0;
  }

  /**
   * @brief The coarse-to-fine particle ghost-mask widths a deposition needs between two levels.
   * @details Callers register these through AmrMesh::registerParticleDeposition before regridding. Interp and HaloNGP
   * need none; Halo needs getHaloGhostWidth; Transition needs getTransitionMaskWidth and getTransitionGhostWidth.
   * @param[in] a_depositionType       Deposition kernel.
   * @param[in] a_coarseFineDeposition Coarse-fine handling.
   * @param[in] a_refRat               Refinement ratio between the level and the finer one.
   * @return The widths (all > 0).
   */
  static std::vector<int>
  getCoarseToFineMaskWidths(const DepositionType       a_depositionType,
                            const CoarseFineDeposition a_coarseFineDeposition,
                            const int                  a_refRat) noexcept
  {
    std::vector<int> widths;

    if (a_coarseFineDeposition == CoarseFineDeposition::Halo) {
      widths.push_back(getHaloGhostWidth(a_depositionType, a_refRat));
    }
    else if (a_coarseFineDeposition == CoarseFineDeposition::Transition) {
      widths.push_back(getTransitionMaskWidth(a_depositionType, a_refRat));
      widths.push_back(getTransitionGhostWidth(a_depositionType, a_refRat));
    }

    widths.erase(std::remove(widths.begin(), widths.end(), 0), widths.end());

    return widths;
  }

  ///@}

protected:
  /**
   * @brief Per-level EB level grids.
   */
  Vector<RefCountedPtr<EBLevelGrid>> m_eblgs;

  /**
   * @brief Per-level refinement ratios.
   */
  Vector<int> m_refRat;

  /**
   * @brief Per-level grid spacing (isotropic).
   */
  Vector<Real> m_dx;

  /**
   * @brief Lower-left corner of the physical domain.
   */
  RealVect m_probLo;

  /**
   * @brief Number of ghost cells in the mesh data holders.
   */
  int m_ghost;

  /**
   * @brief Finest AMR level index.
   */
  int m_finestLevel;

  /**
   * @brief Whether define() has been called.
   */
  bool m_isDefined;

  /**
   * @brief Per-level per-patch deposition/interpolation leaves.
   */
  Vector<RefCountedPtr<LayoutData<EBParticleMesh>>> m_ebParticleMesh;

  /**
   * @brief Per-level valid+ghost -> valid Copiers (fold ghost-deposited mass into the owner patch).
   */
  Vector<Copier> m_levelCopiers;

  /**
   * @brief Per-level coarse-fine mesh operators (null on level 0).
   */
  Vector<RefCountedPtr<EBCoarseFineParticleMesh>> m_coarseFinePM;

  /**
   * @brief Outer-halo masks keyed by halo width: per-level coarse cells just outside the fine patch.
   */
  std::map<int, Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>> m_outerHaloMasks;

  /**
   * @brief Abort on an unrecognized coarse-fine deposition strategy.
   */
  void
  errorUnsupported() const
  {
    MayDay::Error("EBAMRParticleMesh: unsupported CoarseFineDeposition strategy");
  }

  /**
   * @brief Build the per-level per-patch EBParticleMesh leaves.
   */
  void
  defineEBParticleMesh()
  {
    CH_TIME("EBAMRParticleMesh::defineEBParticleMesh");

    m_ebParticleMesh.resize(1 + m_finestLevel);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const ProblemDomain&     domain = m_eblgs[lvl]->getDomain();
      const DisjointBoxLayout& dbl    = m_eblgs[lvl]->getDBL();
      const EBISLayout&        ebisl  = m_eblgs[lvl]->getEBISL();
      const DataIterator&      dit    = dbl.dataIterator();

      m_ebParticleMesh[lvl] = RefCountedPtr<LayoutData<EBParticleMesh>>(new LayoutData<EBParticleMesh>(dbl));

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        (*m_ebParticleMesh[lvl])[din].define(domain, dbl[din], ebisl[din], m_dx[lvl] * RealVect::Unit, m_probLo);
      }
    }
  }

  /**
   * @brief Build the per-level valid+ghost -> valid Copiers.
   */
  void
  defineLevelMotion()
  {
    CH_TIME("EBAMRParticleMesh::defineLevelMotion");

    m_levelCopiers.resize(1 + m_finestLevel);
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const EBLevelGrid&       eblg = *m_eblgs[lvl];
      const DisjointBoxLayout& dbl  = eblg.getDBL();

      m_levelCopiers[lvl].define(dbl, dbl, eblg.getDomain(), m_ghost * IntVect::Unit, true);
      m_levelCopiers[lvl].reverse();
    }
  }

  /**
   * @brief Build the per-level coarse-fine mesh operators (null on level 0).
   */
  void
  defineCoarseFineMotion()
  {
    CH_TIME("EBAMRParticleMesh::defineCoarseFineMotion");

    m_coarseFinePM.resize(1 + m_finestLevel);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      if (lvl > 0) {
        m_coarseFinePM[lvl] = RefCountedPtr<EBCoarseFineParticleMesh>(
          new EBCoarseFineParticleMesh(*m_eblgs[lvl - 1], *m_eblgs[lvl], m_refRat[lvl - 1], m_ghost * IntVect::Unit));
      }
      else {
        m_coarseFinePM[lvl] = RefCountedPtr<EBCoarseFineParticleMesh>(nullptr);
      }
    }
  }

  /**
   * @brief The CoarseFineDeposition::Interp deposition algorithm, parameterized by the per-patch deposit.
   * @details Resets the mesh, deposits each level's valid particles (via a_patchDeposit), folds
   * ghost-deposited mass into the owner patch (level Copier + EBAddOp), then moves the over-boundary
   * mass between levels with the coarse-fine operators. Uses only valid particles.
   * @tparam Container    ParticleContainer type (deduced).
   * @tparam PatchDeposit Callable (const EBParticleMesh&, EBCellFAB&, const ParticleSoA&, Real widthScale).
   * @param[in,out] a_meshData     AMR mesh data (RESET then filled).
   * @param[in]     a_particles    Particle container.
   * @param[in]     a_numComp      Number of mesh components written.
   * @param[in]     a_cutCellsHoldDensity True if the cut cells hold the density rather than kappa times it (mirror
   * deposition).
   * @param[in]     a_patchDeposit Per-patch deposition functor.
   */
  template <typename Container, typename PatchDeposit>
  void
  depositInterpCore(EBAMRCellData&   a_meshData,
                    const Container& a_particles,
                    const int        a_numComp,
                    const bool       a_cutCellsHoldDensity,
                    PatchDeposit     a_patchDeposit) const
  {
    CH_TIME("EBAMRParticleMesh::depositInterpCore");

    CH_assert(m_isDefined);
    CH_assert(a_meshData[0]->nComp() == a_numComp);

    DataOps::setValue(a_meshData, 0.0);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
      const bool               hasCoar = (lvl > 0);
      const DataIterator&      dit     = dbl.dataIterator();

      // 1. Deposit this level's valid particles.
      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&      din    = dit[mybox];
        const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

        EBCellFAB& rho = (*a_meshData[lvl])[din];

        a_patchDeposit(interp, rho, a_particles[lvl][din], 1.0);
      }

      // 2. Fold the mass deposited into ghost cells into the owning patch's valid region.
      a_meshData[lvl]->exchange(Interval(0, a_numComp - 1), m_levelCopiers[lvl], EBAddOp());

      // 3. Move over-refinement-boundary mass between this level and the coarser one (mesh-side).
      if (hasCoar) {
        for (int comp = 0; comp < a_numComp; comp++) {
          LevelData<EBCellFAB> coarAlias;
          LevelData<EBCellFAB> fineAlias;

          aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(comp, comp));
          aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(comp, comp));

          m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);
          m_coarseFinePM[lvl]->addInvalidCoarseToFine(fineAlias, coarAlias, a_cutCellsHoldDensity);
        }
      }
    }
  }

  /**
   * @brief Build the per-width outer-halo masks: per level, the coarse cells just outside the fine patch.
   * @details For each ghost width in [1, m_ghost] and each level with a finer level, marks the coarse
   * cells within that width of the refinement boundary (on the coarse side), excluding the
   * coarsened-fine valid region and its neighbours. The mask drives copyMaskParticles in the Halo
   * deposit. Storage-agnostic (operates on the grids only). Empty BaseFabs are undefined.
   */
  void
  defineOuterHaloMasks()
  {
    CH_TIME("EBAMRParticleMesh::defineOuterHaloMasks");

    constexpr int comp    = 0;
    constexpr int numComp = 1;

    m_outerHaloMasks.clear();

    for (int ighost = 1; ighost <= m_ghost; ighost++) {
      Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> mask(1 + m_finestLevel);

      for (int lvl = 0; lvl < m_finestLevel; lvl++) {
        const DisjointBoxLayout& grids     = m_eblgs[lvl]->getDBL();
        const DisjointBoxLayout& gridsFine = m_eblgs[lvl + 1]->getDBL();
        const ProblemDomain&     domain    = m_eblgs[lvl]->getDomain();

        DisjointBoxLayout gridsCoFi;
        coarsen(gridsCoFi, gridsFine, m_refRat[lvl]);

        mask[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(
          new LevelData<BaseFab<bool>>(grids, numComp, IntVect::Zero));

        LevelData<BaseFab<bool>>& levelMask = *mask[lvl];
        {
          const DataIterator& dit = grids.dataIterator();

          const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
          for (int mybox = 0; mybox < nbox; mybox++) {
            const DataIndex& din = dit[mybox];

            levelMask[din].setVal(false);
          }
        }

        // The coarse-grid cells surrounding the fine level, viewed from the coarsened-fine grids.
        IntVectSet halo;
        {
          const DataIterator& dit  = gridsCoFi.dataIterator();
          const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : halo)
          for (int mybox = 0; mybox < nbox; mybox++) {
            const DataIndex& din      = dit[mybox];
            const Box        coFiBox  = gridsCoFi[din];
            const Box        grownBox = grow(coFiBox, ighost) & domain;
            IntVectSet       myHalo(grownBox);
            myHalo -= coFiBox;

            NeighborIterator nit(gridsCoFi);
            for (nit.begin(din); nit.ok(); ++nit) {
              myHalo -= gridsCoFi[nit()];
            }
            halo |= myHalo;
          }
        }

        // Re-express the halo on the coarse grid via a CoFi data holder whose ghosts add onto the coarse grid.
        LevelData<FArrayBox> coFiMask(gridsCoFi, numComp, ighost * IntVect::Unit);
        LevelData<FArrayBox> coarMask(grids, numComp, IntVect::Zero);
        {
          const DataIterator& dit  = gridsFine.dataIterator();
          const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
          for (int mybox = 0; mybox < nbox; mybox++) {
            const DataIndex& din = dit[mybox];

            coFiMask[din].setVal(0.0);
          }
        }
        {
          const DataIterator& dit  = grids.dataIterator();
          const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
          for (int mybox = 0; mybox < nbox; mybox++) {
            const DataIndex& din = dit[mybox];

            coarMask[din].setVal(0.0);
          }
        }
        {
          const DataIterator& dit  = gridsFine.dataIterator();
          const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
          for (int mybox = 0; mybox < nbox; mybox++) {
            const DataIndex& din     = dit[mybox];
            const Box        region  = coFiMask[din].box();
            const IntVectSet curHalo = halo & region;

            for (IVSIterator ivsit(curHalo); ivsit.ok(); ++ivsit) {
              coFiMask[din](ivsit(), comp) = 1.0;
            }
          }
        }

        Copier copier;
        copier.ghostDefine(gridsCoFi, grids, domain, ighost * IntVect::Unit);
        coFiMask.copyTo(Interval(comp, comp), coarMask, Interval(comp, comp), copier, LDaddOp<FArrayBox>());

        // Make the boolean mask; undefine empty BaseFabs (so copyMaskParticles skips them).
        const DataIterator& dit  = grids.dataIterator();
        const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex& din      = dit[mybox];
          const Box        box      = grids[din];
          const FArrayBox& realMask = coarMask[din];

          BaseFab<bool>& boolMask = levelMask[din];

          bool emptyMask = true;
          for (BoxIterator bit(box); bit.ok(); ++bit) {
            if (realMask(bit(), comp) > 0.0) {
              boolMask(bit(), comp) = true;
              emptyMask             = false;
            }
          }
          if (emptyMask) {
            boolMask.clear();
          }
        }
      }

      mask[m_finestLevel] = RefCountedPtr<LevelData<BaseFab<bool>>>(nullptr);
      m_outerHaloMasks.emplace(ighost, mask);
    }
  }

  /**
   * @brief Assemble a per-level coarse-to-fine particle ghost mask from the realm's masks.
   * @details Entry lvl, for each level with a finer level, is entry lvl of the realm's mask of width
   * a_width(m_refRat[lvl]); it is null where that width is 0 and on the finest level.
   * @tparam MaskGetter Callable (int width) -> const AMRParticleGhostMask&.
   * @tparam WidthFn    Callable (int refRat) -> int.
   * @param[in] a_coarToFineMasks Realm mask getter.
   * @param[in] a_width           Mask width of a level pair, from its refinement ratio.
   * @return The per-level mask.
   */
  template <typename MaskGetter, typename WidthFn>
  AMRParticleGhostMask
  getCoarToFineMask(const MaskGetter& a_coarToFineMasks, const WidthFn& a_width) const
  {
    AMRParticleGhostMask mask(1 + m_finestLevel);

    for (int lvl = 0; lvl < m_finestLevel; lvl++) {
      const int width = a_width(m_refRat[lvl]);

      if (width > 0) {
        mask[lvl] = a_coarToFineMasks(width)[lvl];
      }
    }

    return mask;
  }

  /**
   * @brief Copy the mask-holder particles into the finer patches as coarse ghosts.
   * @details When NGP is forced in cut cells, particles in cut cells deposit by NGP on their own level only, so they
   * are not copied.
   * @tparam Container ParticleContainer type (deduced).
   * @param[in,out] a_particles     Particle container.
   * @param[in]     a_coarToFine    Per-level coarse-to-fine particle ghost mask.
   * @param[in]     a_forceIrregNGP Force NGP in cut cells.
   */
  template <typename Container>
  void
  fillCoarseGhosts(Container& a_particles, const AMRParticleGhostMask& a_coarToFine, const bool a_forceIrregNGP) const
  {
    CH_TIME("EBAMRParticleMesh::fillCoarseGhosts");

    a_particles.fillCoarseGhostParticlesFromMask(a_coarToFine,
                                                 [&](const int a_lvl, const DataIndex& a_din, const IntVect& a_iv) {
                                                   return !(a_forceIrregNGP &&
                                                            m_eblgs[a_lvl]->getEBISL()[a_din].isIrregular(a_iv));
                                                 });
  }

  /**
   * @brief The CoarseFineDeposition::Halo deposition algorithm.
   * @details Like Interp for the on-level deposit and the fine-ghost transfer to the coarse level, but instead of
   * interpolating the coarse mass under the fine grid, the coarse-side halo particles (width-1 outer-halo mask) are
   * copied into the finer patches as coarse ghosts and deposit there with the coarse width, clipped to the valid box.
   * Their coarse-level deposit is unchanged; its part under the finer grid is invalid.
   * @tparam NCOMP        Number of mesh components (compile-time).
   * @tparam Container    ParticleContainer type (deduced).
   * @tparam MaskGetter   Callable (int width) -> const AMRParticleGhostMask&.
   * @tparam PatchDeposit Callable (const EBParticleMesh&, EBCellFAB&, const ParticleSoA&, Real widthScale).
   * @tparam PatchGhosts  Same signature, depositing only the coarse ghosts.
   * @param[in,out] a_meshData        AMR mesh data (RESET then filled).
   * @param[in]     a_particles       Particle container (its mask holder and ghosts are filled then cleared).
   * @param[in]     a_depositionType  Deposition kernel.
   * @param[in]     a_forceIrregNGP   Force NGP in cut cells.
   * @param[in]     a_coarToFineMasks Realm coarse-to-fine mask getter.
   * @param[in]     a_patchDeposit    Per-patch deposition functor.
   * @param[in]     a_patchGhosts     Per-patch coarse-ghost deposition functor.
   */
  template <int NCOMP, typename Container, typename MaskGetter, typename PatchDeposit, typename PatchGhosts>
  void
  depositHaloCore(EBAMRCellData&       a_meshData,
                  const Container&     a_particles,
                  const DepositionType a_depositionType,
                  const bool           a_forceIrregNGP,
                  const MaskGetter&    a_coarToFineMasks,
                  PatchDeposit         a_patchDeposit,
                  PatchGhosts          a_patchGhosts) const
  {
    CH_TIME("EBAMRParticleMesh::depositHaloCore");

    CH_assert(m_isDefined);
    CH_assert(a_meshData[0]->nComp() == NCOMP);

    DataOps::setValue(a_meshData, 0.0);

    // Copy the coarse-side halo particles (width 1) into the container's mask holder, and from there into the finer
    // patches their clouds reach.
    constexpr int coarseMaskWidth = 1;
    Container&    particles       = const_cast<Container&>(a_particles);

    particles.copyMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

    this->fillCoarseGhosts(particles,
                           this->getCoarToFineMask(a_coarToFineMasks,
                                                   [&](const int a_refRat) {
                                                     return getHaloGhostWidth(a_depositionType, a_refRat);
                                                   }),
                           a_forceIrregNGP);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
      const bool               hasCoar = (lvl > 0);

      const DataIterator& dit  = dbl.dataIterator();
      const int           nbox = dit.size();

      // 1. Deposit this level's own particles.
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&      din    = dit[mybox];
        const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

        EBCellFAB& rho = (*a_meshData[lvl])[din];

        a_patchDeposit(interp, rho, a_particles[lvl][din], 1.0);
      }

      // 2. Fold ghost-deposited mass into the owner patch.
      a_meshData[lvl]->exchange(Interval(0, NCOMP - 1), m_levelCopiers[lvl], EBAddOp());

      // 3. Mass deposited into fine ghosts across the refinement boundary belongs to the coarse level.
      if (hasCoar) {
        for (int comp = 0; comp < NCOMP; comp++) {
          LevelData<EBCellFAB> coarAlias;
          LevelData<EBCellFAB> fineAlias;

          aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(comp, comp));
          aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(comp, comp));

          m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);
        }
      }

      // 4. The coarse halo particles deposit into the valid cells of this level with the coarse width.
      if (hasCoar) {
        const Real widthScale = static_cast<Real>(m_refRat[lvl - 1]);

#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex&      din    = dit[mybox];
          const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

          EBCellFAB& rho = (*a_meshData[lvl])[din];

          a_patchGhosts(interp, rho, a_particles[lvl][din], widthScale);
        }
      }
    }

    particles.clearGhostParticles(GhostType::Coarse);
    particles.clearMaskParticles();
  }

  /**
   * @brief The CoarseFineDeposition::HaloNGP deposition algorithm.
   * @details Moves the coarse-side halo particles (width-1 mask) out of the valid holder, deposits
   * the remaining (non-halo) particles with the requested kernel and the halo particles with NGP onto
   * the same level mesh, folds fine-ghost mass to the coarse level, then transfers the halo particles
   * back to the valid holder. NGP keeps the near-boundary clouds from spreading over the refinement
   * boundary.
   * @tparam NCOMP           Number of mesh components (compile-time).
   * @tparam Container       ParticleContainer type (deduced).
   * @tparam PatchDeposit    Per-patch deposit with the requested kernel.
   * @tparam PatchDepositNGP Per-patch deposit forced to NGP.
   * @param[in,out] a_meshData     AMR mesh data (RESET then filled).
   * @param[in]     a_particles    Particle container (halo particles are moved out then back).
   * @param[in]     a_patch        Per-patch deposit functor (requested kernel).
   * @param[in]     a_patchNGP     Per-patch deposit functor (NGP).
   */
  template <int NCOMP, typename Container, typename PatchDeposit, typename PatchDepositNGP>
  void
  depositHaloNGPCore(EBAMRCellData&   a_meshData,
                     const Container& a_particles,
                     PatchDeposit     a_patch,
                     PatchDepositNGP  a_patchNGP) const
  {
    CH_TIME("EBAMRParticleMesh::depositHaloNGPCore");

    CH_assert(m_isDefined);
    CH_assert(a_meshData[0]->nComp() == NCOMP);

    DataOps::setValue(a_meshData, 0.0);

    // Move the coarse-side halo particles (width 1) out of the valid holder into the mask holder.
    constexpr int coarseMaskWidth = 1;
    Container&    particles       = const_cast<Container&>(a_particles);

    particles.transferMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
      const bool               hasCoar = (lvl > 0);

      // 1. Deposit the non-halo particles with the requested kernel, the halo particles with NGP.
      const DataIterator& dit  = dbl.dataIterator();
      const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&      din    = dit[mybox];
        const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

        EBCellFAB& rho = (*a_meshData[lvl])[din];

        a_patch(interp, rho, a_particles[lvl][din], 1.0);
        a_patchNGP(interp, rho, (*a_particles.getMaskParticles()[lvl])[din], 1.0);
      }

      // 2. Fold ghost-deposited mass into the owner patch.
      a_meshData[lvl]->exchange(Interval(0, NCOMP - 1), m_levelCopiers[lvl], EBAddOp());

      // 3. Mass deposited into fine ghosts across the refinement boundary belongs to the coarse level.
      if (hasCoar) {
        for (int comp = 0; comp < NCOMP; comp++) {
          LevelData<EBCellFAB> coarAlias;
          LevelData<EBCellFAB> fineAlias;

          aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(comp, comp));
          aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(comp, comp));

          m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);
        }
      }
    }

    // Return the halo particles to the valid holder, each where it was.
    particles.returnMaskParticles();
  }

  /**
   * @brief The CoarseFineDeposition::Transition deposition algorithm.
   * @details NGP short-circuits to the Interp result. Otherwise, the coarser particles within the transition band of a
   * finer grid deposit with the finer-level width: on their own level, into the same data as the other particles, and
   * as coarse ghosts in the finer patches their clouds reach, clipped to the valid box. The part of the coarse-level
   * deposit under the finer grid is invalid, so the composite mass is conserved.
   * @tparam NCOMP        Number of mesh components (compile-time).
   * @tparam Container    ParticleContainer type (deduced).
   * @tparam MaskGetter   Callable (int width) -> const AMRParticleGhostMask&.
   * @tparam PatchDeposit Per-patch deposit with the requested kernel.
   * @tparam PatchGhosts  Per-patch coarse-ghost deposit with the requested kernel.
   * @param[in,out] a_meshData        AMR mesh data (RESET then filled).
   * @param[in]     a_particles       Particle container.
   * @param[in]     a_depositionType  Deposition kernel.
   * @param[in]     a_forceIrregNGP   Force NGP in cut cells.
   * @param[in]     a_cutCellsHoldDensity True if the cut cells hold the density (mirror deposition); passed to the
   * Interp fallback for NGP.
   * @param[in]     a_coarToFineMasks Realm coarse-to-fine mask getter.
   * @param[in]     a_patch           Per-patch deposit functor (requested kernel).
   * @param[in]     a_patchGhosts     Per-patch coarse-ghost deposit functor.
   */
  template <int NCOMP, typename Container, typename MaskGetter, typename PatchDeposit, typename PatchGhosts>
  void
  depositTransitionCore(EBAMRCellData&       a_meshData,
                        const Container&     a_particles,
                        const DepositionType a_depositionType,
                        const bool           a_forceIrregNGP,
                        const bool           a_cutCellsHoldDensity,
                        const MaskGetter&    a_coarToFineMasks,
                        PatchDeposit         a_patch,
                        PatchGhosts          a_patchGhosts) const
  {
    CH_TIME("EBAMRParticleMesh::depositTransitionCore");

    CH_assert(m_isDefined);
    CH_assert(a_meshData[0]->nComp() == NCOMP);

    if (a_depositionType == DepositionType::NGP) {
      this->depositInterpCore(a_meshData, a_particles, NCOMP, a_cutCellsHoldDensity, a_patch);

      return;
    }

    DataOps::setValue(a_meshData, 0.0);

    // Move the transition particles into the mask holder, and copy them into the finer patches their finer-width clouds
    // reach.
    Container& particles = const_cast<Container&>(a_particles);

    particles.transferMaskParticlesCoarToFine(this->getCoarToFineMask(a_coarToFineMasks, [&](const int a_refRat) {
      return getTransitionMaskWidth(a_depositionType, a_refRat);
    }));

    this->fillCoarseGhosts(particles,
                           this->getCoarToFineMask(a_coarToFineMasks,
                                                   [&](const int a_refRat) {
                                                     return getTransitionGhostWidth(a_depositionType, a_refRat);
                                                   }),
                           a_forceIrregNGP);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
      const bool               hasCoar = (lvl > 0);
      const bool               hasFine = (lvl < m_finestLevel);

      const DataIterator& dit  = dbl.dataIterator();
      const int           nbox = dit.size();

      // 1. Deposit this level's own particles, and its transition particles with the finer-level width.
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&      din    = dit[mybox];
        const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

        EBCellFAB& rho = (*a_meshData[lvl])[din];

        a_patch(interp, rho, a_particles[lvl][din], 1.0);

        if (hasFine) {
          a_patch(interp, rho, (*a_particles.getMaskParticles()[lvl])[din], 1.0 / m_refRat[lvl]);
        }
      }

      // 2. Fold ghost-deposited mass into the owner patch; fine-ghost mass over the boundary goes to the coarse level.
      a_meshData[lvl]->exchange(Interval(0, NCOMP - 1), m_levelCopiers[lvl], EBAddOp());

      if (hasCoar) {
        for (int comp = 0; comp < NCOMP; comp++) {
          LevelData<EBCellFAB> coarAlias;
          LevelData<EBCellFAB> fineAlias;

          aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(comp, comp));
          aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(comp, comp));

          m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);
        }

        // 3. The coarser transition particles deposit into the valid cells of this level with this level's width.
#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex&      din    = dit[mybox];
          const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];

          EBCellFAB& rho = (*a_meshData[lvl])[din];

          a_patchGhosts(interp, rho, a_particles[lvl][din], 1.0);
        }
      }
    }

    // Return the transition particles to the valid holder, each where it was.
    particles.clearGhostParticles(GhostType::Coarse);
    particles.returnMaskParticles();
  }
};

#include <CD_NamespaceFooter.H>

#endif
