/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_EBAMRParticleMesh.H
  @brief  Declaration of an AMR driver that deposits/interpolates ParticleContainer particles.
  @author Robert Marskar
*/

#ifndef CD_EBAMRPARTICLEMESH_H
#define CD_EBAMRPARTICLEMESH_H

// Std includes
#include <cstddef>
#include <map>

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
#include <CD_DepositionType.H>
#include <CD_CoarseFineDeposition.H>
#include <CD_EBCoarseFineParticleMesh.H>
#include <CD_EBParticleMesh.H>
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>

/**
  @brief AMR driver that deposits/interpolates ParticleContainer particles across the hierarchy.
  @details
  EBAMRParticleMesh deposits the particles of a ParticleContainer onto AMR mesh data and interpolates
  mesh data back onto particle columns, by driving the per-patch EBParticleMesh leaf on every box of
  every level. define() builds the per-level leaves (on the valid grids and the refined-coarse grids),
  the level valid+ghost -> valid Copiers, the EBCoarseFineParticleMesh operators, and the halo/transition
  masks; it must be called before any deposit or interpolate.

  <h3>Interpolation</h3>
  interpolate() and interpolateWeight() are purely per-level/per-patch: the mesh is assumed to be
  valid+ghosted by the caller, so no coarse-fine coupling is needed.

  <h3>Deposition</h3>
  deposit(), depositWeight(), and depositGathered() conserve mass across refinement boundaries, and
  select how the coarse-fine overlap is handled through @c CoarseFineDeposition:
  - <b>Interp</b>     -- move the over-boundary mass between levels with the coarse-fine mesh operators.
  - <b>Halo</b>       -- re-deposit the coarse-side halo particles on the refined-coarse grid (via the
                         container's mask holder and the FiCo buffer) and add them to the fine level.
  - <b>HaloNGP</b>    -- deposit the near-boundary halo particles with NGP so their clouds do not spread
                         across the boundary (no refined-coarse buffer is needed).
  - <b>Transition</b> -- deposit the coarse-side transition particles on the refined-coarse grid at the
                         fine-grid width, distributing the result to both levels.
  Halo, HaloNGP, and Transition require at least one ghost cell.
*/
class EBAMRParticleMesh
{
public:
  /**
    @name Construction and lifetime
  */
  ///@{

  /**
    @brief Default constructor. Leaves the object undefined; call define().
  */
  EBAMRParticleMesh() : m_ghost(0), m_finestLevel(-1), m_isDefined(false)
  {}

  /**
    @brief Destructor.
  */
  ~EBAMRParticleMesh() = default;

  /**
    @brief Copy construction is deleted (shared RefCountedPtr state).
  */
  EBAMRParticleMesh(const EBAMRParticleMesh&) = delete;

  /**
    @brief Copy assignment is deleted.
    @return Reference to this.
  */
  EBAMRParticleMesh&
  operator=(const EBAMRParticleMesh&) = delete;

  /**
    @brief Define the per-level leaves, coarse-fine mesh operators, and halo/transition masks.
    @param[in] a_eblgs       Per-level EB level grids.
    @param[in] a_refRat      Per-level refinement ratio (entry l is the ratio between level l and l+1).
    @param[in] a_dx          Per-level (isotropic) grid spacing.
    @param[in] a_probLo      Lower-left corner of the physical domain.
    @param[in] a_ghost       Number of ghost cells in the mesh data holders.
    @param[in] a_finestLevel Finest AMR level index.
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
    this->defineTransitionMasks();

    m_isDefined = true;
  }

  ///@}

  /**
    @name Interpolation
  */
  ///@{

  /**
    @brief Interpolate a mesh field onto one or more payload columns, across all levels.
    @details Per level, per patch: calls EBParticleMesh::interpolate. The mesh is assumed to be
    valid+ghosted by the caller (no coarse-fine coupling is needed for interpolation).
    @tparam Members Pack of payload data-member pointers selecting the target column(s).
    @tparam P       Payload type (deduced).
    @tparam Traits  Column descriptor (deduced).
    @param[in,out] a_particles     Particles whose target columns are overwritten.
    @param[in]     a_meshData      AMR mesh data (>= sizeof...(Members) components).
    @param[in]     a_interpType    Interpolation kernel (NGP/CIC/TSC).
    @param[in]     a_forceIrregNGP Force NGP in cut cells.
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
      const DisjointBoxLayout& dbl  = m_eblgs[lvl]->getDBL();
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();
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
    @brief Interpolate a mesh field onto the container-owned weight column, across all levels.
    @details Per level, per patch: calls EBParticleMesh::interpolateWeight. The mesh is assumed to
    be valid+ghosted by the caller (no coarse-fine coupling is needed for interpolation). This is the
    interpolation mirror of depositWeight.
    @tparam P      Payload type (deduced).
    @tparam Traits Column descriptor (deduced).
    @param[in,out] a_particles     Particles whose weight column is overwritten.
    @param[in]     a_meshData      AMR mesh data (>= 1 component).
    @param[in]     a_interpType    Interpolation kernel (NGP/CIC/TSC).
    @param[in]     a_forceIrregNGP Force NGP in cut cells.
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
      const DisjointBoxLayout& dbl  = m_eblgs[lvl]->getDBL();
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();
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
    @name Deposition (all coarse-fine strategies)
  */
  ///@{

  /**
    @brief Deposit the container-owned weight onto the AMR mesh.
    @details Dispatches on the coarse-fine strategy (Interp / Halo / HaloNGP / Transition); see the
    class documentation. Halo/HaloNGP/Transition need ghost >= 1.
    @tparam P      Payload type (deduced).
    @tparam Traits Column descriptor (deduced).
    @param[in,out] a_meshData             AMR mesh data (one component; RESET then filled).
    @param[in]     a_particles            Particles to deposit.
    @param[in]     a_depositionType       Deposition kernel (NGP/CIC/TSC).
    @param[in]     a_coarseFineDeposition Coarse-fine handling (Interp / Halo / HaloNGP / Transition).
    @param[in]     a_forceIrregNGP        Force NGP in cut cells.
  */
  template <typename P, typename Traits>
  void
  depositWeight(EBAMRCellData&                      a_meshData,
                const ParticleContainer<P, Traits>& a_particles,
                const DepositionType                a_depositionType,
                const CoarseFineDeposition          a_coarseFineDeposition,
                const bool                          a_forceIrregNGP)
  {
    CH_TIME("EBAMRParticleMesh::depositWeight");

    CH_assert(m_isDefined);

    const auto patch =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.depositWeight(a_rho, a_soa, a_depositionType, a_widthScale, a_forceIrregNGP);
      };
    const auto patchNGP =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.depositWeight(a_rho, a_soa, DepositionType::NGP, a_widthScale, a_forceIrregNGP);
      };

    switch (a_coarseFineDeposition) {
    case CoarseFineDeposition::Interp: {
      this->depositInterpCore(a_meshData, a_particles, 1, patch);
      break;
    }
    case CoarseFineDeposition::Halo: {
      this->template depositHaloCore<1>(a_meshData, a_particles, patch);
      break;
    }
    case CoarseFineDeposition::HaloNGP: {
      this->template depositHaloNGPCore<1>(a_meshData, a_particles, patch, patchNGP);
      break;
    }
    case CoarseFineDeposition::Transition: {
      this->template depositTransitionCore<1>(a_meshData, a_particles, a_depositionType, patch);
      break;
    }
    default: {
      this->errorUnsupported();
      break;
    }
    }
  }

  /**
    @brief Deposit a custom per-particle scalar (computed by a_cellGather) onto the AMR mesh.
    @details Generalizes depositWeight() to an arbitrary per-particle value (a product of columns such as
    weight*mobility), dispatching on the coarse-fine strategy exactly like depositWeight(). Deposits one
    mesh component.
    @tparam P          Payload type (deduced).
    @tparam Traits     Column descriptor (deduced).
    @tparam CellGather Callable (const ParticleSoA<P, Traits>& a_leaf, std::size_t a_i) -> Real.
    @param[in,out] a_meshData             AMR mesh data (one component; RESET then filled).
    @param[in]     a_particles            Particles to deposit.
    @param[in]     a_depositionType       Deposition kernel (NGP/CIC/TSC).
    @param[in]     a_coarseFineDeposition Coarse-fine handling (Interp / Halo / HaloNGP / Transition).
    @param[in]     a_forceIrregNGP        Force NGP in cut cells.
    @param[in]     a_cellGather           Per-particle value gatherer (leaf, index) -> Real.
  */
  template <typename P, typename Traits, typename CellGather>
  void
  depositGathered(EBAMRCellData&                      a_meshData,
                  const ParticleContainer<P, Traits>& a_particles,
                  const DepositionType                a_depositionType,
                  const CoarseFineDeposition          a_coarseFineDeposition,
                  const bool                          a_forceIrregNGP,
                  CellGather                          a_cellGather)
  {
    CH_TIME("EBAMRParticleMesh::depositGathered");

    CH_assert(m_isDefined);

    const auto patch =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.depositGathered(a_rho,
                                 0,
                                 a_soa,
                                 a_depositionType,
                                 a_widthScale,
                                 a_forceIrregNGP,
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
                                 a_forceIrregNGP,
                                 [&](const std::size_t a_i, Real* a_out) {
                                   a_out[0] = a_cellGather(a_soa, a_i);
                                 });
      };

    switch (a_coarseFineDeposition) {
    case CoarseFineDeposition::Interp: {
      this->depositInterpCore(a_meshData, a_particles, 1, patch);
      break;
    }
    case CoarseFineDeposition::Halo: {
      this->template depositHaloCore<1>(a_meshData, a_particles, patch);
      break;
    }
    case CoarseFineDeposition::HaloNGP: {
      this->template depositHaloNGPCore<1>(a_meshData, a_particles, patch, patchNGP);
      break;
    }
    case CoarseFineDeposition::Transition: {
      this->template depositTransitionCore<1>(a_meshData, a_particles, a_depositionType, patch);
      break;
    }
    default: {
      this->errorUnsupported();
      break;
    }
    }
  }

  /**
    @brief Deposit one or more payload columns onto the AMR mesh.
    @details See depositWeight for the coarse-fine dispatch.
    @tparam Members Pack of payload data-member pointers selecting the column(s) to deposit.
    @tparam P       Payload type (deduced).
    @tparam Traits  Column descriptor (deduced).
    @param[in,out] a_meshData             AMR mesh data (sizeof...(Members) components; RESET then filled).
    @param[in]     a_particles            Particles to deposit.
    @param[in]     a_depositionType       Deposition kernel (NGP/CIC/TSC).
    @param[in]     a_coarseFineDeposition Coarse-fine handling (Interp / Halo / HaloNGP / Transition).
    @param[in]     a_forceIrregNGP        Force NGP in cut cells.
  */
  template <auto... Members, typename P, typename Traits>
  void
  deposit(EBAMRCellData&                      a_meshData,
          const ParticleContainer<P, Traits>& a_particles,
          const DepositionType                a_depositionType,
          const CoarseFineDeposition          a_coarseFineDeposition,
          const bool                          a_forceIrregNGP)
  {
    CH_TIME("EBAMRParticleMesh::deposit");

    CH_assert(m_isDefined);

    const auto patch =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.template deposit<Members...>(a_rho, a_soa, a_depositionType, a_widthScale, a_forceIrregNGP);
      };
    const auto patchNGP =
      [&](const EBParticleMesh& a_interp, EBCellFAB& a_rho, const auto& a_soa, const Real a_widthScale) {
        a_interp.template deposit<Members...>(a_rho, a_soa, DepositionType::NGP, a_widthScale, a_forceIrregNGP);
      };

    switch (a_coarseFineDeposition) {
    case CoarseFineDeposition::Interp: {
      this->depositInterpCore(a_meshData, a_particles, static_cast<int>(sizeof...(Members)), patch);
      break;
    }
    case CoarseFineDeposition::Halo: {
      this->template depositHaloCore<static_cast<int>(sizeof...(Members))>(a_meshData, a_particles, patch);
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
                                                                                 patch);
      break;
    }
    default: {
      this->errorUnsupported();
      break;
    }
    }
  }

  ///@}

protected:
  /**
    @brief Per-level EB level grids.
  */
  Vector<RefCountedPtr<EBLevelGrid>> m_eblgs;

  /**
    @brief Per-level refinement ratios.
  */
  Vector<int> m_refRat;

  /**
    @brief Per-level grid spacing (isotropic).
  */
  Vector<Real> m_dx;

  /**
    @brief Lower-left corner of the physical domain.
  */
  RealVect m_probLo;

  /**
    @brief Number of ghost cells in the mesh data holders.
  */
  int m_ghost;

  /**
    @brief Finest AMR level index.
  */
  int m_finestLevel;

  /**
    @brief Whether define() has been called.
  */
  bool m_isDefined;

  /**
    @brief Per-level per-patch deposition/interpolation leaves.
  */
  Vector<RefCountedPtr<LayoutData<EBParticleMesh>>> m_ebParticleMesh;

  /**
    @brief Per-level valid+ghost -> valid Copiers (fold ghost-deposited mass into the owner patch).
  */
  Vector<Copier> m_levelCopiers;

  /**
    @brief Per-level coarse-fine mesh operators (null on level 0).
  */
  Vector<RefCountedPtr<EBCoarseFineParticleMesh>> m_coarseFinePM;

  /**
    @brief Per-level per-patch leaves on the refined-coarse (FiCo) grid (null on level 0).
  */
  Vector<RefCountedPtr<LayoutData<EBParticleMesh>>> m_ebParticleMeshFiCo;

  /**
    @brief Outer-halo masks keyed by halo width: per-level coarse cells just outside the fine patch.
  */
  std::map<int, Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>> m_outerHaloMasks;

  /**
    @brief Transition masks keyed by width: per-level fine-side transition band, on the refined-coarse grid.
  */
  std::map<int, Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>> m_transitionMasks;

  /**
    @brief Abort on an unrecognized coarse-fine deposition strategy.
  */
  void
  errorUnsupported() const
  {
    MayDay::Error("EBAMRParticleMesh: unsupported CoarseFineDeposition strategy");
  }

  /**
    @brief Build the per-level per-patch EBParticleMesh leaves (valid grids + refined-coarse grids).
  */
  void
  defineEBParticleMesh()
  {
    CH_TIME("EBAMRParticleMesh::defineEBParticleMesh");

    m_ebParticleMesh.resize(1 + m_finestLevel);
    m_ebParticleMeshFiCo.resize(1 + m_finestLevel);
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const ProblemDomain&     domain = m_eblgs[lvl]->getDomain();
      const DisjointBoxLayout& dbl    = m_eblgs[lvl]->getDBL();
      const EBISLayout&        ebisl  = m_eblgs[lvl]->getEBISL();

      m_ebParticleMesh[lvl]    = RefCountedPtr<LayoutData<EBParticleMesh>>(new LayoutData<EBParticleMesh>(dbl));
      const DataIterator& dit  = dbl.dataIterator();
      const int           nbox = dit.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];
        (*m_ebParticleMesh[lvl])[din].define(domain, dbl[din], ebisl[din], m_dx[lvl] * RealVect::Unit, m_probLo);
      }

      // Refined-coarse leaves: deposit coarse-level halo particles at fine resolution (dx of this level).
      if (lvl > 0) {
        const EBLevelGrid&       eblgFiCo   = m_coarseFinePM[lvl]->getEblgFiCo();
        const ProblemDomain&     domainFiCo = eblgFiCo.getDomain();
        const DisjointBoxLayout& dblFiCo    = eblgFiCo.getDBL();
        const EBISLayout&        ebislFiCo  = eblgFiCo.getEBISL();

        m_ebParticleMeshFiCo[lvl] = RefCountedPtr<LayoutData<EBParticleMesh>>(new LayoutData<EBParticleMesh>(dblFiCo));
        const DataIterator& ditFiCo  = dblFiCo.dataIterator();
        const int           nboxFiCo = ditFiCo.size();
#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nboxFiCo; mybox++) {
          const DataIndex& din = ditFiCo[mybox];
          (*m_ebParticleMeshFiCo[lvl])[din].define(domainFiCo,
                                                   dblFiCo[din],
                                                   ebislFiCo[din],
                                                   m_dx[lvl] * RealVect::Unit,
                                                   m_probLo);
        }
      }
      else {
        m_ebParticleMeshFiCo[lvl] = RefCountedPtr<LayoutData<EBParticleMesh>>(nullptr);
      }
    }
  }

  /**
    @brief Build the per-level valid+ghost -> valid Copiers.
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
    @brief Build the per-level coarse-fine mesh operators (null on level 0).
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
    @brief The CoarseFineDeposition::Interp deposition algorithm, parameterized by the per-patch deposit.
    @details Resets the mesh, deposits each level's valid particles (via a_patchDeposit), folds
    ghost-deposited mass into the owner patch (level Copier + EBAddOp), then moves the over-boundary
    mass between levels with the coarse-fine operators. Uses only valid particles.
    @tparam Container    ParticleContainer type (deduced).
    @tparam PatchDeposit Callable (const EBParticleMesh&, EBCellFAB&, const ParticleSoA&, Real widthScale).
    @param[in,out] a_meshData     AMR mesh data (RESET then filled).
    @param[in]     a_particles    Particle container.
    @param[in]     a_numComp      Number of mesh components written.
    @param[in]     a_patchDeposit Per-patch deposition functor.
  */
  template <typename Container, typename PatchDeposit>
  void
  depositInterpCore(EBAMRCellData&   a_meshData,
                    const Container& a_particles,
                    const int        a_numComp,
                    PatchDeposit     a_patchDeposit) const
  {
    CH_TIME("EBAMRParticleMesh::depositInterpCore");

    CH_assert(m_isDefined);
    CH_assert(a_meshData[0]->nComp() == a_numComp);

    DataOps::setValue(a_meshData, 0.0);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
      const bool               hasCoar = (lvl > 0);

      // 1. Deposit this level's valid particles.
      const DataIterator& dit  = dbl.dataIterator();
      const int           nbox = dit.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&      din    = dit[mybox];
        const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];
        EBCellFAB&            rho    = (*a_meshData[lvl])[din];
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
          m_coarseFinePM[lvl]->addInvalidCoarseToFine(fineAlias, coarAlias);
        }
      }
    }
  }

  /**
    @brief Build the per-width outer-halo masks: per level, the coarse cells just outside the fine patch.
    @details For each ghost width in [1, m_ghost] and each level with a finer level, marks the coarse
    cells within that width of the refinement boundary (on the coarse side), excluding the
    coarsened-fine valid region and its neighbours. The mask drives copyMaskParticles in the Halo
    deposit. Storage-agnostic (operates on the grids only). Empty BaseFabs are undefined.
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
          const DataIterator& dit  = grids.dataIterator();
          const int           nbox = dit.size();
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
          // serial: accumulates into the shared IntVectSet halo (not indexed by box).
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
          BaseFab<bool>&   boolMask = levelMask[din];
          const FArrayBox& realMask = coarMask[din];

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
    @brief The CoarseFineDeposition::Halo deposition algorithm.
    @details Like Interp for the on-level deposit + fine-ghost->coarse transfer, but instead of
    interpolating the coarse mass under the fine grid, it re-deposits the coarse-side halo particles
    (selected by the width-1 outer-halo mask, copied into the container mask holder) on the
    refined-coarse grid with widthScale = refinement ratio, and adds that to the fine level.
    @tparam NCOMP        Number of mesh components (compile-time; selects the FiCo buffer).
    @tparam Container    ParticleContainer type (deduced).
    @tparam PatchDeposit Callable (const EBParticleMesh&, EBCellFAB&, const ParticleSoA&, Real widthScale).
    @param[in,out] a_meshData     AMR mesh data (RESET then filled).
    @param[in]     a_particles    Particle container (its mask holder is filled then cleared).
    @param[in]     a_patchDeposit Per-patch deposition functor.
  */
  template <int NCOMP, typename Container, typename PatchDeposit>
  void
  depositHaloCore(EBAMRCellData& a_meshData, const Container& a_particles, PatchDeposit a_patchDeposit) const
  {
    CH_TIME("EBAMRParticleMesh::depositHaloCore");

    CH_assert(m_isDefined);
    CH_assert(a_meshData[0]->nComp() == NCOMP);

    DataOps::setValue(a_meshData, 0.0);

    // Copy the coarse-side halo particles (width 1) into the container's mask holder.
    constexpr int coarseMaskWidth = 1;
    Container&    particles       = const_cast<Container&>(a_particles);
    particles.copyMaskParticles(m_outerHaloMasks.at(coarseMaskWidth));

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl     = m_eblgs[lvl]->getDBL();
      const bool               hasCoar = (lvl > 0);

      // 1. Deposit this level's valid particles.
      const DataIterator& dit  = dbl.dataIterator();
      const int           nbox = dit.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&      din    = dit[mybox];
        const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];
        EBCellFAB&            rho    = (*a_meshData[lvl])[din];
        a_patchDeposit(interp, rho, a_particles[lvl][din], 1.0);
      }

      // 2. Fold ghost-deposited mass into the owner patch.
      a_meshData[lvl]->exchange(Interval(0, NCOMP - 1), m_levelCopiers[lvl], EBAddOp());

      // 3. Mass deposited into fine ghosts across the refinement boundary belongs to the coarse level.
      if (hasCoar) {
        m_coarseFinePM[lvl]->addFineGhostsToCoarse(*a_meshData[lvl - 1], *a_meshData[lvl]);
      }

      // 4. Re-deposit the coarse halo particles on the refined-coarse grid and add to the fine level.
      if (hasCoar) {
        const auto&              coarHalo   = a_particles.getMaskParticles()[lvl - 1];
        const Real               widthScale = static_cast<Real>(m_refRat[lvl - 1]);
        const EBLevelGrid&       eblgFiCo   = m_coarseFinePM[lvl]->getEblgFiCo();
        const DisjointBoxLayout& dblFiCo    = eblgFiCo.getDBL();
        LevelData<EBCellFAB>&    bufferFiCo = m_coarseFinePM[lvl]->template getBufferFiCo<NCOMP>();

        const DataIterator& ditFiCo  = dblFiCo.dataIterator();
        const int           nboxFiCo = ditFiCo.size();
#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nboxFiCo; mybox++) {
          const DataIndex& din      = ditFiCo[mybox];
          EBCellFAB&       dataFiCo = bufferFiCo[din];
          dataFiCo.setVal(0.0);
          const EBParticleMesh& interp = (*m_ebParticleMeshFiCo[lvl])[din];
          a_patchDeposit(interp, dataFiCo, (*coarHalo)[din], widthScale);
        }

        for (int comp = 0; comp < NCOMP; comp++) {
          LevelData<EBCellFAB> meshAlias;
          LevelData<EBCellFAB> bufAlias;
          aliasLevelData<EBCellFAB>(meshAlias, &(*a_meshData[lvl]), Interval(comp, comp));
          aliasLevelData<EBCellFAB>(bufAlias, &bufferFiCo, Interval(comp, comp));
          m_coarseFinePM[lvl]->addFiCoDataToFine(meshAlias, bufAlias);
        }
      }
    }

    particles.clearMaskParticles();
  }

  /**
    @brief The CoarseFineDeposition::HaloNGP deposition algorithm.
    @details Moves the coarse-side halo particles (width-1 mask) out of the valid holder, deposits
    the remaining (non-halo) particles with the requested kernel and the halo particles with NGP onto
    the same level mesh, folds fine-ghost mass to the coarse level, then transfers the halo particles
    back to the valid holder. NGP keeps the near-boundary clouds from spreading over the refinement
    boundary, so no refined-coarse buffer is needed.
    @tparam NCOMP           Number of mesh components (compile-time).
    @tparam Container       ParticleContainer type (deduced).
    @tparam PatchDeposit    Per-patch deposit with the requested kernel.
    @tparam PatchDepositNGP Per-patch deposit forced to NGP.
    @param[in,out] a_meshData     AMR mesh data (RESET then filled).
    @param[in]     a_particles    Particle container (halo particles are moved out then back).
    @param[in]     a_patch        Per-patch deposit functor (requested kernel).
    @param[in]     a_patchNGP     Per-patch deposit functor (NGP).
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
        EBCellFAB&            rho    = (*a_meshData[lvl])[din];
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

    // Return the halo particles to the valid holder.
    particles.transferParticles(particles.getMaskParticles());
  }

  /**
    @brief Transition-mask width for a kernel: CIC -> refRat/2, TSC -> refRat (NGP is handled upstream).
    @param[in] a_depositionType Deposition kernel.
    @param[in] a_refRat         Refinement ratio between the level and the finer one.
    @return The transition band width (in fine cells).
  */
  int
  getTransitionMaskWidth(const DepositionType a_depositionType, const int a_refRat) const
  {
    int maskWidth = 0;
    if (a_depositionType == DepositionType::CIC) {
      maskWidth = a_refRat / 2;
    }
    else if (a_depositionType == DepositionType::TSC) {
      maskWidth = a_refRat;
    }
    else {
      MayDay::Abort("EBAMRParticleMesh::getTransitionMaskWidth - logic bust");
    }
    if (maskWidth > m_ghost) {
      MayDay::Abort("EBAMRParticleMesh::getTransitionMaskWidth - not enough ghost cells");
    }
    return maskWidth;
  }

  /**
    @brief Build the per-width transition masks: the fine-side transition band, on the refined-coarse grid.
    @details For each width in [1, m_ghost] and each level with a finer level, builds the fine-grid
    coarse-fine interface band (fine ghost cells just outside the fine boxes), transfers it to the
    refined-coarse (FiCo) grid of the finer level, and stores it as a boolean mask there. Storage-
    agnostic. Empty BaseFabs are undefined.
  */
  void
  defineTransitionMasks()
  {
    CH_TIME("EBAMRParticleMesh::defineTransitionMasks");

    constexpr int  comp    = 0;
    constexpr int  numComp = 1;
    const Interval interv  = Interval(comp, comp);

    m_transitionMasks.clear();

    for (int ighost = 1; ighost <= m_ghost; ighost++) {
      Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> mask(1 + m_finestLevel);

      for (int lvl = 0; lvl < m_finestLevel; lvl++) {
        const DisjointBoxLayout& gridsFine  = m_eblgs[lvl + 1]->getDBL();
        const EBLevelGrid&       eblgFiCo   = m_coarseFinePM[lvl + 1]->getEblgFiCo();
        const DisjointBoxLayout& gridsFiCo  = eblgFiCo.getDBL();
        const ProblemDomain&     domainFine = m_eblgs[lvl + 1]->getDomain();

        mask[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(
          new LevelData<BaseFab<bool>>(gridsFiCo, numComp, IntVect::Zero));
        LevelData<BaseFab<bool>>& levelMask = *mask[lvl];
        LevelData<FArrayBox>      cfivsFine(gridsFine, numComp, ighost * IntVect::Unit);
        LevelData<FArrayBox>      cfivsFiCo(gridsFiCo, numComp, IntVect::Zero);

        {
          const DataIterator& dit  = gridsFiCo.dataIterator();
          const int           nbox = dit.size();
#pragma omp parallel for schedule(runtime)
          for (int mybox = 0; mybox < nbox; mybox++) {
            const DataIndex& din = dit[mybox];
            levelMask[din].setVal(false);
            cfivsFiCo[din].setVal(0.0);
          }
        }

        // Fine-grid coarse-fine interface band: ghost cells outside the fine boxes (minus neighbours).
        {
          const DataIterator& dit  = gridsFine.dataIterator();
          const int           nbox = dit.size();
#pragma omp parallel for schedule(runtime)
          for (int mybox = 0; mybox < nbox; mybox++) {
            const DataIndex& din      = dit[mybox];
            const Box        cellBox  = gridsFine[din];
            const Box        ghostBox = grow(cellBox, ighost) & domainFine;

            cfivsFine[din].setVal(1.0, ghostBox, comp);
            cfivsFine[din].setVal(0.0, cellBox, comp);

            NeighborIterator nit(gridsFine);
            for (nit.begin(din); nit.ok(); ++nit) {
              cfivsFine[din].setVal(0.0, ghostBox & gridsFine[nit()], comp);
            }
          }
        }

        Copier copier;
        copier.ghostDefine(gridsFine, gridsFiCo, domainFine, ighost * IntVect::Unit);
        cfivsFine.copyTo(interv, cfivsFiCo, interv, copier, LDaddOp<FArrayBox>());

        const DataIterator& dit  = gridsFiCo.dataIterator();
        const int           nbox = dit.size();
#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex& din      = dit[mybox];
          const Box        cellBox  = gridsFiCo[din];
          BaseFab<bool>&   boolMask = levelMask[din];
          const FArrayBox& realMask = cfivsFiCo[din];

          bool emptyMask = true;
          for (BoxIterator bit(cellBox); bit.ok(); ++bit) {
            if (realMask(bit(), comp) > 0.5) {
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
      m_transitionMasks.emplace(ighost, mask);
    }
  }

  /**
    @brief Move the coarse-side transition particles into the container mask holder (using the FiCo mask).
    @details The transition mask lives on the refined-coarse grid, so this is done directly here
    (the container's transferMaskParticles works on the valid grids). Coarse particles whose
    fine-resolution cell is marked are moved (swap-and-pop) from the valid to the mask holder.
    @tparam Container ParticleContainer type (deduced).
    @param[in,out] a_particles      Particle container.
    @param[in]     a_depositionType Deposition kernel (selects the mask width; NGP is a no-op).
  */
  template <typename Container>
  void
  transferMaskParticlesTransition(Container& a_particles, const DepositionType a_depositionType) const
  {
    CH_TIME("EBAMRParticleMesh::transferMaskParticlesTransition");

    if (a_depositionType == DepositionType::NGP) {
      return;
    }

    for (int lvl = 0; lvl < m_finestLevel; lvl++) {
      const int                maskWidth = this->getTransitionMaskWidth(a_depositionType, m_refRat[lvl]);
      const DisjointBoxLayout& dbl       = m_eblgs[lvl]->getDBL();
      const Real               dxFine    = m_dx[lvl + 1];

      const DataIterator& dit  = dbl.dataIterator();
      const int           nbox = dit.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&     din  = dit[mybox];
        const BaseFab<bool>& mask = (*m_transitionMasks.at(maskWidth)[lvl])[din];
        if (!mask.isUsable()) {
          continue;
        }
        const Box maskBox = mask.box();

        auto& valid    = a_particles[lvl][din];
        auto& maskLeaf = (*a_particles.getMaskParticles()[lvl])[din];

        std::size_t i = 0;
        while (i < valid.size()) {
          const RealVect pos = valid.position(i);
          IntVect        iv;
          for (int dir = 0; dir < SpaceDim; dir++) {
            iv[dir] = static_cast<int>(std::floor((pos[dir] - m_probLo[dir]) / dxFine));
          }
          if (maskBox.contains(iv) && mask(iv, 0)) {
            maskLeaf.append(pos, valid.weight(i), valid.gather(i));
            const std::size_t k    = maskLeaf.size() - 1;
            maskLeaf.particleID(k) = valid.particleID(i);
            maskLeaf.rankID(k)     = valid.rankID(i);
            valid.remove(i);
          }
          else {
            i++;
          }
        }
      }
    }
  }

  /**
    @brief The CoarseFineDeposition::Transition deposition algorithm.
    @details NGP short-circuits to the Interp result. Otherwise: moves the coarse-side transition
    particles out of the valid holder; deposits the remaining particles with the kernel and folds
    fine-ghost mass to coarse; then deposits the transition particles on the refined-coarse grid with
    the fine-grid width and distributes that buffer to both the fine level (addFiCoDataToFine +
    exchangeAndAddFiCoData) and the coarse level (restrictAndAddFiCoDataToCoar); finally returns the
    transition particles to the valid holder.
    @tparam NCOMP        Number of mesh components (compile-time).
    @tparam Container    ParticleContainer type (deduced).
    @tparam PatchDeposit Per-patch deposit with the requested kernel.
    @param[in,out] a_meshData       AMR mesh data (RESET then filled).
    @param[in]     a_particles      Particle container.
    @param[in]     a_depositionType Deposition kernel.
    @param[in]     a_patch          Per-patch deposit functor (requested kernel).
  */
  template <int NCOMP, typename Container, typename PatchDeposit>
  void
  depositTransitionCore(EBAMRCellData&       a_meshData,
                        const Container&     a_particles,
                        const DepositionType a_depositionType,
                        PatchDeposit         a_patch) const
  {
    CH_TIME("EBAMRParticleMesh::depositTransitionCore");

    CH_assert(m_isDefined);
    CH_assert(a_meshData[0]->nComp() == NCOMP);

    if (a_depositionType == DepositionType::NGP) {
      this->depositInterpCore(a_meshData, a_particles, NCOMP, a_patch);
      return;
    }

    DataOps::setValue(a_meshData, 0.0);

    Container& particles = const_cast<Container&>(a_particles);
    this->transferMaskParticlesTransition(particles, a_depositionType);

    // 1. Deposit the non-transition particles on their levels.
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl  = m_eblgs[lvl]->getDBL();
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&      din    = dit[mybox];
        const EBParticleMesh& interp = (*m_ebParticleMesh[lvl])[din];
        EBCellFAB&            rho    = (*a_meshData[lvl])[din];
        a_patch(interp, rho, a_particles[lvl][din], 1.0);
      }
    }

    // 2. Fold ghost mass into the owner patch; fine-ghost mass over the boundary goes to coarse.
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      a_meshData[lvl]->exchange(Interval(0, NCOMP - 1), m_levelCopiers[lvl], EBAddOp());
      if (lvl > 0) {
        for (int comp = 0; comp < NCOMP; comp++) {
          LevelData<EBCellFAB> coarAlias;
          LevelData<EBCellFAB> fineAlias;
          aliasLevelData<EBCellFAB>(coarAlias, &(*a_meshData[lvl - 1]), Interval(comp, comp));
          aliasLevelData<EBCellFAB>(fineAlias, &(*a_meshData[lvl]), Interval(comp, comp));
          m_coarseFinePM[lvl]->addFineGhostsToCoarse(coarAlias, fineAlias);
        }
      }
    }

    // 3a. Deposit the transition particles on the refined-coarse grid (fine-grid width).
    for (int lvl = 0; lvl < m_finestLevel; lvl++) {
      const int                maskWidth  = this->getTransitionMaskWidth(a_depositionType, m_refRat[lvl]);
      const DisjointBoxLayout& dblFiCo    = m_coarseFinePM[lvl + 1]->getEblgFiCo().getDBL();
      LevelData<EBCellFAB>&    bufferFiCo = m_coarseFinePM[lvl + 1]->template getBufferFiCo<NCOMP>();

      const DataIterator& ditFiCo  = dblFiCo.dataIterator();
      const int           nboxFiCo = ditFiCo.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nboxFiCo; mybox++) {
        const DataIndex& din      = ditFiCo[mybox];
        EBCellFAB&       dataFiCo = bufferFiCo[din];
        dataFiCo.setVal(0.0);

        const BaseFab<bool>& mask = (*m_transitionMasks.at(maskWidth)[lvl])[din];
        if (mask.isUsable()) {
          const EBParticleMesh& interp = (*m_ebParticleMeshFiCo[lvl + 1])[din];
          a_patch(interp, dataFiCo, (*a_particles.getMaskParticles()[lvl])[din], 1.0);
        }
      }
    }

    // 3b/c/d. Distribute the refined-coarse buffer to the fine level and (restricted) to the coarse level.
    for (int lvl = 0; lvl < m_finestLevel; lvl++) {
      LevelData<EBCellFAB>& bufferFiCo = m_coarseFinePM[lvl + 1]->template getBufferFiCo<NCOMP>();
      m_coarseFinePM[lvl + 1]->addFiCoDataToFine(*a_meshData[lvl + 1], bufferFiCo);
      m_coarseFinePM[lvl + 1]->exchangeAndAddFiCoData(bufferFiCo);
      m_coarseFinePM[lvl + 1]->restrictAndAddFiCoDataToCoar(*a_meshData[lvl],
                                                            bufferFiCo,
                                                            EBCoarseFineParticleMesh::Average::Conservative);
    }

    // Return the transition particles to the valid holder.
    particles.transferParticles(particles.getMaskParticles());
  }
};

#include <CD_NamespaceFooter.H>

#endif
