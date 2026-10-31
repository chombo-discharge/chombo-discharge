/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBStateRedistribution.H
 * @brief  Declaration of weighted state redistribution as an EBRedistribution.
 * @author Robert Marskar
 */

#ifndef CD_EBSTATEREDISTRIBUTION_H
#define CD_EBSTATEREDISTRIBUTION_H

// Std includes
#include <vector>

// Chombo includes
#include <EBLevelGrid.H>
#include <ProblemDomain.H>
#include <RealVect.H>
#include <LevelData.H>
#include <FArrayBox.H>
#include <BaseFab.H>
#include <LayoutData.H>
#include <RefCountedPtr.H>
#include <BaseIVFAB.H>
#include <Copier.H>

// Our includes
#include <CD_EBRedistribution.H>
#include <CD_EBFluxRedistribution.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Weighted state redistribution (Berger and Giuliani, J. Comput. Phys. 428, 2021; Giuliani et al., J. Comput.
 * Phys. 464, 2022), as an EBRedistribution.
 * @details The algorithm follows the implementation in AMReX (AMReX_EB_StateRedistribute.cpp and the files it uses).
 * A cut cell i with volume fraction kappa_i below the target volume fraction is merged with up to 2^D - 1 of its
 * neighbours: first the neighbour in the direction of the largest component of the EB normal, then, while the merged
 * volume is below the target or the normal components tie, the neighbours that complete a 2x2 (and in 3D a 2x2x2)
 * block. Each cell r is in N_r neighbourhoods, its own included. The neighbourhood of i has the weights
 *
 *    alpha_1(i) = (target - kappa_i)/sum(kappa_r)  on the merged neighbours r,
 *    alpha_0(r) = 1 - sum(alpha_1(i)/N_r)           on each cell r, over the neighbourhoods that contain r,
 *
 * the volume V_i = alpha_0(i)*kappa_i + sum(alpha_1(i)*kappa_r/N_r) and the centroid c_i (the same weighted average of
 * the cell centroids). An update phi^{n+1} = phi^n - dt*D starts from the provisional state U = phi^n - dt*D, forms the
 * neighbourhood averages Q_i with the same weights, a limited least-squares slope of Q over the neighbouring
 * neighbourhoods, and gives each cell
 *
 *    phi_r^{n+1} = (alpha_0(r)*N_r*q_r(x_r) + sum(alpha_1(i)*q_i(x_r)))/N_r,
 *
 * where q_i is the linear reconstruction of neighbourhood i and the sum runs over the merged neighbourhoods that
 * contain r. The right-hand side returned by redistribute() is R = (phi^n - phi^{n+1})/dt, which depends on phi^n and
 * dt. A cell in no merged neighbourhood keeps R = D.
 *
 * The neighbourhoods reach one cell, and the slopes up to two cells, so the state and the divergence must have two
 * ghost cells, filled on refinement boundaries. Cells near small cells must be single-valued.
 *
 * On AMR hierarchies the redistribution follows Barrio Sanchez et al. (arXiv:2309.06372, 2023) and AMReX
 * (ApplyMLRedistribution). Each level is redistributed separately. The neighbourhoods are centred on the valid cells of
 * the level, those inside its grids and not covered by the finer level, but they can also contain cells that stand in
 * for another level: cells covered by the finer level, and ghost cells outside the grids of the level, which overlay
 * valid cells of the coarser level. The state of a ghost cell is that of the coarse cell under it. The stand-ins also
 * enter the slopes and the limiter. The redistribution conserves mass in each neighbourhood. Its update is split into
 * transfers between pairs of cells, the part from the slopes included, and the transfers between valid cells and
 * stand-ins are moved to the cells that the stand-ins represent:
 *
 *   - a transfer to a ghost cell goes to the coarse cell under it. As in the paper, a coarse cell with volume fraction
 *     kappa takes the fraction kappa of it, and the rest is redistributed to its neighbours with the
 *     EBFluxRedistribution operators.
 *   - a transfer to a covered cell goes to the fine cells under it, as the same change of the state in each of them.
 *     The paper instead gives it back to the valid coarse cell, which undoes the stabilization of a small coarse cell
 *     whose neighbours are all covered.
 *
 * The update is then conservative over the composite grid when the fluxes at refinement boundaries are consistent
 * between the levels.
 */
class EBStateRedistribution : public EBRedistribution
{
public:
  /**
   * @brief Disallowed weak construction.
   */
  EBStateRedistribution() = delete;

  /**
   * @brief Number of ghost cells of the neighbourhood data and of the cell masks. Also the number of ghost cells needed
   * in the state and the divergence.
   */
  static constexpr int s_numGhost = 2;

  /**
   * @brief Constructor. Builds the neighbourhoods on all levels.
   * @param[in] a_eblg             Grids on each level
   * @param[in] a_eblgCoFi         Coarsening of the grids of each level: entry l is the coarsening of the grids on
   * level l + 1.
   * @param[in] a_refinementRatios Refinement ratios between levels
   * @param[in] a_finestLevel      Finest grid level
   * @param[in] a_order            Order of the reconstruction in the neighbourhoods: 0 (constant) or 2 (linear,
   * limited)
   * @param[in] a_targetVolFrac    Target volume fraction. Cut cells with a smaller volume fraction are merged.
   * @param[in] a_fluxRedistribution Flux redistribution operators on each level, for the corrections on refinement
   * boundaries. Only used if there is more than one level.
   * @param[in] a_coveredByFiner   Cells covered by the finer level, on each level, with s_numGhost ghost cells.
   * @param[in] a_coarseFineGhost  Cells inside the domain but outside the grids of each level, with s_numGhost ghost
   * cells. They overlay valid cells of the coarser level.
   */
  EBStateRedistribution(const Vector<RefCountedPtr<EBLevelGrid>>&              a_eblg,
                        const Vector<RefCountedPtr<EBLevelGrid>>&              a_eblgCoFi,
                        const Vector<int>&                                     a_refinementRatios,
                        const int                                              a_finestLevel,
                        const int                                              a_order,
                        const Real                                             a_targetVolFrac,
                        const Vector<RefCountedPtr<EBFluxRedistribution>>&     a_fluxRedistribution,
                        const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_coveredByFiner,
                        const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_coarseFineGhost) noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~EBStateRedistribution() noexcept;

  /**
   * @brief Turn kappa*D into the right-hand side of the update stabilized by state redistribution.
   * @details See EBRedistribution::redistribute. a_dt must not be zero.
   * @param[in,out] a_divG On input kappa*D, on output R.
   * @param[in]     a_phi  State that the update starts from.
   * @param[in]     a_dt   Time step of the update, signed.
   */
  void
  redistribute(EBAMRCellData& a_divG, const EBAMRCellData& a_phi, const Real a_dt) noexcept override;

protected:
  /**
   * @brief Number of cells in the 3^D block around a cell, which indexes the neighbour offsets.
   */
  static constexpr int s_numOffsets = D_TERM(3, *3, *3);

  /**
   * @brief Number of components of the tracker: the number of merged neighbours, and up to 2^D - 1 neighbours.
   */
  static constexpr int s_trackerComps = D_TERM(2, *2, *2);

  /**
   * @brief Cell type: outside the domain, or not reachable.
   */
  static constexpr int s_outside = 0;

  /**
   * @brief Cell type: valid cell of the level.
   */
  static constexpr int s_valid = 1;

  /**
   * @brief Cell type: covered by the finer level.
   */
  static constexpr int s_covered = 2;

  /**
   * @brief Cell type: inside the domain but outside the grids of the level, over a valid cell of the coarser level.
   */
  static constexpr int s_ghost = 3;

  /**
   * @brief Finest grid level.
   */
  int m_finestLevel;

  /**
   * @brief Order of the reconstruction: 0 or 2.
   */
  int m_order;

  /**
   * @brief Target volume fraction.
   */
  Real m_targetVolFrac;

  /**
   * @brief Grids on each level.
   */
  Vector<RefCountedPtr<EBLevelGrid>> m_eblg;

  /**
   * @brief Coarsening of the grids: entry l is the coarsening of the grids on level l + 1.
   */
  Vector<RefCountedPtr<EBLevelGrid>> m_eblgCoFi;

  /**
   * @brief Cells covered by the finer level, with s_numGhost ghost cells.
   */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_coveredByFiner;

  /**
   * @brief Cells inside the domain but outside the grids of the level, with s_numGhost ghost cells.
   */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_coarseFineGhost;

  /**
   * @brief Refinement ratios between levels.
   */
  Vector<int> m_refinementRatios;

  /**
   * @brief Flux redistribution operators on each level, which redistribute the corrections on refinement boundaries.
   */
  Vector<RefCountedPtr<EBFluxRedistribution>> m_fluxRedistribution;

  /**
   * @brief Neighbourhoods. Component 0 is the number of merged neighbours of the cell, components 1, 2, ... are their
   * offsets, as indices into the 3^D block around the cell (see offset()).
   */
  Vector<RefCountedPtr<LevelData<BaseFab<int>>>> m_tracker;

  /**
   * @brief Number of neighbourhoods that each cell is in, its own included.
   */
  Vector<RefCountedPtr<LevelData<BaseFab<int>>>> m_numNbhd;

  /**
   * @brief Weights. Component 0 is alpha_0 of the cell, component 1 is alpha_1 of the neighbourhood of the cell.
   */
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_alpha;

  /**
   * @brief Volume of the neighbourhood of each cell.
   */
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_nbhdVol;

  /**
   * @brief Centroid of the neighbourhood of each cell, relative to the cell centre in units of the cell size.
   */
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_nbhdCentroid;

  /**
   * @brief Neighbours in the 3^D block around each merged cell that are connected to it, as a bit mask over the
   * offsets. Used in the slopes.
   */
  Vector<RefCountedPtr<LevelData<BaseFab<int>>>> m_connected;

  /**
   * @brief Per box, true if the box has a cell whose state is changed by the redistribution.
   */
  Vector<RefCountedPtr<LayoutData<bool>>> m_activeBox;

  /**
   * @brief Scratch: neighbourhood averages Q.
   */
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_nbhdAverage;

  /**
   * @brief Scratch: limited slopes of the merged neighbourhoods.
   */
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_slopes;

  /**
   * @brief Scratch: corrections of the valid cells, as mass in units of the cell volume of the level.
   */
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_syncMass;

  /**
   * @brief Scratch: the provisional state U = phi - dt*D. In the ghost cells it is U of the coarse cell under the ghost
   * cell.
   */
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_state;

  /**
   * @brief Scratch on the coarsening of the grids of a level, with one ghost cell, for data of the coarser level: U of
   * the coarse cells, the corrections from the ghost cells of the level (as mass in units of the cell volume of the
   * coarser level), and the corrections of the coarse cells covered by the level.
   */
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_ghostMass;

  /**
   * @brief Copiers that add m_ghostMass of a level, ghost cells included, to m_syncMass of the coarser level.
   */
  Vector<Copier> m_ghostCopier;

  /**
   * @brief Copiers from the coarser level to m_ghostMass of a level, ghost cells included.
   */
  Vector<Copier> m_coarseCopier;

  /**
   * @brief Scratch: the parts of the corrections that are redistributed with the flux redistribution, as rates.
   */
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real>>>> m_massDiff;

  /**
   * @brief Get the cell offset that an index into the 3^D block stands for.
   * @param[in] a_index Index into the 3^D block, 0 <= a_index < 3^D.
   * @return The offset, with components in {-1, 0, 1}.
   */
  static IntVect
  offset(const int a_index) noexcept;

  /**
   * @brief Get the index into the 3^D block of a cell offset.
   * @param[in] a_offset Offset with components in {-1, 0, 1}.
   * @return The index.
   */
  static int
  offsetIndex(const IntVect& a_offset) noexcept;

  /**
   * @brief Get the volume fraction of a cell from the EB geometry.
   * @param[in] a_ebisbox EB geometry of the grid box, defined at least s_numGhost cells around it.
   * @param[in] a_domain  Domain of the level.
   * @param[in] a_iv      Cell
   * @return The volume fraction; 0 for cells covered by the EB and outside the domain, -1 for multi-valued cells.
   */
  static Real
  volumeFraction(const EBISBox& a_ebisbox, const ProblemDomain& a_domain, const IntVect& a_iv) noexcept;

  /**
   * @brief Get the centroid of a cell from the EB geometry.
   * @param[in] a_ebisbox EB geometry of the grid box, defined at least s_numGhost cells around it.
   * @param[in] a_domain  Domain of the level.
   * @param[in] a_iv      Cell
   * @return The centroid as an offset from the cell centre in units of the cell size; zero except in single-valued cut
   * cells.
   */
  static RealVect
  cellCentroid(const EBISBox& a_ebisbox, const ProblemDomain& a_domain, const IntVect& a_iv) noexcept;

  /**
   * @brief Get the type of a cell.
   * @param[in] a_coveredByFiner  Cells covered by the finer level, around the grid box.
   * @param[in] a_coarseFineGhost Cells inside the domain but outside the grids of the level, around the grid box.
   * @param[in] a_domain          Domain of the level.
   * @param[in] a_iv              Cell
   * @return s_ghost, s_outside, s_covered or s_valid.
   */
  static int
  cellType(const BaseFab<bool>& a_coveredByFiner,
           const BaseFab<bool>& a_coarseFineGhost,
           const ProblemDomain& a_domain,
           const IntVect&       a_iv) noexcept;

  /**
   * @brief Build the neighbourhood data on one level.
   * @param[in] a_level           Grid level
   */
  void
  defineLevel(const int a_level) noexcept;

  /**
   * @brief Find the neighbours that a small cut cell merges with.
   * @param[out] a_offsets   Offsets of the merged neighbours.
   * @param[in]  a_iv        Cell
   * @param[in]  a_normal    EB normal of the cell, pointing into the fluid.
   * @param[in]  a_ebisbox   EB geometry of the grid box.
   * @param[in]  a_coveredByFiner  Cells covered by the finer level, around the grid box.
   * @param[in]  a_coarseFineGhost Cells inside the domain but outside the grids of the level, around the grid box.
   * @param[in]  a_domain    Domain of the level.
   */
  void
  mergeNeighbours(std::vector<IntVect>& a_offsets,
                  const IntVect&        a_iv,
                  const RealVect&       a_normal,
                  const EBISBox&        a_ebisbox,
                  const BaseFab<bool>&  a_coveredByFiner,
                  const BaseFab<bool>&  a_coarseFineGhost,
                  const ProblemDomain&  a_domain) const noexcept;

  /**
   * @brief Compute the limited slope of the neighbourhood average at a merged cell.
   * @param[in] a_iv          Merged cell
   * @param[in] a_average     Neighbourhood averages
   * @param[in] a_centroid    Neighbourhood centroids
   * @param[in] a_ebisbox     EB geometry of the grid box.
   * @param[in] a_coveredByFiner  Cells covered by the finer level, around the grid box.
   * @param[in] a_coarseFineGhost Cells inside the domain but outside the grids of the level, around the grid box.
   * @param[in] a_domain      Domain of the level.
   * @param[in] a_connected   Connected neighbours of the merged cell (bit mask over the 3^D block)
   * @return The slope, in units of the state per cell size.
   */
  RealVect
  computeSlope(const IntVect&       a_iv,
               const FArrayBox&     a_average,
               const FArrayBox&     a_centroid,
               const EBISBox&       a_ebisbox,
               const BaseFab<bool>& a_coveredByFiner,
               const BaseFab<bool>& a_coarseFineGhost,
               const ProblemDomain& a_domain,
               const int            a_connected) const noexcept;

  /**
   * @brief Add the corrections on refinement boundaries to the right-hand side.
   * @details On input, m_syncMass holds the corrections of the valid cells and m_ghostMass those from the ghost cells.
   * @param[in,out] a_divG Right-hand side R on all levels.
   * @param[in]     a_dt   Time step of the update.
   */
  void
  synchronize(EBAMRCellData& a_divG, const Real a_dt) noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
