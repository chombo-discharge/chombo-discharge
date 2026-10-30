/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CdrCTU.H
 * @brief  Declaration of a class which implements CdrMultigrid using MUSCL for advection.
 * @author Robert Marskar
 */

#ifndef CD_CDRCTU_H
#define CD_CDRCTU_H

// Our includes
#include <CD_CdrMultigrid.H>
#include <CD_CdrGodunov.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Class that uses a slope limited method for advection, in combination with corner transport upwind (CTU).
 */
class CdrCTU : public CdrMultigrid
{
public:
  /**
   * @brief Default constructor.
   */
  CdrCTU();

  /**
   * @brief Destructor (does nothing)
   */
  virtual ~CdrCTU();

  /**
   * @brief Parse class options to put object in usable state.
   */
  virtual void
  parseOptions() override;

  /**
   * @brief Parse runtime options
   */
  virtual void
  parseRuntimeOptions() override;

  /**
   * @brief MUSCL advection to faces
   * @param[out] a_facePhi  Phi on face centers
   * @param[in]  a_cellPhi  Phi on cell centers
   * @param[in]  a_extrapDt Time step (i.e. extrapolation) of the face-centered states.
   * @param[in]  a_dt       Time step of the update that the face states enter, or 0 if they enter no update. Used for
   * the EB inflow state on outgoing faces of cut cells and for the flux correction.
   * @details The flux correction keeps the update phi - a_dt*div(F) non-negative in single-valued cells away from
   * refinement boundaries if the first-order update is, which needs a_extrapDt = a_dt (CTU) or a_extrapDt = 0 (donor
   * cell). It counts the EB flux and the computed fluxes on domain faces, whether or not the caller adds them to the
   * divergence, so it can limit more than needed but not less.
   */
  virtual void
  advectToFaces(EBAMRFluxData&       a_facePhi,
                const EBAMRCellData& a_cellPhi,
                const Real           a_extrapDt,
                const Real           a_dt) override;

  /**
   * @brief Compute the largest possible advective time step (for explicit methods).
   * @details This computes dt = dx/max(|vx|,|vy|,|vz|), minimized over all grid levels and patches.
   * @return Maximum stable advection time step.
   */
  virtual Real
  computeAdvectionDt() override;

  /**
   * @brief Allocate internal storage and build pre-built VoF iterators.
   */
  virtual void
  allocate() override;

  /**
   * @brief Regrid method — also rebuilds the grown VoF iterators.
   * @param[in] a_lbase          Coarsest level that changes
   * @param[in] a_oldFinestLevel Finest level before regrid
   * @param[in] a_newFinestLevel Finest level after regrid
   */
  virtual void
  regrid(const int a_lbase, const int a_oldFinestLevel, const int a_newFinestLevel) override;

  /**
   * @brief Pre-regrid hook — releases the grown VoF iterators.
   * @param[in] a_lbase          Coarsest level that changes
   * @param[in] a_oldFinestLevel Finest level before regrid
   */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel) override;

protected:
  /**
   * @brief Pre-built VoFIterators over grow(validBox,1) for the irregular kernel in computeNormalSlopes.
   */
  Vector<RefCountedPtr<LayoutData<VoFIterator>>> m_grownVofIter;

  /**
   * @brief 1 in the cells of each level's grids (with one ghost cell, filled where other boxes cover it), 0 elsewhere.
   * The EB inflow state is only used where the EB flux is known, which excludes ghost cells on refinement boundaries.
   */
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_levelCells;

  /**
   * @brief Get the EB normal of a cut cell, pointing into the fluid.
   * @details This is EBISBox::normal, which points into the fluid.
   * @param[in] a_vof     Cut cell
   * @param[in] a_ebisbox EBISBox of the grid box
   * @return The unit normal.
   */
  RealVect
  getEbNormalIntoFluid(const VolIndex& a_vof, const EBISBox& a_ebisbox) const noexcept;

  /**
   * @brief Get the state of the fluid that enters a cut cell through the EB.
   * @details Where the flow leaves the EB (v.n > 0, with n the EB normal pointing into the fluid), the state is
   * phi_EB = max(0, -F_EB/(v.n)), with F_EB the EB flux (positive out of the fluid) and v the cell velocity. An EB flux
   * out of the fluid where the flow leaves the EB gives a zero inflow state.
   * @param[out] a_phiEB   The EB inflow state, if the flow leaves the EB.
   * @param[in]  a_vof     Cut cell
   * @param[in]  a_ebisbox EBISBox of the grid box
   * @param[in]  a_level   Grid level
   * @param[in]  a_dit     Grid index
   * @return True if the cell has an EB and the flow leaves it, false otherwise.
   */
  bool
  getEbInflowState(Real&            a_phiEB,
                   const VolIndex&  a_vof,
                   const EBISBox&   a_ebisbox,
                   const int        a_level,
                   const DataIndex& a_dit) const noexcept;

  /**
   * @brief Flux-corrected transport that keeps the update non-negative.
   * @details The first-order states (CTU without slopes) give a non-negative update for Courant numbers up to 1 in
   * regular cells. The difference between the second-order and the first-order states is a correction that can make
   * the update negative. For each cell, the corrections that take mass out of it are summed, and limited to the mass
   * that the first-order update leaves in the cell: the fraction R = min(1, mass/removed). The correction on a face is
   * multiplied by R of the cell that it takes mass from. The fractions are exchanged so that faces on box boundaries
   * are corrected the same way on both sides, and each face state stays single-valued, so the update is conservative.
   * Ghost cells on refinement boundaries and faces of multi-valued cells are not limited.
   * @param[in,out] a_facePhi    On input the second-order face states, on output the corrected states.
   * @param[in]     a_facePhiLow First-order face states.
   * @param[in]     a_cellPhi    Cell states.
   * @param[in]     a_dt         Time step of the update that the face states enter.
   */
  void
  fluxCorrection(EBAMRFluxData&       a_facePhi,
                 const EBAMRFluxData& a_facePhiLow,
                 const EBAMRCellData& a_cellPhi,
                 const Real           a_dt) noexcept;

  /**
   * @brief Correct the states on the outgoing faces of cut cells for the fluid that enters through the EB.
   * @details Where the flow leaves the EB, the fluid that crosses an outgoing face of a cut cell during the time step
   * comes from the cell only until the flow has carried the EB (as a plane through the EB centroid) to the face. After
   * that it is fluid that entered through the EB, whose state is phi_EB = -F_EB/(v.n), with n the EB normal pointing
   * into the fluid and F_EB the EB flux (positive out of the fluid). With t_EB the time it takes the flow to reach the
   * face from the EB plane, the face state becomes the time average over the step,
   *
   *    f * phi_face + (1-f) * phi_EB,    f = min(1, t_EB/dt).
   *
   * The correction only changes the face state, which is the same for the cells on both sides of the face, so the
   * update stays conservative.
   * @param[in,out] a_facePhi Face states in direction a_dir.
   * @param[in]     a_cellVel Cell-centered velocity.
   * @param[in]     a_faceVel Face-centered velocity in direction a_dir.
   * @param[in]     a_cellBox Grid box.
   * @param[in]     a_dir     Face direction.
   * @param[in]     a_level   Grid level.
   * @param[in]     a_dit     Grid index.
   * @param[in]     a_dt      Time step of the update that the face states enter.
   */
  void
  ebInflowState(EBFaceFAB&       a_facePhi,
                const EBCellFAB& a_cellVel,
                const EBFaceFAB& a_faceVel,
                const Box&       a_cellBox,
                const int        a_dir,
                const int        a_level,
                const DataIndex& a_dit,
                const Real       a_dt) noexcept;

  /**
   * @brief Build m_grownVofIter and m_levelCells for all AMR levels.
   */
  void
  defineIterators();

  /**
   * @brief Supported limiters
   */
  enum class Limiter
  {
    None,
    MinMod,
    Superbee,
    MonotonizedCentral,
  };

  /**
   * @brief If true, slopes are limited in the normal direction.
   */
  Limiter m_limiter;

  /**
   * @brief If true, the discretization uses CTU. Otherwise it's DCU (donor cell upwind)
   */
  bool m_useCTU;

  /**
   * @brief If true, the EB is an inflow boundary where the flow leaves it. The slope of a cut cell without a face on
   * one side uses the EB inflow state as the neighbour on that side, and the state on an outgoing face of a cut cell
   * accounts for the fluid that enters through the EB during the time step of the update. See CdrCTU::ebInflowState.
   */
  bool m_ebInflowState;

  /**
   * @brief If true, the second-order part of the face states is limited by flux-corrected transport so that the update
   * stays non-negative. See fluxCorrection.
   */
  bool m_fluxCorrection;

  /**
   * @brief Compute slopes.
   * @param[out] a_normalSlopes Cell slopes. Must have SpaceDim components (for slopes in each direction).
   * @param[in]  a_cellPhi      Cell-centered states
   * @param[in]  a_cellBox      Cell-centered box identifying the cells where we compute the slopes.
   * @param[in]  a_domain       Domain on grid level
   * @param[in]  a_level        AMR level
   * @param[in]  a_dit          Grid index
   */
  virtual void
  computeNormalSlopes(EBCellFAB&           a_normalSlopes,
                      const EBCellFAB&     a_cellPhi,
                      const Box&           a_cellBox,
                      const ProblemDomain& a_domain,
                      const int            a_level,
                      const DataIndex&     a_dit);

  /**
   * @brief Upwind/Riemann solve. This extrapolates the cell-centered data to the inside face centers and selects the
   * upwind side.
   * @details This uses the Colella CTU method where we extrapolate with the corner cells as well. a_normalSlopes holds
   * the slopes normal to the face direction, and the transverse slopes are added inside this method. This routine also
   * solves the corresponding Riemann problem (which is just the upwind face).
   * @param[out] a_facePhi      Face-centered states.
   * @param[in]  a_normalSlopes Normal cell slopes
   * @param[in]  a_cellPhi      Cell-centered states
   * @param[in]  a_cellVel      Cell-centered velocities.
   * @param[in]  a_faceVel      Face-centered velocities.
   * @param[in]  a_domain       Domain on grid level.
   * @param[in]  a_cellBox      Grid box/patch
   * @param[in]  a_level        AMR level
   * @param[in]  a_dit          Grid index.
   * @param[in]  a_dt           Time step.
   */
  virtual void
  upwind(EBFluxFAB&           a_facePhi,
         const EBCellFAB&     a_normalSlopes,
         const EBCellFAB&     a_cellPhi,
         const EBCellFAB&     a_cellVel,
         const EBFluxFAB&     a_faceVel,
         const ProblemDomain& a_domain,
         const Box&           a_cellBox,
         const int&           a_level,
         const DataIndex&     a_dit,
         const Real&          a_dt);

  /**
   * @brief Parse slope limiting on/off
   */
  virtual void
  parseSlopeLimiter();

  /**
   * @brief minmod slope function.
   * @param[in] dwl Left slope.
   * @param[in] dwr Right slope.
   * @return Limited slope.
   */
  static Real
  minmod(const Real& dwl, const Real& dwr) noexcept;

  /**
   * @brief Superbee slope limiter.
   * @param[in] dwl Left slope.
   * @param[in] dwr Right slope.
   * @return Limited slope.
   */
  static Real
  superbee(const Real& dwl, const Real& dwr) noexcept;

  /**
   * @brief Monotonized central difference slope limiter.
   * @param[in] dwl Left slope.
   * @param[in] dwr Right slope.
   * @return Limited slope.
   */
  static Real
  monotonizedCentral(const Real& dwl, const Real& dwr) noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
