/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_MirrorDeposition.H
 * @brief  Namespace containing the geometry used by mirrored cut-cell deposition.
 * @author Robert Marskar
 */

#ifndef CD_MIRRORDEPOSITION_H
#define CD_MIRRORDEPOSITION_H

// Chombo includes
#include <REAL.H>
#include <RealVect.H>
#include <SPACE.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
 * @brief Geometry for mirrored cut-cell deposition, i.e. the even extension of the density about an embedded boundary.
 * @details Deposition divides by dx^D and never by the volume fraction, so a cut cell holds kappa*n where the fluid
 * solvers hold n. The mirror removes that by depositing each particle's cloud AND the cloud of its image reflected
 * across the boundary. This namespace holds the per-particle arithmetic for that reflection -- the local quadratic
 * surface patch, the image position, and the Jacobian that reweights it -- with no mesh and no I/O, so that it can be
 * unit tested directly against the prototype harnesses in Prototypes/MirrorDeposition.
 * @details The surface patch is stored per band cell as (status, x_c, n_c, S_c), where S_c is the shape operator
 * lifted into WORLD coordinates. Storing it frame-free is deliberate: a 2x2 shape operator "in the tangent frame"
 * only means something together with the frame it was fitted in, and a mismatch between the fitting frame and the
 * using frame moves the image by up to 2.24 dx while leaving tr S, det S, the Jacobian and any mass check exactly
 * right -- an error invisible to every check the scheme has. See Prototypes/MirrorDeposition/mirror_frame.py.
 */
namespace MirrorDeposition {

/**
 * @brief Number of tangent directions on the surface. One in 2-D, two in 3-D.
 */
constexpr int numTangents = SpaceDim - 1;

/**
 * @brief Number of independent entries in a symmetric SpaceDim x SpaceDim tensor. Three in 2-D, six in 3-D.
 */
constexpr int numShapeComp = (SpaceDim * (SpaceDim + 1)) / 2;

/**
 * @brief Component holding the band cell's status. See MirrorDeposition::Status.
 */
constexpr int compStatus = 0;

/**
 * @brief First of the SpaceDim components holding the boundary centroid, as a PHYSICAL position.
 */
constexpr int compCentroid = 1;

/**
 * @brief First of the SpaceDim components holding the unit normal, pointing INTO the fluid.
 */
constexpr int compNormal = compCentroid + SpaceDim;

/**
 * @brief First of the numShapeComp components holding the world-frame shape operator's upper triangle.
 */
constexpr int compShape = compNormal + SpaceDim;

/**
 * @brief Total number of components in the surface-data holder. Eight in 2-D, thirteen in 3-D.
 */
constexpr int numComp = compShape + numShapeComp;

/**
 * @brief What a band cell's surface data is good for.
 * @details Kept as three states rather than a plain ok/not-ok flag because the two failure modes are different
 * problems with different responses. A cell that merely failed the curvature fit still has a valid centroid and
 * normal, so it can and must still reflect -- against a plane. A cell that found no cut cell at all cannot reflect
 * at any quality. Counting them separately is what makes the diagnostic meaningful.
 */
enum class Status
{
  /**
   * @brief No surface data. The band cell found no cut cell whose data was actually delivered to this level.
   * @details Not an error in itself -- it says the level's grids do not cover the embedded-boundary band, which is
   * a tagging property. Particles in such a cell must not be reflected.
   */
  None = 0,

  /**
   * @brief Full quadratic patch: centroid, normal, and a fitted shape operator.
   */
  Fitted = 1,

  /**
   * @brief Centroid and normal, but the curvature fit was refused, so the shape operator is zero.
   * @details A zero shape operator gives a Jacobian of exactly one through jacobian() below, i.e. the surface is
   * treated as locally flat. No branch is needed on the per-particle path; this is a diagnostic label only.
   */
  Planar = 2
};

/**
 * @brief Why reflect() refused to trust its Jacobian.
 */
enum class Refusal
{
  /**
   * @brief No refusal; the Jacobian was accepted.
   */
  None,

  /**
   * @brief The Jacobian's denominator came too close to zero.
   * @details 1 + 2*H*d + K*d^2 vanishes at d = -1/c_i, i.e. at the centre of curvature of a concave surface. The
   * image is then on top of the focus and its weight is unbounded.
   */
  SmallDenominator,

  /**
   * @brief The Jacobian exceeded the permitted magnitude.
   */
  LargeJacobian,

  /**
   * @brief The Jacobian came out non-positive, i.e. the reflection is past the surface's centre of curvature.
   */
  NonPositiveJacobian
};

/**
 * @brief Lift a shape operator from the tangent frame into world coordinates.
 * @details Computes S_c = P S P^T with P = [t_1 ... t_(D-1)], which is the frame-free form: the tangent basis is
 * consumed here and never stored, so no later reader has to agree with the fit about which basis was used.
 * @param[in]  a_tangents Orthonormal tangent basis used by the fit. Must be perpendicular to the normal.
 * @param[in]  a_S        Shape operator in that tangent basis, symmetric.
 * @param[out] a_shape    Upper triangle of the world-frame shape operator, numShapeComp entries.
 */
inline void
liftShapeOperator(const RealVect (&a_tangents)[numTangents],
                  const Real (&a_S)[numTangents][numTangents],
                  Real (&a_shape)[numShapeComp]) noexcept;

/**
 * @brief Index into the stored upper triangle of a symmetric SpaceDim x SpaceDim tensor.
 * @details The triangle is stored row-major: (0,0),(0,1),...,(0,D-1),(1,1),...  Indices may be given in either
 * order.
 * @param[in] a_i Row index.
 * @param[in] a_j Column index.
 * @return Offset of the (a_i, a_j) entry within the numShapeComp stored values.
 */
inline int
shapeIndex(const int a_i, const int a_j) noexcept;

/**
 * @brief Read one entry of a world-frame shape operator stored as an upper triangle.
 * @param[in] a_shape Upper triangle, numShapeComp entries.
 * @param[in] a_i     Row index.
 * @param[in] a_j     Column index.
 * @return The (a_i, a_j) entry. The tensor is symmetric, so the order of the indices does not matter.
 */
inline Real
shapeEntry(const Real* a_shape, const int a_i, const int a_j) noexcept;

/**
 * @brief Apply a world-frame shape operator to a vector.
 * @param[in] a_shape Upper triangle, numShapeComp entries.
 * @param[in] a_vec   Vector to apply it to.
 * @return S_c * a_vec.
 */
inline RealVect
applyShapeOperator(const Real* a_shape, const RealVect& a_vec) noexcept;

/**
 * @brief Twice the mean curvature, i.e. the first invariant of the shape operator.
 * @param[in] a_shape Upper triangle of the world-frame shape operator.
 * @return tr(S_c), which is c_1 + c_2 in 3-D and c_1 in 2-D.
 */
inline Real
meanCurvatureTimesTwo(const Real* a_shape) noexcept;

/**
 * @brief Gaussian curvature, i.e. the second invariant of the shape operator.
 * @details Computed as (1/2)[(tr S_c)^2 - tr(S_c^2)], NOT as det(S_c). The world-frame shape operator annihilates
 * the normal, so it has rank at most SpaceDim - 1 and its determinant is identically zero -- using det here would
 * silently set K = 0 and degrade the exact Jacobian to the linearized one, whose error runs to 48-77% at the radii
 * this scheme targets. The second-invariant form is correct in both dimensions and needs no 2-D special case: in
 * 2-D, S_c = c*t*t^T gives tr(S_c) = c and tr(S_c^2) = c^2, so K comes out zero on its own.
 * @param[in] a_shape Upper triangle of the world-frame shape operator.
 * @return c_1*c_2 in 3-D, zero in 2-D.
 */
inline Real
gaussianCurvature(const Real* a_shape) noexcept;

/**
 * @brief The exact area/volume Jacobian of the reflection across a curved surface.
 * @details J = (1 - 2*H*d + K*d^2)/(1 + 2*H*d + K*d^2) = prod_i (1 - c_i*d)/(1 + c_i*d) over the principal
 * curvatures. This is a product over principal curvatures, not a power (r'/r)^(D-1) -- the power form over-corrects
 * a torus by 44%.
 * @param[in]  a_twoH        Twice the mean curvature.
 * @param[in]  a_K           Gaussian curvature.
 * @param[in]  a_d           Signed distance to the surface patch.
 * @param[out] a_denominator The denominator, returned so the caller can guard on it directly.
 * @return The Jacobian. Positivity is conditional, not structural, so the caller must guard.
 */
inline Real
jacobian(const Real a_twoH, const Real a_K, const Real a_d, Real& a_denominator) noexcept;

/**
 * @brief Reflect a position across the quadratic surface patch stored for its cell, and weight the image.
 * @details Implements the first iterate of the signed distance to the patch,
 *
 *   w    = x - x_c,  eta = n_c . w,  d = eta + (1/2) w . (S_c w),
 *   nhat = normalize(n_c + S_c w),
 *   R(x) = x - 2*d*nhat,
 *
 * followed by the Jacobian above. Note that d is the FIRST ITERATE of the signed distance, not the signed distance
 * itself; that is what the error budget of the scheme was measured against.
 * @details On refusal the image is still produced, with a Jacobian of one. A refused image is deliberately NOT
 * dropped: dropping it removes the very correction the mirror exists to apply, in exactly the cells where the
 * correction is largest.
 * @param[in]  a_pos            Particle position, physical.
 * @param[in]  a_surfaceData    The cell's numComp surface components, gathered contiguously. See the note below.
 * @param[in]  a_maxJacobian    Refuse Jacobians above this magnitude.
 * @param[in]  a_minDenominator Refuse denominators below this magnitude.
 * @param[out] a_image          The reflected position.
 * @param[out] a_jacobian       The Jacobian, or one if refused.
 * @param[out] a_refusal        Which guard fired, or Refusal::None.
 * @return True if the Jacobian was accepted, false if a guard fired.
 * @note a_surfaceData must be a CONTIGUOUS array of numComp Reals, so a caller holding an EBCellFAB has to gather
 * the components rather than pass a pointer into it. A BaseFab stores the component index last -- all cells of
 * component 0, then all cells of component 1 -- so &fab(iv, 0) walks into the NEXT CELL, not the next component.
 */
inline bool
reflect(const RealVect& a_pos,
        const Real*     a_surfaceData,
        const Real      a_maxJacobian,
        const Real      a_minDenominator,
        RealVect&       a_image,
        Real&           a_jacobian,
        Refusal&        a_refusal) noexcept;

} // namespace MirrorDeposition

#include <CD_NamespaceFooter.H>

#include <CD_MirrorDepositionImplem.H>

#endif
