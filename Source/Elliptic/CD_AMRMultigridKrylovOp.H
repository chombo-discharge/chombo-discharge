/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_AMRMultigridKrylovOp.H
  @brief  Linear-operator adapter that exposes an AMR-multigrid V-cycle as a Krylov preconditioner
  @author Robert Marskar
*/

#ifndef CD_AMRMultigridKrylovOp_H
#define CD_AMRMultigridKrylovOp_H

// Std includes
#include <cmath>

// Chombo includes
#include <LinearSolver.H>
#include <AMRMultiGrid.H>
#include <LevelData.H>
#include <BaseFab.H>
#include <RefCountedPtr.H>
#include <Vector.H>
#include <EBCellFAB.H>
#include <MFCellFAB.H>

// Our includes
#include <CD_EBHelmholtzOp.H>
#include <CD_MFHelmholtzOp.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Trait mapping a FAB type to the concrete chombo-discharge Helmholtz operator that owns it.
  @details Used by AMRMultigridKrylovOp to reach the operator's dotProductMasked method (not part of the Chombo
           AMRLevelOp interface). Specialized for EBCellFAB and MFCellFAB.
*/
template <class T>
struct AMRMultigridKrylovOpTraits;

/**
  @brief EBCellFAB specialization -> EBHelmholtzOp.
*/
template <>
struct AMRMultigridKrylovOpTraits<EBCellFAB>
{
  /**
    @brief Concrete operator type for EBCellFAB.
  */
  using OpType = EBHelmholtzOp;
};

/**
  @brief MFCellFAB specialization -> MFHelmholtzOp.
*/
template <>
struct AMRMultigridKrylovOpTraits<MFCellFAB>
{
  /**
    @brief Concrete operator type for MFCellFAB.
  */
  using OpType = MFHelmholtzOp;
};

/**
  @brief Adapter presenting an already-defined AMRMultiGrid as a LinearOp over the AMR hierarchy.
  @details This lets an outer Krylov solver (GMRES/BiCGStab, operating on Vector<LevelData<T>*>) use the
           multigrid V-cycle as a preconditioner. The matrix-vector product and the preconditioner delegate to
           the supplied, fully set-up AMRMultiGrid (its AMR operators and AMRVCycle). The vector operations
           delegate to the per-level AMR operators, and the AMR inner product is computed allocation-free using
           a per-level fine-covered mask precomputed in define().
  @note T is the FAB type (EBCellFAB or MFCellFAB); the LinearOp acts on Vector<LevelData<T>*>.
*/
template <class T>
class AMRMultigridKrylovOp : public LinearOp<Vector<LevelData<T>*>>
{
public:
  /**
    @brief Multilevel data type the LinearOp acts on.
  */
  using Vec = Vector<LevelData<T>*>;

  /**
    @brief Default constructor. Call define() before use.
  */
  AMRMultigridKrylovOp() = default;

  /**
    @brief Destructor. Frees preconditioner scratch storage.
  */
  virtual ~AMRMultigridKrylovOp();

  /**
    @brief Define the adapter. Call once per regrid (after the AMRMultiGrid has been defined and init'd).
    @param[in] a_mg         The fully set-up multigrid solver to reuse for matvec and preconditioning
    @param[in] a_grids      Grids on each level (index 0..a_lmax)
    @param[in] a_refRat     Refinement ratio to the next finer level (index 0..a_lmax)
    @param[in] a_dx         Grid spacing on each level (index 0..a_lmax)
    @param[in] a_lbase      Coarsest level in the solve
    @param[in] a_lmax       Finest level in the solve
    @param[in] a_numVCycles Number of V-cycles applied per preconditioner invocation
  */
  void
  define(AMRMultiGrid<LevelData<T>>*      a_mg,
         const Vector<DisjointBoxLayout>& a_grids,
         const Vector<int>&               a_refRat,
         const Vector<Real>&              a_dx,
         const int                        a_lbase,
         const int                        a_lmax,
         const int                        a_numVCycles);

  /**
    @brief Matrix-vector product L(phi) over the AMR hierarchy (always homogeneous for the Krylov path).
    @param[out] a_lhs        Result L(phi)
    @param[in]  a_phi        Input vector
    @param[in]  a_homogeneous Use homogeneous boundary conditions
  */
  void
  applyOp(Vec& a_lhs, const Vec& a_phi, bool a_homogeneous = false) override;

  /**
    @brief AMR residual a_lhs = L(a_phi) - a_rhs.
    @param[out] a_lhs        Residual
    @param[in]  a_phi        Solution estimate
    @param[in]  a_rhs        Right-hand side
    @param[in]  a_homogeneous Use homogeneous boundary conditions
  */
  void
  residual(Vec& a_lhs, const Vec& a_phi, const Vec& a_rhs, bool a_homogeneous = false) override;

  /**
    @brief Apply the preconditioner: a_cor approx L^{-1} a_residual via m_numVCycles multigrid V-cycles.
    @param[out] a_cor      Preconditioned correction
    @param[in]  a_residual Residual to precondition
  */
  void
  preCond(Vec& a_cor, const Vec& a_residual) override;

  /**
    @brief Allocate a_lhs to mirror a_rhs.
    @param[out] a_lhs New holder
    @param[in]  a_rhs Template holder
  */
  void
  create(Vec& a_lhs, const Vec& a_rhs) override;

  /**
    @brief Free a holder created by create().
    @param[in,out] a_lhs Holder to free
  */
  void
  clear(Vec& a_lhs) override;

  /**
    @brief Copy a_rhs into a_lhs.
    @param[out] a_lhs Destination
    @param[in]  a_rhs Source
  */
  void
  assign(Vec& a_lhs, const Vec& a_rhs) override;

  /**
    @brief AMR-composite inner product (allocation-free, fine-covered cells masked out).
    @param[in] a_1 First operand
    @param[in] a_2 Second operand
    @return Inner product
  */
  Real
  dotProduct(const Vec& a_1, const Vec& a_2) override;

  /**
    @brief Batched inner products a_mdots[k] = dotProduct(a_1, a_2[k]) for k=0..a_sz-1.
    @details Overrides the LinearOp default (which loops dotProduct, one MPI reduction each). Here all a_sz*nLevel
             rank-local partials are reduced in a single MPI_Allreduce. This is the per-iteration cost of Krylov
             orthogonalization (Gram-Schmidt), so the single fused reduction is the dominant scaling win.
    @param[in]  a_1     Fixed operand
    @param[in]  a_sz    Number of right operands
    @param[in]  a_2     Array of a_sz right operands
    @param[out] a_mdots Array of a_sz inner products
  */
  void
  mDotProduct(const Vec& a_1, const int a_sz, const Vec a_2[], Real a_mdots[]) override;

  /**
    @brief Increment a_lhs += a_scale*a_x.
    @param[in,out] a_lhs   Accumulator
    @param[in]     a_x     Increment
    @param[in]     a_scale Scale factor
  */
  void
  incr(Vec& a_lhs, const Vec& a_x, Real a_scale) override;

  /**
    @brief Set a_lhs = a_a*a_x + a_b*a_y.
    @param[out] a_lhs Result
    @param[in]  a_x   First operand
    @param[in]  a_y   Second operand
    @param[in]  a_a   Scale of a_x
    @param[in]  a_b   Scale of a_y
  */
  void
  axby(Vec& a_lhs, const Vec& a_x, const Vec& a_y, Real a_a, Real a_b) override;

  /**
    @brief Scale a_lhs *= a_scale.
    @param[in,out] a_lhs   Vector
    @param[in]     a_scale Scale factor
  */
  void
  scale(Vec& a_lhs, const Real& a_scale) override;

  /**
    @brief Norm of a_rhs. For a_ord==2 this is the inner-product-induced norm sqrt(dotProduct(x,x)) so that it is
           consistent with dotProduct (a requirement of GMRES); other orders use the per-level max norm.
    @param[in] a_rhs Vector
    @param[in] a_ord Norm order
    @return Norm
  */
  Real
  norm(const Vec& a_rhs, int a_ord) override;

  /**
    @brief Set a_lhs = 0.
    @param[out] a_lhs Vector
  */
  void
  setToZero(Vec& a_lhs) override;

protected:
  /**
    @brief Reused multigrid solver (not owned).
  */
  AMRMultiGrid<LevelData<T>>* m_mg = nullptr;

  /**
    @brief Per-level AMR operators borrowed from m_mg (not owned).
  */
  Vector<AMRLevelOp<LevelData<T>>*> m_amrOps;

  /**
    @brief Per-level dx^SpaceDim, the AMR inner-product weight.
  */
  Vector<Real> m_dxScale;

  /**
    @brief Per-level mask: false where the level is covered by a finer level (excluded from the inner product).
  */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_validMask;

  /**
    @brief Per-level reduced total volume (geometry-only, vector-independent); cached on first dotProduct.
  */
  mutable Vector<Real> m_volume;

  /**
    @brief Whether m_volume has been computed and reduced. Reset to false in define() (regrid invalidates it).
  */
  mutable bool m_volumeCached = false;

  /**
    @brief Coarsest level in the solve.
  */
  int m_lbase = 0;

  /**
    @brief Finest level in the solve.
  */
  int m_lmax = 0;

  /**
    @brief Number of V-cycles per preconditioner application.
  */
  int m_numVCycles = 1;

  /**
    @brief Whether the preconditioner scratch has been allocated.
  */
  bool m_precondInit = false;

  /**
    @brief Preconditioner scratch (residual), lazily allocated on first preCond.
  */
  mutable Vec m_precondTmp;

  /**
    @brief Preconditioner scratch (V-cycle correction), lazily allocated on first preCond.
  */
  mutable Vec m_precondDcor;
};

#include <CD_NamespaceFooter.H>

#include <CD_AMRMultigridKrylovOpImplem.H>

#endif
