#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


/******************************************************************************/
/**
 * \file
 *
 * \brief Non-inline definitions for classes in FourthOrderMappedFineInterp.H
 *
 *//*+*************************************************************************/

#include "FourthOrderMappedFineInterp.H"
#include "FourthOrderMappedFineInterpF_F.H"
#include "FourthOrderMappedFineInterpSup.H"
#include "RealVect.H"
#include "InsertionSort.H"
#include "BoxIterator.H"
#include "ProblemDomain.H"
#include "FourthOrderUtil.H"
#include "LevelGridMetrics.H"

#include "NamespaceHeader.H"


/*******************************************************************************
 *
 * Class FourthOrderMappedFineInterp: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//   Non-defining constructor
/**  \param[in]  a_degree
 *                      Degree of the polynomial
 *   \param[in]  a_nComp
 *                      Number of components in \f$<U>\f$ and
 *                      \f$<JU>\f$
 *//*-----------------------------------------------------------------*/

FourthOrderMappedFineInterp::FourthOrderMappedFineInterp(const int a_degree,
                                                         const int a_nComp)
  :
  m_degree(a_degree),
  // This is hard-coded because the code for marking the cells in the stencil
  // is quite dependent on a value of 2.  An extent of 2 means that the
  // interior stencil uses a maximum of 2 extra cells in a direction.
  m_stencilExtent(2),
  m_nComp(a_nComp),
  m_willFillGhostsWithCrFnJU(false),
  m_CrFnNumGhost(-1*IntVect::Unit),
  m_CrFnStatus(CrFnStatusInvalid),
  m_isMultiblock(false)
{
  m_defined[0] = false;
  m_defined[1] = false;
  m_defined[2] = false;
  m_defined[3] = false;
}

/*--------------------------------------------------------------------*/
//   Constructor that initializes the stencils
/**  \param[in]  a_degree
 *                      Degree of the polynomial
 *   \param[in]  a_nComp
 *                      Number of components in \f$<U>\f$ and
 *                      \f$<JU>\f$
 *   \param[in]  a_refRatio
 *                      Refinement ratio in each direction
 *   \param[in]  a_h    Vector of mesh spacing in each direction on
 *                      the coarse mesh in computation space
 *//*-----------------------------------------------------------------*/

FourthOrderMappedFineInterp::FourthOrderMappedFineInterp(
  const int       a_degree,
  const int       a_nComp,
  const IntVect&  a_refRatio,
  const RealVect& a_h)
  :
  m_degree(a_degree),
  // This is hard-coded because the code for marking the cells in the stencil
  // is quite dependent on a value of 2.  An extent of 2 means that the
  // interior stencil uses a maximum of 2 extra cells in a direction.
  m_stencilExtent(2),
  m_nComp(a_nComp),
  m_willFillGhostsWithCrFnJU(false),
  m_CrFnNumGhost(-1*IntVect::Unit),
  m_CrFnStatus(CrFnStatusInvalid),
  m_isMultiblock(false)
{
  m_defined[0] = false;
  m_defined[1] = false;
  m_defined[2] = false;
  m_defined[3] = false;
  defineStencils(a_refRatio, a_h);
}

/*--------------------------------------------------------------------*/
//   Destructor
/*--------------------------------------------------------------------*/

FourthOrderMappedFineInterp::~FourthOrderMappedFineInterp()
{
  if (isStencilDefined())
    {
      delete[] m_fnXiLSx;
      delete[] m_fnGradXiLSx;
    }
}


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//   Define that initializes the stencils
/**  Stencils are selected.  The number is minimized by considering
 *   symmetries.  Then, A^T is computed for each stencil.  Also, the
 *   powers of displacments and gradients are computed for all the
 *   fine cells in a coarse cell.
 *   \param[in]  a_refRatio
 *                      Refinement ratio in each direction
 *   \param[in]  a_h    Vector of mesh spacing in each direction on
 *                      the coarse mesh in computation space.
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::defineStencils(
  const IntVect&  a_refRatio,
  const RealVect& a_h)
{
  m_refRatio = a_refRatio;
  m_mapStencilIndex.clear();

  // Number of powers.  This is number of columns in 'A'.
  const int lenN = binomial(m_degree + SpaceDim, SpaceDim);
  // Maximum number of cells in a stencil
  const int lenMMax = ((int)std::pow((float)3, SpaceDim)) + 2*SpaceDim;
  int seqIndex = 0;

//--Loop through all the stencils.

  // Mark which indices share the same dx
  m_uniqueDx[0] = 0;
  m_numUniqueDx = 1;
  for (int i = 1; i != SpaceDim; ++i)
    {
      int j = 0;
      while (j != i && a_h[i] != a_h[j]) ++j;
      m_uniqueDx[i] = (j == i) ? m_numUniqueDx++ : m_uniqueDx[j];
    }

  // Predict the number of stencils
  m_numStencil = 1;
  for (int i = 0; i != m_numUniqueDx; ++i)
    {
      int c = 0;
      for (int j = 0; j != SpaceDim; ++j)
        {
          if (m_uniqueDx[j] == i) ++c;
        }
      m_numStencil *= binomial(c + m_stencilExtent, m_stencilExtent);
    }
  
  // Set the maximum matrix size in the allocator
  m_At.getAllocator().define(lenN, lenMMax);
  // Now define the array of matrices
  m_At.define(m_numStencil);

//--Compute powers of displacements for all fine cells in a coarse cell

  RealVect fineDx;
  for (int i = 0; i != SpaceDim; ++i)
    {
      fineDx[i] = 1./m_refRatio[i];
    }
  // The actual box where we need <xi>^p and the gradients
  Box fineCellBox(IntVect::Zero, m_refRatio - IntVect::Unit);
  m_avgXipFine.define(lenN, fineCellBox);
  m_avgXiKpFine.define(lenN, fineCellBox);
  m_avgGradXipFine.define(lenN*SpaceDim, fineCellBox);
  BoxIterator fineCellBIt1(fineCellBox);
  for (fineCellBIt1.begin(); fineCellBIt1.ok(); ++fineCellBIt1)
    {
      const IntVect& fineCellIV1 = fineCellBIt1();
      RealVect rdispl;
      for (int i = 0; i != SpaceDim; ++i)
        {
          // The displacement from the center of the coarse cell
          rdispl[i] = (fineCellIV1[i] - 0.5*(m_refRatio[i] - 1))*fineDx[i];
        }
      loadAvgXipA(&m_avgXipFine(0, fineCellIV1), rdispl, fineDx);
      loadAvgXiKpA(&m_avgXiKpFine(0, fineCellIV1), rdispl, fineDx);
      loadAvgGradXip(&m_avgGradXipFine(0, fineCellIV1), rdispl, fineDx);
    }

//--Prepare matrices for each stencil

  // 'At' with constant K is temporary.  This is the transpose of matrix A
  // for solving \f$A\mathbf{x} = \mathbf{b}\f$.  I.e., the powers of the
  // displacements \f$\delta\xi\f$ to the neighbour cells.  This matrix is
  // created from XiK^p on the coarse cells.
  CHMatrix AKt(lenN, lenMMax);

  m_fnXiLSx = new CHArray<Real, SpaceDim+1, ArRangeCol>[m_numStencil];
  for (int i = 0; i != m_numStencil; ++i)
    {
      m_fnXiLSx[i].define(lenMMax, fineCellBox);
    }
  m_fnGradXiLSx = new CHArray<Real, SpaceDim+1, ArRangeCol>[m_numStencil];
  for (int i = 0; i != m_numStencil; ++i)
    {
      m_fnGradXiLSx[i].define(lenMMax*SpaceDim, fineCellBox);
    }

  // The maximum stencil extent is 1 beyond the standard interior
  // stencil (see Fig 3c for 2D in "Fourth-order conservation law example")
  const int maxStencilExtent = m_stencilExtent + 1;
  const Box maxStencilBox(-maxStencilExtent*IntVect::Unit,
                          maxStencilExtent*IntVect::Unit);
  m_ivsStencilMask.define(m_numStencil);

  // soi[] (Stencil Offset from Interior) is the displacement of the center
  // of the stencil from the standard interior stencil
  int soi[SpaceDim];

  // This describes the dependence between the dimensions.  Dimensions with
  // different dx are independent and loop to the full stencil extent.
  // Dimensions with the same dx can be interchanged so that order is not
  // important --- they are "combinations" so only consider
  // i0 >= i1 >= i2 >= i3 >= i4 >= i5
  // for all dimensions with the same dx
  const int *dimDependence[SpaceDim];
  dimDependence[0] = &m_stencilExtent;
  for (int i = 1; i != SpaceDim; ++i)
    {
      int j = i-1;
      while (j >= 0 && a_h[i] != a_h[j]) --j;
      dimDependence[i] = (j < 0) ? &m_stencilExtent : &soi[j];
    }
  D_TERM6(for (soi[0] = 0; soi[0] <= (*dimDependence[0]); ++soi[0]),
          for (soi[1] = 0; soi[1] <= (*dimDependence[1]); ++soi[1]),
          for (soi[2] = 0; soi[2] <= (*dimDependence[2]); ++soi[2]),
          for (soi[3] = 0; soi[3] <= (*dimDependence[3]); ++soi[3]),
          for (soi[4] = 0; soi[4] <= (*dimDependence[4]); ++soi[4]),
          for (soi[5] = 0; soi[5] <= (*dimDependence[5]); ++soi[5]))
    {
      // Store the index to this matrix
      m_mapStencilIndex.insert(std::make_pair(offset2Key(soi), seqIndex));

//--Mark cells in this stencil.

      // Walls appear in the positive direction and push the stencil the
      // other way.

      // Use an IntVectSet to save the stencil
      m_ivsStencilMask(seqIndex).define(maxStencilBox);
      m_ivsStencilMask(seqIndex).makeEmptyBits();
      // Inner set
      IntVect innerCenter = IntVect::Zero;
      for (int iDir = 0; iDir != SpaceDim; ++iDir)
        {
          if (soi[iDir] == 2)
            {
              innerCenter.shift(iDir, -1);
            }
        }
      m_ivsStencilMask(seqIndex) |= Box(diagShift(innerCenter, -1),
                                        diagShift(innerCenter,  1));
      // Outer set
      for (int iDir = 0; iDir != SpaceDim; ++iDir)
        {
          const int innerShiftDir = innerCenter[iDir];
          // Always add to the -ve direction
          m_ivsStencilMask(seqIndex) |=
            IntVect(IntVect::Zero).shift(iDir, innerShiftDir - 2);
          // Only to +ve direction if no wall
          if (soi[iDir] == 0)
            {
              m_ivsStencilMask(seqIndex) |=
                IntVect(IntVect::Zero).shift(iDir, innerShiftDir + 2);
            }
        }

//--Compute matrix A^T for this stencil

      int idxCenCrCell = -1;  // Seq. index in M of center coarse cell
      // Load A^T (row = powers, col = cells, and column-ordered storage)
      IVSIterator displ(m_ivsStencilMask(seqIndex));
      int iCol = 0;
      for (displ.begin(); displ.ok(); ++displ)
        {
          if (displ() == IntVect::Zero)
            {
              idxCenCrCell = iCol;
            }
          const RealVect rdispl(displ());
          const RealVect rdx(IntVect::Unit);
          int iColIdx = iCol++;
          loadAvgXipA(&m_At(seqIndex)(0, iColIdx), rdispl, rdx);
          loadAvgXiKpA(&AKt(0, iColIdx), rdispl, rdx);
        }
      // Note lenM = iCol

//--Compute Xi (A^TA) AT (constrained least squares intepolation of U).  The
//--constraint is on U, not JU for filling ghost cells.

      int ier;
      FORT_MAPLSATAINVAT(
        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_fnXiLSx[seqIndex]),
        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_fnGradXiLSx[seqIndex]),
        CHF_MATRIX(AKt),
        CHF_BOX(fineCellBox),
        CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_avgXiKpFine),
        CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_avgGradXipFine),
        CHF_CONST_INT(idxCenCrCell),
        CHF_CONST_INT(iCol),
        CHF_INT(ier));
      if (ier != 0)
        {
          pout() << "Lapack error code: " << ier << std::endl;
          MayDay::Error("Failure in matrix inversion in FORT_MAPINVERSEATA");
        }

      ++seqIndex;
    }  // Loop over stencils
  // Sanity check
  CH_assert(seqIndex == m_numStencil);

  m_defined[0] = true;
}

/*--------------------------------------------------------------------*/
//   Define that allocates data for coarse representations of the fine
//   grid
/**  Note that this data hangs around forever.
 *   \param[in]  a_FnGrid
 *                      Layout of the fine boxes.
 *   \param[in]  a_CrGrid
 *                      Layout of the coarse boxes.
 *   \param[in]  a_CrProblemDomain
 *                      Problem domain on the coarse mesh
 *   \param[in]  a_FnInterpRadVec
 *                      Number of fine cells to fill in each direction
 *                      (default 0)
 *   \param[in]  a_willFillGhostsWithCrFnJU
 *                      T - \<JU\> will be provided on the CrFn mesh
 *                          for filling ghost cells.  This implies
 *                          that \<U\> needs to be computed and
 *                          therefore one extra ghost is needed on
 *                          the CrFn layout (default true).  You
 *                          don't _have_ to fill with \<JU\> if set
 *                          true but you definitely cannot if set to
 *                          false.
 *                      F - Memory for \<JU\> on the CrFn mesh is not
 *                          even allocated
 *   \param[in]  a_mbgeoPtr
 *                      Information about the multiblock geometry.  If
 *                      provided, this triggers a hack where the CrFn
 *                      mesh for \<JU\> is actually the coarse mesh.
 *                      This allows us to compute \<U\> before
 *                      invoking a special multiblock copyTo.  Much of
 *                      this hack, relating to time interpolation
 *                      instead of space interpolation, is within this
 *                      class simply because the relevant data
 *                      structures are here.  And since it's probably
 *                      just a temporary hack.
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::defineGrids(
  const DisjointBoxLayout&         a_FnGrid,
  const DisjointBoxLayout&         a_CrGrid,
  const ProblemDomain&             a_CrProblemDomain,
  const IntVect&                   a_FnInterpRadVec,
  const bool                       a_willFillGhostsWithCrFnJU,
  const MultiBlockLevelGeom *const a_mbgeoPtr)
{
  CH_TIME("FourthOrderMappedFineInterp::defineGrids");
  m_isMultiblock = (a_mbgeoPtr != 0);
  // Ghost radius
  m_FnInterpRadVec = a_FnInterpRadVec;
  D_TERM6(m_CrInterpRadVec[0] = (a_FnInterpRadVec[0] + m_refRatio[0] - 1)/m_refRatio[0];,
          m_CrInterpRadVec[1] = (a_FnInterpRadVec[1] + m_refRatio[1] - 1)/m_refRatio[1];,
          m_CrInterpRadVec[2] = (a_FnInterpRadVec[2] + m_refRatio[2] - 1)/m_refRatio[2];,
          m_CrInterpRadVec[3] = (a_FnInterpRadVec[3] + m_refRatio[3] - 1)/m_refRatio[3];,
          m_CrInterpRadVec[4] = (a_FnInterpRadVec[4] + m_refRatio[4] - 1)/m_refRatio[4];,
          m_CrInterpRadVec[5] = (a_FnInterpRadVec[5] + m_refRatio[5] - 1)/m_refRatio[5];)
  m_willFillGhostsWithCrFnJU = a_willFillGhostsWithCrFnJU;
  //** The following requires an int but we generally assume that refRatio is
  //** a vector
  m_CrFnGrid = DisjointBoxLayout();
  coarsen(m_CrFnGrid, a_FnGrid, m_refRatio[0]);
  int numExtraGhosts = m_stencilExtent;
  if (a_willFillGhostsWithCrFnJU && !isMultiblock())
    {
      // The CrFn JU needs to be 1 ghost cell larger than U because when we are
      // filling ghosts we need to first compute U.  Note however that
      // m_CrFnLevU is allocated with the same number of ghosts so that we can
      // use the same copier.
      // However, if multiblock, we *always* copy U from the coarse grid to the
      // coarsened fine mesh.  Hence, we don't need this extra ghost.  THIS MAY
      // CHANGE ONCE THE WORKAROUND FOR TIME INTERPOLATION IN MULTIBLOCK IS
      // RESOLVED!!!
      numExtraGhosts += 1;
    }
  m_CrFnNumGhost = (m_CrInterpRadVec + numExtraGhosts*IntVect::Unit);
  m_CrFnLevU.define(m_CrFnGrid, m_nComp, m_CrFnNumGhost);
  if (a_willFillGhostsWithCrFnJU)
    {
      if (a_mbgeoPtr)
        {
          // Special hack for multiblock.  Since this covers the entire coarse
          // grid, only 1 ghost cell is need to allow for computing <U>
          m_CrFnLevJU.define(a_CrGrid, m_nComp, IntVect::Unit);
          // Define the multiblock copier
          // m_mbcp.define(a_mbgeoPtr, m_CrFnGrid, m_degree + 1);
          m_mbcp.define(a_mbgeoPtr, m_CrFnGrid, m_CrFnNumGhost, m_degree + 1);
          m_mbcp.defineVector(); // added by petermc, 30 Aug 2012
        }
      else
        {
          // Assuming single block
          m_CrFnLevJU.define(m_CrFnGrid, m_nComp, m_CrFnNumGhost);
        }
    }
  m_copier.define(a_CrGrid, m_CrFnGrid, a_CrProblemDomain, m_CrFnNumGhost);
  m_defined[1] = true;
  m_defined[2] = false;              // Invalidates CFInterface.
  m_defined[3] = false;              // Invalidates m_CrFnLevJ.,
  m_CrFnStatus = CrFnStatusInvalid;  // Invalidate CrFn Data
}

/*--------------------------------------------------------------------*/
//   Define that saves the coarse <J> on the coarsened-fine layout
/**  This is useful if we are given \<JU\> on the coarsened-fine
 *   and need to convert it to </U\> before interpolating to fill any
 *   ghost cells.  The grid must be defined so we know the coarsened-
 *   fine layout.  It is necessary to call this routine before using
 *   'fillGhostsFromCrJU'
 *   \param[in]  a_CrJ  \<J\> from the coarse layout
 *   \note
 *   <ul>
 *     <li> You must have indicated that you will fill ghosts with
 *          CrFnJU when defineGrids was called.  This ensures there
 *          are sufficient layers of ghosts to calculate \<U\> on all
 *          cells in a stencil.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::defineCrFnJ(const LevelData<FArrayBox>& a_CrJ)
{
  CH_TIME("FourthOrderMappedFineInterp::defineCrFnJ");
  CH_assert(isGridDefined());
  CH_assert(m_willFillGhostsWithCrFnJU);
  if (isMultiblock())
    {
      // Special hack for multiblock.
      LevelData<FArrayBox> *const fullCrJPtr =
        const_cast<LevelData<FArrayBox>*>(&a_CrJ);
      aliasLevelData(m_CrFnLevJ, fullCrJPtr, a_CrJ.interval());
    }
  else
    {
      m_CrFnLevJ.define(m_CrFnLevJU.getBoxes(), 1, m_CrFnLevJU.ghostVect());
      a_CrJ.copyTo(m_CrFnLevJ, m_copier);
    }
  m_defined[3] = true;
}

/*--------------------------------------------------------------------*/
//   Define that initializes the coarse-fine interface for filling
//   ghosts
/**  \param[in]  a_CrMetricsPtr
 *                      Used to get problem domain on the coarse mesh
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::defineCFInterface(
  LevelGridMetrics *const a_CrMetricsPtr)
{
  CH_TIME("FourthOrderMappedFineInterp::defineCFInterface");
  // Need some stencil data
  CH_assert(isStencilDefined());
  // The grid has to be defined to use this functionality
  CH_assert(isGridDefined());

//--Initializations

  DataIterator dit = m_CrFnGrid.dataIterator();

  bool isSorted = m_CrFnGrid.isSorted();

  // Re-definitions
  m_numCrFnFillGhost = 0;
  m_CrFnFillInterp.define(m_CrFnGrid);
  m_ArIdxFillInterp.define(m_CrFnGrid);
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_ArIdxFillInterp[dit].define(grow(m_CrFnGrid[dit], m_CrInterpRadVec), 1);
    }

  // A box which will determine whether a given box adjoins a periodic boundary
  Box crPeriodicTestBox;
  bool isPeriodicDomain = false;
  if (!isMultiblock())
    {
      const ProblemDomain& crProblemDomain = a_CrMetricsPtr->problemDomain();
      crPeriodicTestBox = crProblemDomain.domainBox();
      if (crProblemDomain.isPeriodic())
        {
          isPeriodicDomain = true;
          for (int iDir = 0; iDir < SpaceDim; ++iDir)
            {
              if (crProblemDomain.isPeriodic(iDir))
                {
                  crPeriodicTestBox.grow(iDir, -1);
                }
            }
        }
    }

//--Mark all coarse cells on the coarsened fine layout that will be used to fill
//--ghost cells on the fine mesh

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box crFnBox = m_CrFnGrid[dit];
      const ProblemDomain& crBlockDomain =
        a_CrMetricsPtr->blockDomain(crFnBox, m_CrInterpRadVec[0]);

      const Box crFnGhostBox =
        grow(crFnBox, m_CrInterpRadVec) & crBlockDomain;
      IntVectSet& crFnFillGhostsIVS = m_CrFnFillInterp[dit];

      // Initialize the IVS as the whole crFnGhostBox
      crFnFillGhostsIVS.define(crFnGhostBox);

      // Check if any of the ghost cells extend across a periodic domain
      bool notAcrossPeriodic = true;
      if (isPeriodicDomain)
        {
          notAcrossPeriodic =
            crBlockDomain.domainBox().contains(crFnGhostBox);
        }

      // Now we iterate over all the boxes in m_CrFnGrid and subtract off all
      // the real cells leaving only the invalid ghost cells
      LayoutIterator lit = m_CrFnGrid.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
        {
          const Box crFnRealBox = m_CrFnGrid[lit];

          if (notAcrossPeriodic && isSorted &&
              (crFnRealBox.bigEnd(0) < crFnGhostBox.smallEnd(0)))
            {
              // Nothing from this box
              continue;
            }

          if (notAcrossPeriodic && isSorted &&
              (crFnRealBox.smallEnd(0) > crFnGhostBox.bigEnd(0)))
            {
              // Can break out of loop, since we know that the smallEnd of all
              // the remaining boxes are lexigraphically beyond the ghost box.
              break;
            }

          crFnFillGhostsIVS -= crFnRealBox;

          // We also need to remove any periodic images of real cells
          if (isPeriodicDomain
              && !crPeriodicTestBox.contains(crFnRealBox)
              && !crPeriodicTestBox.contains(crFnGhostBox))
            {
              ShiftIterator shiftIt = crBlockDomain.shiftIterator();
              IntVect shiftMult(crBlockDomain.domainBox().size());
              Box shiftedBox(crFnRealBox);
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect = shiftMult*shiftIt();
                  shiftedBox.shift(shiftVect);
                  crFnFillGhostsIVS -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                }
            }
        }  // Second loop over all boxes

      // Set a sequential index for all coarse ghost cells in the IVS:
      BaseFab<int>& arIdxFillInterp = m_ArIdxFillInterp[dit];
      arIdxFillInterp.setVal(-1);
      IVSIterator crFnFillGhostsIt(crFnFillGhostsIVS);
      for (crFnFillGhostsIt.begin(); crFnFillGhostsIt.ok(); ++crFnFillGhostsIt)
        {
          arIdxFillInterp(crFnFillGhostsIt()) = m_numCrFnFillGhost++;
        }
    }  // First loop over boxes on the processor

//--Set a box of fine ghost cells to be filled for each coarse cell

  m_FnBoxFillInterp.define(m_numCrFnFillGhost);
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box fnGhostBox =
        grow(refine(m_CrFnGrid[dit], m_refRatio), m_FnInterpRadVec);
      BaseFab<int>& arIdxFillInterp = m_ArIdxFillInterp[dit];
      IVSIterator crFnFillGhostsIt(m_CrFnFillInterp[dit]);
      for (crFnFillGhostsIt.begin(); crFnFillGhostsIt.ok(); ++crFnFillGhostsIt)
        {
          const IntVect crFnFillGhostIV = crFnFillGhostsIt();
          int arIdx = arIdxFillInterp(crFnFillGhostIV);
          m_FnBoxFillInterp(arIdx) = fnGhostBox &
            refine(Box(crFnFillGhostIV, crFnFillGhostIV), m_refRatio);
        }
    }

  m_defined[2] = true;
  m_CrFnStatus = CrFnStatusInvalid;  // Invalidate CrFn Data
}

/*--------------------------------------------------------------------*/
//   Interpolate from a coarse to fine level on a mapped grid
/**  The copier will be used to get the coarse data onto a coarse
 *   representation of the fine mesh.  If the grids have been defined,
 *   the coarse data will be placed there.  Otherwise it will be
 *   placed in a temporary and deallocated.
 *   \param[out] a_FnLevJU
 *                      \<JU\> interpolated on the fine grid
 *   \param[in]  a_FnLevJ
 *                      The metrics \<J\> on the fine grid
 *   \param[in]  a_CrLevU
 *                      \<U\> on the coarse grid
 *   \param[in]  a_CrLevJU
 *                      \<JU\> on the coarse grid
 *   \param[in]  a_CrMetricsPtr
 *                      Used to get problem domain on the coarse grid
 *  \param[in]  a_vectorIntv
 *                      An interval in components of U consisting of
 *                      vector data
 *  \note
 *   <ul>
 *     <li> Stencils must have been defined.
 *     <li> If the grids have been defined, the coarse data will be
 *          copied there.
 *     <li> If the grids have not been defined, the coarse data will
 *          be copied to a temporary and deleted upon exit.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::interpToFine(
  LevelData<FArrayBox>&       a_FnLevJU,
  const LevelData<FArrayBox>& a_FnLevJ,
  const LevelData<FArrayBox>& a_CrLevU,
  const LevelData<FArrayBox>& a_CrLevJU,
  LevelGridMetrics *const     a_CrMetricsPtr,
  const Interval&             a_vectorIntv)
{
  CH_TIME("FourthOrderMappedFineInterp::interpToFine");
  CH_assert(isStencilDefined());
  //**FIXME m_nComp and m_CrFnNumGhost are defined, maybe use those?
  const int nComp   = a_FnLevJU.nComp();
  const IntVect numGhost(m_stencilExtent*IntVect::Unit);

  LevelData<FArrayBox> *CrFnLevU;
  LevelData<FArrayBox> *CrFnLevJU;
  Copier *copier;
  // Modifiable
  DisjointBoxLayout *CrFnGridPtr = 0;
  // Unmodifiable if we only reference
  const DisjointBoxLayout *CrFnGrid = 0;
  if (!isGridDefined())
    {
      // Proceeding without a grid defined is only possible for single-block
      CH_assert(!isMultiblock());
      CrFnGridPtr = new DisjointBoxLayout;
      CrFnLevU = new LevelData<FArrayBox>;
      //** The following requires an int but we generally assume that refRatio
      //** is a vector
      coarsen(*CrFnGridPtr, a_FnLevJU.getBoxes(), m_refRatio[0]);
      CrFnGrid = CrFnGridPtr;
      CrFnLevU->define(*CrFnGrid, nComp, numGhost);
      if (!isMultiblock()) // Special hack for multiblock below.
        {
          CrFnLevJU = new LevelData<FArrayBox>;
          CrFnLevJU->define(*CrFnGrid, nComp, numGhost);
        }
      const ProblemDomain& crProblemDomain = a_CrMetricsPtr->problemDomain();
      copier = new Copier(a_CrLevU.getBoxes(),
                          *CrFnGrid,
                          crProblemDomain,
                          numGhost);
    }
  else
    {
      CrFnGrid = &m_CrFnGrid;
      CrFnLevU = &m_CrFnLevU;
      CrFnLevJU = &m_CrFnLevJU;
      copier = &m_copier;
    }
  if (isMultiblock())
    {
      // Special hack for multiblock.
      CrFnLevJU = new LevelData<FArrayBox>;
      CrFnLevJU->define(*CrFnGrid, m_nComp, m_CrFnNumGhost);
    }

  // Check for proper nesting
  CH_assert(properNesting(a_CrLevJU.getBoxes(),
                          *CrFnGrid,
                          a_CrMetricsPtr));

  // Copy data to the coarsened fine grid
  // For multiblock, <JU> can be copied using the regular copier because we only
  // need it in valid cells.  To copy <U>, we use the presetCrFnLevU routine
  // which copies to m_CrFnLevU.  For multiblock, the grid must be defined (see
  // assertion above) so this is okay.
  if (isMultiblock())
    {
      presetCrFnLevU(a_CrLevU, a_vectorIntv);
    }
  else
    {
      a_CrLevU.copyTo(*CrFnLevU, *copier);
    }
  a_CrLevJU.copyTo(*CrFnLevJU, *copier);

  // Now continue with the least-squares interpolation
  interpCrFnToFine(
    a_FnLevJU,
    a_FnLevJ,
    *CrFnGrid,
    *CrFnLevU,
    *CrFnLevJU,
    a_CrMetricsPtr);

  // Clean up memory if we allocated it
  if (!isGridDefined())
    {
      delete CrFnGridPtr;
      delete CrFnLevU;
      delete CrFnLevJU;
      delete copier;
    }
  else if (isMultiblock())
    {
      // Special hack for multiblock.
      delete CrFnLevJU;
    }
}

/*--------------------------------------------------------------------*/
//   Fill the invalid ghost cells on a mapped grid (OBSOLETE)
/**  \param[out] a_FnLevU
 *                      \<U\> interpolated on the fine grid
 *   \param[out] a_FnLevJU
 *                      \<JU\> interpolated on the fine grid.  Can
 *                      have less ghosts than a_FnLevU -- only updated
 *                      for cells that exist
 *   \param[in]  a_FnLevJ
 *                      The metrics \<J\> on the fine grid
 *   \param[in]  a_CrLevData
 *                      \<U\> or \<JU\> on the coarse grid
 *   \param[in]  a_CrDataType
 *                      Type of data, \<U\> or \<JU\> on the coarse
 *                      grid.
 *                      "U"  = \<U\>
 *                      "JU" = \<JU\>
 *   \param[in]  a_CrProblemDomain
 *                      Problem domain on the coarse grid
 *
 *   \note
 *   <ul>
 *     <li> Stencil, grids, and the coarse-fine interface must have
 *          been defined.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::fillGhosts(
  LevelData<FArrayBox>&       a_FnLevU,
  LevelData<FArrayBox>&       a_FnLevJU,
  const LevelData<FArrayBox>& a_FnLevJ,
  const LevelData<FArrayBox>& a_CrLevData,
  const char *const           a_CrDataType,
  const ProblemDomain&        a_CrProblemDomain)
{
  // This routine is becoming cumbersome to maintain and is unused.  In the
  // hope of deleting it soon...
  MayDay::Error("Don't use this FourthOrderMappedFineInterp::fillGhosts");
  CH_assert(isStencilDefined());
  CH_assert(isGridDefined());
  CH_assert(isCFInterfaceDefined());
  CH_assert(a_CrDataType != NULL);

//**Disabled, only check during a regrid
//   // Check for proper nesting
//   CH_assert(properNesting(a_CrLevJU.getBoxes(),
//                           *CrFnGrid,
//                           a_CrProblemDomain));

  // Copy data to the coarsened fine grid
  switch (a_CrDataType[0])
    {
    // case 'U':
    //   presetCrFnLevU(a_CrLevData);
    //   break;
    // case 'J':  // Meaning "JU"
    //   presetCrFnLevJU(a_CrLevData);
    //   break;
    default:
      MayDay::Error("FourthOrderMappedFineInterp::fillGhosts: invalid argument "
                    "for a_CrDataType");
      break;
    }

  // Now continue with the least-squares interpolation
  // fillGhosts(
  //   a_FnLevU,
  //   a_FnLevJU,
  //   a_FnLevJ,
  //   a_CrProblemDomain);
}

/*--------------------------------------------------------------------*/
//   Fill invalid ghost cells from preset coarsened-fine data on a
//   mapped grid
/**  To call this, one of the "::presetCrFnLev*" routines must have
 *   been used to get the coarse data on a coarse representation of
 *   the fine mesh.  Also, the coarse-fine interface must have been
 *   identified.
 *   \param[out] a_FnLevU
 *                      \<U\> interpolated on the fine grid
 *   \param[out] a_FnLevJU
 *                      \<JU\> interpolated on the fine grid.  Can
 *                      have less ghosts than a_FnLevU -- only updated
 *                      for cells that exist
 *   \param[in]  a_FnLevJ
 *                      The metrics \<J\> on the fine grid
 *   \param[in]  a_CrMetricsPtr
 *                      Used to get problem domain on the coarse grid
 *
 *  \note
 *   <ul>
 *     <li> Stencils must have been defined.
 *     <li> The internal coarsened fine representation of \<U\> must
 *          have been set either by 'fillGhosts()' or using
 *          'presetCrFnLevU()' or 'presetCrFnLevJU()' to set it
 *          externally.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::fillGhosts(
  LevelData<FArrayBox>&       a_FnLevU,
  LevelData<FArrayBox>&       a_FnLevJU,
  const LevelData<FArrayBox>& a_FnLevJ,
  LevelGridMetrics *const     a_CrMetricsPtr)
{
  CH_TIME("FourthOrderMappedFineInterp::fillGhosts");
  CH_assert(isStencilDefined());
  CH_assert(isGridDefined());
  CH_assert(isCFInterfaceDefined());
  if (m_CrFnStatus == CrFnStatusHaveJU)
    {
      // Deferred calculations of <U> is only possible for single-block
      CH_assert(!isMultiblock());
      const ProblemDomain& crProblemDomain = a_CrMetricsPtr->problemDomain();
      // This means we have <JU>, but not <U> so compute <U>
      for (DataIterator dit = m_CrFnLevU.dataIterator(); dit.ok(); ++dit)
        {
          // Compute CrFn <U> in the domain
          FArrayBox& crFnFabU = m_CrFnLevU[dit];
          const FArrayBox& crFnFabJU = m_CrFnLevJU[dit];
          const Box UBox = grow(crFnFabJU.box(), -1) & crProblemDomain;
          cellFGToCellF(crFnFabU,
                        crFnFabJU,
                        m_CrFnLevJ[dit],
                        UBox,
                        crProblemDomain,
//**CHECKME Why do we need 'true' here?  Why use one-sided gradients of <J>
//**        near boundaries?
                        true);
        }
      m_CrFnStatus |= CrFnStatusHaveU;
    }
  CH_assert(m_CrFnStatus & CrFnStatusHaveU);

  // Number of ghosts required for full interpolation stencil
  int numGhostStencil = m_CrInterpRadVec[0] + m_stencilExtent;

  // Loop over the boxes and call fillGhosts
  for (DataIterator dit = m_CrFnGrid.dataIterator(); dit.ok(); ++dit)
    {
      const ProblemDomain& crBlockDomain =
        a_CrMetricsPtr->blockDomain(m_CrFnGrid[dit], numGhostStencil);
      fillGhosts(a_FnLevU[dit],
                 a_FnLevJU[dit],
                 a_FnLevJ[dit],
                 m_CrFnFillInterp[dit],
                 m_ArIdxFillInterp[dit],
                 m_CrFnLevU[dit],
                 crBlockDomain);
    }
}

/*--------------------------------------------------------------------*/
//   Interpolate from a coarsened-fine to fine level on a mapped grid
/**  To call this, the copier must have already been used to get the
 *   coarse data on a coarse representation of the fine mesh.
 *   \param[out] a_FnLevJU
 *                      \<JU\> interpolated on the fine grid
 *   \param[in]  a_FnLevJ
 *                      The metrics \<J\> on the fine grid
 *   \param[in]  a_CrFnGrid
 *                      A coarse representation of the fine grid.
 *                      This is required instead of m_CrFnGrid
 *                      directly because the grid may have been
 *                      temporarily defined in routine interpToFine
 *   \param[in]  a_CrFnLevU
 *                      \<U\> on a coarse representation of the fine
 *                      grid
 *   \param[in]  a_CrFnLevJU
 *                      \<JU\> on a coarse representation of the fine
 *                      grid
 *   \param[in]  a_CrMetricsPtr
 *                      Used to get problem domain on the coarse grid
 *
 *   \note
 *   <ul>
 *     <li> Stencils must have been defined.
 *     <li> The internal coarsened fine grids are unused and do not
 *          need to be defined.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::interpCrFnToFine(
  LevelData<FArrayBox>&       a_FnLevJU,
  const LevelData<FArrayBox>& a_FnLevJ,
  const DisjointBoxLayout&    a_CrFnGrid,
  const LevelData<FArrayBox>& a_CrFnLevU,
  const LevelData<FArrayBox>& a_CrFnLevJU,
  LevelGridMetrics *const     a_CrMetricsPtr) const
{
  CH_TIME("FourthOrderMappedFineInterp::interpCrFnToFine");
  CH_assert(isStencilDefined());

  // Simply loop over the boxes and call interpolate
  for (DataIterator dit = a_CrFnGrid.dataIterator(); dit.ok(); ++dit)
    {
      const Box& crFnBox = a_CrFnGrid[dit];
      const ProblemDomain& crBlockDomain =
        a_CrMetricsPtr->blockDomain(crFnBox, m_stencilExtent);
      int ier = interpToFine(a_FnLevJU[dit],
                             a_FnLevJ[dit],
                             crFnBox,
                             a_CrFnLevU[dit],
                             a_CrFnLevJU[dit],
                             crBlockDomain);
      if (ier)
        {
          MayDay::Error("Constrained least-squares interpolation failed");
        }
    }
}

/*--------------------------------------------------------------------*/
//   Interpolate from a coarse to fine box on a mapped grid
/**  An appropriate stencil is selected depending on the proximity of
 *   a coarse cell to a boundary.  Data required for the least-squares
 *   problem is packed for the stencil before handing off to a Fortran
 *   routine.
 *   \param[out] a_FnFabJU
 *                      Interpolation of \<JU\> on the fine cells
 *   \param[in]  a_FnFabJ
 *                      Metrics on the fine mesh
 *   \param[in]  a_CrBox
 *                      Coarse box to interpolate from, a coarsened
 *                      representation of the fine box
 *   \param[in]  a_CrFabU
 *                      Coarse \<U\> on the box.  Mapped to LSbl for
 *                      each LS problem.
 *   \param[in]  a_CrFabJU
 *                      Coarse \<JU\> on the box.  This is the
 *                      constraint 'LSd' for each coarse cell.
 *   \param[in]  a_CrProblemDomain
 *                      The problem domain as specified on the coarse
 *                      level
 *   \return            Error status.  >0 means a Lapack routine
 *                      returned an error code.
 *                       0 - success
 *                      >0 - failure
 *//*-----------------------------------------------------------------*/

int
FourthOrderMappedFineInterp::interpToFine(
  FArrayBox&           a_FnFabJU,
  const FArrayBox&     a_FnFabJ,
  const Box&           a_CrBox,
  const FArrayBox&     a_CrFabU,
  const FArrayBox&     a_CrFabJU,
  const ProblemDomain& a_CrProblemDomain) const
{
  CH_TIME("FourthOrderMappedFineInterp::interpToFine");
  // Some definitions
  const int lenMMax = m_At(0).size(1);
  const int nComp   = a_FnFabJU.nComp();

  // Some temporaries
  CHArray<Real, 2> LSbl(lenMMax, nComp);
  Vector<Real> LSd(nComp);
  int iStencil;             // Stencil index
  CoordTransform ctStoP;    // Transformation from stencil coordinate sytem
                            // to physical coordinate system
  CoordTransform ctPtoS;    // Transformation from physical coordinate sytem
                            // to stencil coordinate system
  IntVect offset;           // Offset between (0,0) in S and P coordinates
  IntVect fineSmlEnd;       // What to add to a fine cell location transformed
                            // from S to P to get the actual location in Z

  int ierBox = 0;  // Accumulated over the box
  // For each coarse cell
  BoxIterator bIt(a_CrBox);
  for (bIt.begin(); bIt.ok(); ++bIt)
    {
      IntVect center = bIt();

      // Select a stencil and get indexes and transformation info
      selectStencil(center, a_CrProblemDomain, iStencil, ctStoP, ctPtoS, offset,
                    fineSmlEnd);

      // Get the refinement ratio in the stencil coordinate system
      IntVect refRatioS = absolute(ctPtoS.transform(m_refRatio));

      // Box of fine cells to update (in stencil space)
      Box fineBoxInCrS(IntVect::Zero, refRatioS - IntVect::Unit);

      // Pack b = <U>_{j^l} according to the stencil
      int iN = 0;
      IVSIterator displ(m_ivsStencilMask(iStencil));
      for (displ.begin(); displ.ok(); ++displ)
        {
          IntVect tdispl = ctStoP.transform(displ());
          tdispl += center;
          for (int iComp = 0; iComp != nComp; ++iComp)
            {
              LSbl(iN, iComp) = a_CrFabU(tdispl, iComp);
            }
          ++iN;
        }

      // Pack d = <JU>_{i^l}
      for (int iComp = 0; iComp != nComp; ++iComp)
        {
          LSd[iComp] = a_CrFabJU(center, iComp);
        }

      // Stencil size
      const int lenM = m_ivsStencilMask(iStencil).numPts();

      int ier;
      FORT_MAPLSINTERP(
        CHF_FRA(a_FnFabJU),
        CHF_CONST_FRA(a_FnFabJ),
        CHF_CONST_I1D(ctStoP.xctm(), 2*SpaceDim),
        CHF_BOX(fineBoxInCrS),
        CHF_CONST_INTVECT(fineSmlEnd),
        CHF_CONST_RCHARRAY(2, LSbl),
        CHF_CONST_R1D(&LSd[0], nComp),
        CHF_CONST_MATRIX(m_At(iStencil)),
        CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_avgXipFine),
        CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_avgGradXipFine),
        CHF_CONST_INT(lenM),
        CHF_INT(ier));
      if (ier != 0)
        {
          if (ierBox == 0)
            {
              ierBox = -1;
              MayDay::Warning("Possible failure in constrained least-squares "
                              "computation in FORT_MAPLSINTERP");
              pout() << "Error code: " << ier << std::endl;
            }
          ierBox = std::max(ierBox, ier);
          // But continue anyways
        }
    }  // Loop over coarse cells
  return ierBox;
}

/*--------------------------------------------------------------------*/
//   Fill any ghost cells on a finer mesh
/**  An appropriate stencil is selected depending on the proximity of
 *   a coarse cell to a boundary.  Data required for the least-squares
 *   problem is packed for the stencil before handing off to a Fortran
 *   routine.
 *   \param[out] a_FnFabU
 *                      Interpolation of \<U\> into the invalid fine
 *                      ghost cells
 *   \param[out] a_FnFabJU
 *                      Interpolation of \<JU\> into the invalid fine
 *                      ghost cells.  Can have less ghosts than
 *                      a_FnLevU -- only updated for cells that exist
 *   \param[in]  a_FnFabJ
 *                      Metrics on the fine mesh
 *   \param[in]  a_CrCells
 *                      Coarse cells underlying the fine ghost cells
 *                      in this box
 *   \param[in]  a_CrCellsArIdx
 *                      A basefab giving an array index for all the
 *                      cells in a_CrCells.  This is the location of
 *                      any auxiliary data associated with the coarse
 *                      cell.
 *   \param[in]  a_CrFabU
 *                      Coarse \<U\> on the box.  Mapped to LSbl for
 *                      each LS problem.
 *   \param[in]  a_CrProblemDomain
 *                      The problem domain as specified on the coarse
 *                      level
 *   \note
 *   <ul>
 *     <li> The coarsened fabs and IVS must exist on a coarsened
 *          representation of the fine mesh
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::fillGhosts(
  FArrayBox&           a_FnFabU,
  FArrayBox&           a_FnFabJU,
  const FArrayBox&     a_FnFabJ,
  const IntVectSet&    a_CrCells,
  const BaseFab<int>&  a_CrCellsArIdx,
  const FArrayBox&     a_CrFabU,
  const ProblemDomain& a_CrProblemDomain) const
{
  CH_TIME("FourthOrderMappedFineInterp::fillGhosts");
  // There is a requirement here that J and JU have the same number of ghost
  // cells.  Both boxes should have equal sizes.
  CH_assert(a_FnFabJU.box() == a_FnFabJ.box());

  // Some definitions
  const int lenMMax = m_At(0).size(1);
  const int nComp   = a_FnFabJU.nComp();

  // Some temporaries
  CHArray<Real, 2> LSbl(lenMMax, nComp);
  int iStencil;             // Stencil index
  CoordTransform ctStoP;    // Transformation from stencil coordinate sytem
                            // to physical coordinate system
  CoordTransform ctPtoS;    // Transformation from physical coordinate sytem
                            // to stencil coordinate system
  IntVect offset;           // Offset between (0,0) in S and P coordinates
  IntVect fineSmlEnd;       // What to add to a fine cell location transformed
                            // from S to P to get the actual location in Z

  int ierBox = 0;  // Accumulated over the box
  // For each coarse cell
  IVSIterator crCellIt(a_CrCells);
  for (crCellIt.begin(); crCellIt.ok(); ++crCellIt)
    {
      IntVect center = crCellIt();

      // Select a stencil and get indexes and transformation info
      selectStencil(center, a_CrProblemDomain, iStencil, ctStoP, ctPtoS, offset,
                    fineSmlEnd);

      // Pack b = <U>_{j^l} according to the stencil
      int iN = 0;
      IVSIterator displ(m_ivsStencilMask(iStencil));
      for (displ.begin(); displ.ok(); ++displ)
        {
          IntVect tdispl = ctStoP.transform(displ());
          tdispl += center;
          for (int iComp = 0; iComp != nComp; ++iComp)
            {
              LSbl(iN, iComp) = a_CrFabU(tdispl, iComp);
            }
          ++iN;
        }

      // Get the box of fine cells to update (in stencil space)
      // Start with the box in physical space
      const int iArIdx = a_CrCellsArIdx(center);
      Box fineBoxInCrP = m_FnBoxFillInterp(iArIdx);
      // Get indicies w.r.t. small corner of coarse cell
      const IntVect loP(fineBoxInCrP.smallEnd() - center*m_refRatio);
      const IntVect hiP(fineBoxInCrP.bigEnd() - center*m_refRatio);
      // Transform to stencil space
      const IntVect loPtoS(ctPtoS.transform(loP - offset));
      const IntVect hiPtoS(ctPtoS.transform(hiP - offset));
      // Find the small and big corners in stencil space
      const IntVect loS(min(loPtoS, hiPtoS));
      const IntVect hiS(max(loPtoS, hiPtoS));
      Box fineBoxInCrS(loS, hiS);

      // Stencil size
      const int lenM = m_ivsStencilMask(iStencil).numPts();

      int ier;
      FORT_MAPLSFILLGHOSTS(
        CHF_FRA(a_FnFabU),
        CHF_FRA(a_FnFabJU),
        CHF_CONST_FRA(a_FnFabJ),
        CHF_CONST_I1D(ctStoP.xctm(), 2*SpaceDim),
        CHF_BOX(fineBoxInCrS),
        CHF_CONST_INTVECT(fineSmlEnd),
        CHF_CONST_RCHARRAY(2, LSbl),
        CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_fnXiLSx[iStencil]),
        CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_fnGradXiLSx[iStencil]),
        CHF_CONST_INT(lenM),
        CHF_INT(ier));
      if (ier != 0)
        {
          if (ierBox == 0)
            {
              ierBox = -1;
              MayDay::Warning("Possible failure in least-squares computation "
                              "in FORT_MAPLSFILLGHOST");
              pout() << "Error code: " << ier << std::endl;
            }
          ierBox = std::max(ierBox, ier);
          // But continue anyways
        }
    }  // Loop over coarse cells
}

/*--------------------------------------------------------------------*/
//   Preset the coarsened-fine \<U\> from coarse level data
/**  \param[in]  a_CrU  \<U\> on the coarse level
 *   \param[in]  a_vectorIntv
 *                      An interval in components of U consisting of
 *                      vector data
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::presetCrFnLevU(
  const LevelData<FArrayBox>& a_CrLevU,
  const Interval&             a_vectorIntv)
{
  CH_TIME("FourthOrderMappedFineInterp::presetCrFnLevU");
  CH_assert(isGridDefined());
  m_CrFnStatus |= CrFnStatusHaveU;
  // Normal copy
  a_CrLevU.copyTo(m_CrFnLevU, m_copier);
  // Across multiblock boundaries
  if (isMultiblock())
    {
      //**FIXME Grumble... vector data in between scalars.  Should we copy it or
      //**      change the indices in LGintergrator.H?
      const Interval Uintv = m_CrFnLevU.interval();
      if (a_vectorIntv.size() == 0)
        { // In this case, all components are scalar, none vector.
          m_mbcp.copyTo(a_CrLevU, m_CrFnLevU);
        }
      else
        {
          // Ugly ... ow ow ow!
          if (Uintv.begin() < a_vectorIntv.begin())
            { // scalar component indices before a_vectorIntv
              Interval scalarIntv(Uintv.begin(), a_vectorIntv.begin()-1);
              m_mbcp.copyTo(a_CrLevU, m_CrFnLevU, scalarIntv);
            }
          m_mbcp.copyToVector(a_CrLevU, m_CrFnLevU, a_vectorIntv);
          if (a_vectorIntv.end() < Uintv.end())
            { // scalar component indices after a_vectorIntv
              Interval scalarIntv(a_vectorIntv.end()+1, Uintv.end());
              m_mbcp.copyTo(a_CrLevU, m_CrFnLevU, scalarIntv);
            }
        }
    }
}

// Should not be used
// /*--------------------------------------------------------------------*/
// //   Preset the coarsened-fine \<JU\> from coarse level data
// /**  \param[in]  a_CrJU \<JU\> on the coarse level
//  *//*-----------------------------------------------------------------*/

// void
// FourthOrderMappedFineInterp::presetCrFnLevJU(
//   const LevelData<FArrayBox>& a_CrLevJU)
// {
//   CH_assert(isGridDefined());
//   CH_assert(!m_mbcp.isDefined());  // Must not be used with multi-block
//   m_CrFnStatus |= CrFnStatusHaveJU;
//   a_CrLevJU.copyTo(m_CrFnLevJU, m_copier);
// }

/*--------------------------------------------------------------------*/
//   Select a stencil
/**  Selects a stencil and returns a number of useful parameters
 *   describing the stencil and associated transformations
 *   \param[in]  a_center
 *                      Location of coarse cell
 *   \param[in]  a_CrProblemDomain
 *                      The problem domain as specified on the coarse
 *                      level
 *   \param[out] a_iStencil
 *                      Index of the stencil to use
 *   \param[out] a_ctStoP
 *                      Transforms from stencil (S) to physical (P)
 *                      space
 *   \param[out] a_ctPtoS
 *                      Transforms from physical (P) to stencil (S)
 *                      space
 *   \param[out] a_offset
 *                      The offset defines the variation in location
 *                      (0,0) between S and P space.  We need the
 *                      offset since the location of (0,0) for a fine
 *                      box of cells in a coarse cells in stencil
 *                      space is a corner and is not the same in S and
 *                      P spaces (unlike the stencil box for which the
 *                      location of (0,0) is at the center of the box
 *                      and does not change in S and P spaces).
 *   \param[out] a_fineSmlEnd
 *                      When we transform an index from S to P, we
 *                      will need to add this to get the actual
 *                      IntVect in Z.
 *
 *   \note
 *   <ul>
 *     <li> For transforming the location of a fine cell, with (0,0)
 *          normalized to the location of the small end in the coarse
 *          cell, the operations are essentially
 *            P = Ts->p(S) + offset
 *          where
 *            offset[i] = 0            if Tp->s[i] is positive
 *                      = n_ref[i] - 1 if Tp->s[i] is negative
 *          and
 *            S = Tp->s(P - offset)
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::selectStencil(
  const IntVect        a_center,
  const ProblemDomain& a_CrProblemDomain,
  int&                 a_iStencil,
  CoordTransform&      a_ctStoP,
  CoordTransform&      a_ctPtoS,
  IntVect&             a_offset,
  IntVect&             a_fineSmlEnd) const
{
  // Most cells will be at the interior so use that as default
  int soi[SpaceDim] =
  {
    D_DECL6(0, 0, 0, 0, 0, 0)
  };
  int transform[SpaceDim] =
  {
    D_DECL6(1, 2, 3, 4, 5, 6)
  };
  Box offsetBox(a_center - m_stencilExtent*IntVect::Unit,
                a_center + m_stencilExtent*IntVect::Unit);
  if (!a_CrProblemDomain.domainBox().contains(offsetBox))
    {
      for (int iDir = 0; iDir != SpaceDim; ++iDir)
        {
          if (!a_CrProblemDomain.isPeriodic(iDir))
            {
              // We assume the box is big enough such that we can only
              // be at one boundary in a given direction.
              const int offsetDirHi =
                (a_center[iDir] + m_stencilExtent) -
                a_CrProblemDomain.domainBox().bigEnd(iDir);
              const int offsetDirLo =
                a_CrProblemDomain.domainBox().smallEnd(iDir) -
                (a_center[iDir] - m_stencilExtent);
              const int maxOffset = std::max(offsetDirHi, offsetDirLo);
              if (maxOffset > 0)
                {
                  soi[iDir] = maxOffset;
                  if (offsetDirLo > offsetDirHi)
                    {
                      transform[iDir] = -transform[iDir];
                    }
                }
            }
        }
      // This gives us both the stencil index and the coordinate
      // transformation
      sortEachSetDx(soi, transform);
    }
  a_iStencil = m_mapStencilIndex.find(offset2Key(soi))->second;

/*** Now we are in a transformed space.  Any offset o(ii, jj, kk) from      ***
 *** cell u(i, j, k) as prescribed by a stencil is now given by             ***
 *** ct.transform(o)                                                        ***/

  a_ctStoP.expand(transform);           // Transformation from stencil
                                        // coordinate sytem to physical
                                        // coordinate system
  a_ctPtoS = a_ctStoP;                  // Transformation from physical
  a_ctPtoS.reverse();                   // coordinate sytem to stencil
                                        // coordinate system

  // Define the offset
  a_offset = IntVect::Zero;
  for (int iDir = 0; iDir != SpaceDim; ++iDir)
    {
      if (a_ctPtoS.signAt(iDir) < 0)
        {
          a_offset[iDir] += m_refRatio[iDir] - 1;
        }
    }

  // When we transform an index from S to P, we will need to add this to
  // get the actual IntVect in Z.
  a_fineSmlEnd = a_center*m_refRatio + a_offset;
}

/*--------------------------------------------------------------------*/
//   Check for proper nesting of the fine level in the coarse level
/**  This is really expensive but does perform a thorough check.  If
 *   ghost cells are to be filled, then defineCFInterface should have
 *   already been called.  Otherwise, only sufficient nesting to
 *   interpolate all interior cells on the fine mesh is guaranteed.
 *   \param[in]  a_CrBoxes
 *                      Layout of the fine boxes.
 *   \param[in]  a_CrFnBoxes
 *                      Layout of the coarsened fine boxes.
 *   \param[in]  a_CrMetricsPtr
 *                      Used to get problem domain on the coarse mesh
 *//*-----------------------------------------------------------------*/

bool
FourthOrderMappedFineInterp::properNesting(
  const DisjointBoxLayout& a_CrBoxes,
  const DisjointBoxLayout& a_CrFnBoxes,
  LevelGridMetrics *const  a_CrMetricsPtr) const
{
  CH_TIME("FourthOrderMappedFineInterp::properNesting");
  IntVect totalRad = m_stencilExtent*IntVect::Unit;
  if (isCFInterfaceDefined())
    {
      totalRad += m_CrInterpRadVec;
    }
  if (m_willFillGhostsWithCrFnJU && !isMultiblock())
    {
      totalRad += IntVect::Unit;
    }

  DataIterator ditCrFn = a_CrFnBoxes.dataIterator();

  // A box which will determine whether a given box adjoins a periodic boundary
  Box crPeriodicTestBox;
  bool isPeriodicDomain = false;
  if (!isMultiblock())
    {
      const ProblemDomain& crProblemDomain = a_CrMetricsPtr->problemDomain();
      crPeriodicTestBox = crProblemDomain.domainBox();
      if (crProblemDomain.isPeriodic())
        {
          isPeriodicDomain = true;
          for (int iDir = 0; iDir < SpaceDim; ++iDir)
            {
              if (crProblemDomain.isPeriodic(iDir))
                {
                  crPeriodicTestBox.grow(iDir, -1);
                }
            }
        }
    }

  for (DataIterator ditCrFn = a_CrFnBoxes.dataIterator(); ditCrFn.ok();
       ++ditCrFn)
    {
      const Box crFnBox = a_CrFnBoxes[ditCrFn];
      //**FIXME currently, there is no testing for if there is sufficient
      //**      nesting across multiblock boundaries
      // const ProblemDomain& crBlockDomain =
      //   a_CrMetricsPtr->blockDomain(crFnBox, totalRad[0]);
      const ProblemDomain& crBlockDomain = a_CrMetricsPtr->blockDomain(crFnBox);
      // This is the number of cells in coarsened fine boxes that need to be
      // filled from real cells in the coarse grid
      Box crFnGhostBox = grow(crFnBox, totalRad);
      //**FIXME this really isn't sufficient.  You either need to be adjacent
      //**      to the boundary or m_CrInterpRadVec cells away.  In meshRefine,
      //**      the block factor generally ensures this.
      crFnGhostBox &= crBlockDomain;

      // Initialize the IVS as the whole crFnGhostBox
      IntVectSet crFnGhostIVS(crFnGhostBox);

      // Now we iterate over all the boxes in the coarse grid and subtract off
      // all the real cells.  All cells in the coarsened-fine box should be
      // removed
      LayoutIterator litCr = a_CrBoxes.layoutIterator();
      for (litCr.begin(); litCr.ok(); ++litCr)
        {
          const Box crRealBox = a_CrBoxes[litCr];

          crFnGhostIVS -= crRealBox;

          // We also need to remove any periodic images of real cells
          if (isPeriodicDomain
              && !crPeriodicTestBox.contains(crFnGhostBox)
              && !crPeriodicTestBox.contains(crRealBox))
            {
              ShiftIterator shiftIt = crBlockDomain.shiftIterator();
              IntVect shiftMult(crBlockDomain.domainBox().size());
              Box shiftedBox(crRealBox);
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect = shiftMult*shiftIt();
                  shiftedBox.shift(shiftVect);
                  crFnGhostIVS -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                }
            }
        }  // Second loop over all boxes

      // crFnGhostIVS should now be empty
      if (!crFnGhostIVS.isEmpty())
        {
          MayDay::Warning("Insufficient nesting detected.");
          return false;
        }
    }  // First loop over boxes on the processor
  return true;
}

/*--------------------------------------------------------------------*/
//   Sort sets of dimensions with the same mesh spacing
/**  Allows us to find a unique stencil index and also gives the
 *   transformation required to orient the physical coordinates in
 *   the coordinate system of this stencil
 *   \param[in]  a_a    The array to sort
 *   \param[out] a_a    Sorted from greatest to lowest for each dx
 *   \param[in]  a_b    An array to rearrange along with the sorting of a_a
 *   \param[out] a_b    Rearranged the same as a_a
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::sortEachSetDx(int *const a_a,
                                           int *const a_b) const
{
  int index[SpaceDim] =
  {
    D_DECL6(0, 1, 2, 3, 4, 5)
  };
  // Sort each unique dx separately
  bool mask[SpaceDim];
  for (int i = 0; i != m_numUniqueDx; ++i)
    {
      for (int j = 0; j != SpaceDim; ++j)
        {
          mask[j] = (m_uniqueDx[j] == i);
        }
      Sort::insertion(SpaceDim,
                      index,
                      Sort::CmpGreaterIndex<int>(a_a),
                      mask);
    }
  // Now rearrange both a and b to the index
  Sort::Move2Array<int, int> mover(a_a, a_b);
  Sort::rearrangeToIndex(SpaceDim, index, mover);
}

/*--------------------------------------------------------------------*/
//   Get the stencil key from a list of offsets
/**  The offsets must be sorted according to 'sortEachSetDx'
 *   \param[in]  a_offsets
 *                      Sorted offsets of the stencil from the
 *                      interior
 *   \return            A unique stencil index
 *//*-----------------------------------------------------------------*/

int
FourthOrderMappedFineInterp::offset2Key(const int *const a_offsets) const
{
  int mult = 1;
  int key = 0;
  for (int i = SpaceDim; i--;)
    {
      key += mult*a_offsets[i];
      mult *= m_stencilExtent + 1;
    }
  return key;
}

/*--------------------------------------------------------------------*/
//   Load <\delta\xi^\mathbf{p}> for all p
/**  We evaluate
 *   \f[
 *     \langle \delta\xi^{\mathbf{p}} \rangle
 *     =
 *     \prod_d \left(\frac
 *       {\displaystyle \int_{\delta\xi_d - \frac{h_d}{2}}
 *                          ^{\delta\xi_d + \frac{h_d}{2}}
 *         (\xi_d)^{p_d} \mathrm{d}\xi_d}
 *       {\displaystyle h_d}\right)
 *     =
 *     \prod_d \left(\frac
 *       {\displaystyle \left[
 *         \left(\delta\xi_d + \frac{h_d}{2}\right)^{p_d + 1} -
 *         \left(\delta\xi_d - \frac{h_d}{2}\right)^{p_d + 1}\right]}
 *       {\displaystyle h_d\left(p_d + 1\right)}\right)
 *   \f]
 *   for all values of \f$\|\mathbf{p}\|_1 \leq \mbox{m\_degree}\f$.
 *   \param[out] a_A    Vector of results
 *   \param[in]  a_dxi  Displacement of \f$\xi\f$ from center cell
 *   \param[in]  a_h    Mesh spacing
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::loadAvgXipA(
  Real*           a_A,
  const RealVect& a_dxi,
  const RealVect  a_h) const
{
  RealVect mhi;
  RealVect mlo;
  for (int iDir = 0; iDir != SpaceDim; ++iDir)
    {
      mhi[iDir] = a_dxi[iDir] + 0.5*a_h[iDir];
      mlo[iDir] = a_dxi[iDir] - 0.5*a_h[iDir];
    }
  // For direction '#'
  int D_DECL6(p0, p1, p2, p3, p4, p5);
                                      // p# is the power
  RealVect hi;                        // hi# is the integrand evaluated at the
                                      // upper bound (= mhi#^p#)
  RealVect lo;                        // lo# is the integrand evaluated at the
                                      // lower bound (= mlo#^p#)
  RealVect num;                       // num# is the factor of the numerator for
                                      // this and all lower dimensions
                                      // (= Prod[(hi# - lo#)(...)(hi0 - lo0)])
  RealVect den;                       // den# is the factor of the denominator
                                      // for this and all lower dimension
                                      // (= Prod[(p# + 1)(...)(p0 + 1)])
  Real hVol = a_h.product();          // cell volume
  // And the result is in num[SpaceDim-1]/den[SpaceDim-1]
  D_TERM6(for (p0 = 0, hi[0] = mhi[0], lo[0] = mlo[0], num[0] = (hi[0] - lo[0])       , den[0] = hVol  ; p0                          <= m_degree; ++p0, hi[0] *= mhi[0], lo[0] *= mlo[0], num[0] = (hi[0] - lo[0])       , den[0] = (p0 + 1)*hVol  ),
          for (p1 = 0, hi[1] = mhi[1], lo[1] = mlo[1], num[1] = (hi[1] - lo[1])*num[0], den[1] = den[0]; p0 + p1                     <= m_degree; ++p1, hi[1] *= mhi[1], lo[1] *= mlo[1], num[1] = (hi[1] - lo[1])*num[0], den[1] = (p1 + 1)*den[0]),
          for (p2 = 0, hi[2] = mhi[2], lo[2] = mlo[2], num[2] = (hi[2] - lo[2])*num[1], den[2] = den[1]; p0 + p1 + p2                <= m_degree; ++p2, hi[2] *= mhi[2], lo[2] *= mlo[2], num[2] = (hi[2] - lo[2])*num[1], den[2] = (p2 + 1)*den[1]),
          for (p3 = 0, hi[3] = mhi[3], lo[3] = mlo[3], num[3] = (hi[3] - lo[3])*num[2], den[3] = den[2]; p0 + p1 + p2 + p3           <= m_degree; ++p3, hi[3] *= mhi[3], lo[3] *= mlo[3], num[3] = (hi[3] - lo[3])*num[2], den[3] = (p3 + 1)*den[2]),
          for (p4 = 0, hi[4] = mhi[4], lo[4] = mlo[4], num[4] = (hi[4] - lo[4])*num[3], den[4] = den[3]; p0 + p1 + p2 + p3 + p4      <= m_degree; ++p4, hi[4] *= mhi[4], lo[4] *= mlo[4], num[4] = (hi[4] - lo[4])*num[3], den[4] = (p4 + 1)*den[3]),
          for (p5 = 0, hi[5] = mhi[5], lo[5] = mlo[5], num[5] = (hi[5] - lo[5])*num[4], den[5] = den[4]; p0 + p1 + p2 + p3 + p4 + p5 <= m_degree; ++p5, hi[5] *= mhi[5], lo[5] *= mlo[5], num[5] = (hi[5] - lo[5])*num[4], den[5] = (p5 + 1)*den[4]))
    {
      *a_A++ = num[SpaceDim-1]/den[SpaceDim-1];
    }
}

/*--------------------------------------------------------------------*/
//   Load <\delta\xi^\mathbf{p}> for all p
/**  We evaluate
 *   \f[
 *     \langle \delta\xi^{\mathbf{p}} \rangle
 *     =
 *     \prod_d \left(\frac
 *       {\displaystyle \int_{\delta\xi_d - \frac{h_d}{2}}
 *                          ^{\delta\xi_d + \frac{h_d}{2}}
 *         (\xi_d)^{p_d} \mathrm{d}\xi_d}
 *       {\displaystyle h_d}\right)
 *     =
 *     \prod_d \left(\frac
 *       {\displaystyle \left[
 *         \left(\delta\xi_d + \frac{h_d}{2}\right)^{p_d + 1} -
 *         \left(\delta\xi_d - \frac{h_d}{2}\right)^{p_d + 1}\right]}
 *       {\displaystyle h_d\left(p_d + 1\right)}\right)
 *   \f]
 *   for all values of \f$\|\mathbf{p}\|_1 \leq \mbox{m\_degree}\f$.
 *   \param[out] a_A    Vector of results
 *   \param[in]  a_dxi  Displacement of \f$\xi\f$ from center cell
 *   \param[in]  a_h    Mesh spacing
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::loadAvgXiKpA(
  Real*           a_A,
  const RealVect& a_dxi,
  const RealVect  a_h) const
{
  RealVect mhi;
  RealVect mlo;
  
  for (int iDir = 0; iDir != SpaceDim; ++iDir)
    {
      mhi[iDir] = a_dxi[iDir] + 0.5*a_h[iDir];
      mlo[iDir] = a_dxi[iDir] - 0.5*a_h[iDir];
    }
  // For direction '#'
  int D_DECL6(p0, p1, p2, p3, p4, p5);
                                      // p# is the power
  RealVect hi;                        // hi# is the integrand evaluated at the
                                      // upper bound (= mhi#^p#)
  RealVect lo;                        // lo# is the integrand evaluated at the
                                      // lower bound (= mlo#^p#)
  RealVect num;                       // num# is the factor of the numerator for
                                      // this and all lower dimensions
                                      // (= Prod[(hi# - lo#)(...)(hi0 - lo0)])
  RealVect den;                       // den# is the factor of the denominator
                                      // for this and all lower dimension
                                      // (= Prod[(p# + 1)(...)(p0 + 1)])
  Real hVol = a_h.product();          // cell volume
  // And the result is in num[SpaceDim-1]/den[SpaceDim-1]
  D_TERM6(for (p0 = 0, hi[0] = mhi[0], lo[0] = mlo[0], num[0] = (hi[0] - lo[0])       , den[0] = hVol  ; p0                          <= m_degree; ++p0, hi[0] *= mhi[0], lo[0] *= mlo[0], num[0] = (hi[0] - lo[0] - (p0 && p0%2 == 0)*std::pow(0.5, p0)*(mhi[0] - mlo[0]))       , den[0] = (p0 + 1)*hVol  ),
          for (p1 = 0, hi[1] = mhi[1], lo[1] = mlo[1], num[1] = (hi[1] - lo[1])*num[0], den[1] = den[0]; p0 + p1                     <= m_degree; ++p1, hi[1] *= mhi[1], lo[1] *= mlo[1], num[1] = (hi[1] - lo[1] - (p1 && p1%2 == 0)*std::pow(0.5, p1)*(mhi[1] - mlo[1]))*num[0], den[1] = (p1 + 1)*den[0]),
          for (p2 = 0, hi[2] = mhi[2], lo[2] = mlo[2], num[2] = (hi[2] - lo[2])*num[1], den[2] = den[1]; p0 + p1 + p2                <= m_degree; ++p2, hi[2] *= mhi[2], lo[2] *= mlo[2], num[2] = (hi[2] - lo[2] - (p2 && p2%2 == 0)*std::pow(0.5, p2)*(mhi[2] - mlo[2]))*num[1], den[2] = (p2 + 1)*den[1]),
          for (p3 = 0, hi[3] = mhi[3], lo[3] = mlo[3], num[3] = (hi[3] - lo[3])*num[2], den[3] = den[2]; p0 + p1 + p2 + p3           <= m_degree; ++p3, hi[3] *= mhi[3], lo[3] *= mlo[3], num[3] = (hi[3] - lo[3] - (p3 && p3%2 == 0)*std::pow(0.5, p3)*(mhi[3] - mlo[3]))*num[2], den[3] = (p3 + 1)*den[2]),
          for (p4 = 0, hi[4] = mhi[4], lo[4] = mlo[4], num[4] = (hi[4] - lo[4])*num[3], den[4] = den[3]; p0 + p1 + p2 + p3 + p4      <= m_degree; ++p4, hi[4] *= mhi[4], lo[4] *= mlo[4], num[4] = (hi[4] - lo[4] - (p4 && p4%2 == 0)*std::pow(0.5, p4)*(mhi[4] - mlo[4]))*num[3], den[4] = (p4 + 1)*den[3]),
          for (p5 = 0, hi[5] = mhi[5], lo[5] = mlo[5], num[5] = (hi[5] - lo[5])*num[4], den[5] = den[4]; p0 + p1 + p2 + p3 + p4 + p5 <= m_degree; ++p5, hi[5] *= mhi[5], lo[5] *= mlo[5], num[5] = (hi[5] - lo[5] - (p5 && p5%2 == 0)*std::pow(0.5, p5)*(mhi[5] - mlo[5]))*num[4], den[5] = (p5 + 1)*den[4]))
    {
      *a_A++ = num[SpaceDim-1]/den[SpaceDim-1];
    }
}

/*--------------------------------------------------------------------*/
//   Load gradients of <\delta\xi^\mathbf{p}> for all p
/**  These gradients are the same for the displacements determined
 *   from either loadAvgXip or loadAvgXiKp
 *   We evaluate
 *   \f[
 *     \langle \delta\frac{\mathrm{d}\xi^{\mathbf{p}}}
 *                        {\mathrm{d}\xi^{d}} \rangle
 *   \f]
 *   for all values of \f$\|\mathbf{p}\|_1 \leq \mbox{m\_degree}\f$.
 *   This is the same as loadAvgXipA except that a new coefficient is
 *   added in front (from the deriviate), p is reduced for
 *   direction d, and the constant K is removed.
 *   \param[out] a_A    Vector of results
 *   \param[in]  a_dxi  Displacement of \f$\xi\f$ from center cell
 *   \param[in]  a_h    Mesh spacing
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::loadAvgGradXip(
  Real*           a_gradXip,
  const RealVect& a_dxi,
  const RealVect  a_h) const
{
  RealVect mhi;
  RealVect mlo;
  
  for (int iDir = 0; iDir != SpaceDim; ++iDir)
    {
      mhi[iDir] = a_dxi[iDir] + 0.5*a_h[iDir];
      mlo[iDir] = a_dxi[iDir] - 0.5*a_h[iDir];
    }
  // For direction '#'
  IntVect p;
  D_TERM6(int& p0 = p[0];,
          int& p1 = p[1];,
          int& p2 = p[2];,
          int& p3 = p[3];,
          int& p4 = p[4];,
          int& p5 = p[5];)            // p# is the power
  RealVect hi;                        // hi# is the integrand evaluated at the
                                      // upper bound (= mhi#^p#)
  RealVect lo;                        // lo# is the integrand evaluated at the
                                      // lower bound (= mlo#^p#)
  RealVect num;                       // num# is the factor of the numerator for
                                      // this and all lower dimensions
                                      // (= Prod[(hi# - lo#)(...)(hi0 - lo0)])
  RealVect den;                       // den# is the factor of the denominator
                                      // for this and all lower dimension
                                      // (= Prod[(p# + 1)(...)(p0 + 1)])
  Real hVol = a_h.product();          // cell volume
  for (int dir = 0; dir != SpaceDim; ++dir)
    {

      // And the result is in p[dir]*num[SpaceDim-1]/den[SpaceDim-1]
      D_TERM6(for (p0 = 0, hi[0] = (dir != 0)*mhi[0] + (dir == 0), lo[0] = (dir != 0)*mlo[0] + (dir == 0), num[0] = (hi[0] - lo[0])       , den[0] = hVol  ; p0                          <= m_degree; ++p0, hi[0] *= mhi[0], lo[0] *= mlo[0], num[0] = (hi[0] - lo[0])       , den[0] = (p0 + (dir != 0))*hVol  ),
              for (p1 = 0, hi[1] = (dir != 1)*mhi[1] + (dir == 1), lo[1] = (dir != 1)*mlo[1] + (dir == 1), num[1] = (hi[1] - lo[1])*num[0], den[1] = den[0]; p0 + p1                     <= m_degree; ++p1, hi[1] *= mhi[1], lo[1] *= mlo[1], num[1] = (hi[1] - lo[1])*num[0], den[1] = (p1 + (dir != 1))*den[0]),
              for (p2 = 0, hi[2] = (dir != 2)*mhi[2] + (dir == 2), lo[2] = (dir != 2)*mlo[2] + (dir == 2), num[2] = (hi[2] - lo[2])*num[1], den[2] = den[1]; p0 + p1 + p2                <= m_degree; ++p2, hi[2] *= mhi[2], lo[2] *= mlo[2], num[2] = (hi[2] - lo[2])*num[1], den[2] = (p2 + (dir != 2))*den[1]),
              for (p3 = 0, hi[3] = (dir != 3)*mhi[3] + (dir == 3), lo[3] = (dir != 3)*mlo[3] + (dir == 3), num[3] = (hi[3] - lo[3])*num[2], den[3] = den[2]; p0 + p1 + p2 + p3           <= m_degree; ++p3, hi[3] *= mhi[3], lo[3] *= mlo[3], num[3] = (hi[3] - lo[3])*num[2], den[3] = (p3 + (dir != 3))*den[2]),
              for (p4 = 0, hi[4] = (dir != 4)*mhi[4] + (dir == 4), lo[4] = (dir != 4)*mlo[4] + (dir == 4), num[4] = (hi[4] - lo[4])*num[3], den[4] = den[3]; p0 + p1 + p2 + p3 + p4      <= m_degree; ++p4, hi[4] *= mhi[4], lo[4] *= mlo[4], num[4] = (hi[4] - lo[4])*num[3], den[4] = (p4 + (dir != 4))*den[3]),
              for (p5 = 0, hi[5] = (dir != 5)*mhi[5] + (dir == 5), lo[5] = (dir != 5)*mlo[5] + (dir == 5), num[5] = (hi[5] - lo[5])*num[4], den[5] = den[4]; p0 + p1 + p2 + p3 + p4 + p5 <= m_degree; ++p5, hi[5] *= mhi[5], lo[5] *= mlo[5], num[5] = (hi[5] - lo[5])*num[4], den[5] = (p5 + (dir != 5))*den[4]))
        {
          // p is from derivative, a_h is to make the derivative "undivided"
          *a_gradXip++ = p[dir]*a_h[dir]*num[SpaceDim-1]/den[SpaceDim-1];
        }
    }
}

#include "NamespaceFooter.H"
