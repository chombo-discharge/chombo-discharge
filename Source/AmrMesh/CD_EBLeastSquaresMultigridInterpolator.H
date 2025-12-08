/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBLeastSquaresMultigridInterpolator.H
  @brief  Declaration of a class that can interpolate more ghost cells near the coarse-fine boundary near the EB. 
  @author Robert Marskar
*/

#ifndef CD_EBLeastSquaresMultigridInterpolator_H
#define CD_EBLeastSquaresMultigridInterpolator_H

// Std includes
#include <map>

// Chombo includes
#include <EBLevelGrid.H>
#include <QuadCFInterp.H>
#include <AggStencil.H>

// Our includes
#include <CD_Location.H>
#include <CD_EBMultigridInterpolator.H>
#include <CD_CoarseInterpQuadCF.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Multigrid interpolator class. 
  @details This interpolator class interpolates with classical multigrid-type interpolation where the interpolation stencils only reach into valid cells, 
  ignoring the region underneath the fine level. Near regular refinement boundaries we use the conventional approach where we interpolate coarse data onto the
  line of cells that runs between the fine-grid cells and the ghost cell. This value is used together with the valid data on the fine-grid cells to interpolate
  the ghost cell to O(h^2). Near EBs this is much more complex since both the fine and coarse regions may be missing cells that make that interpolation possible. 
  Worse, EB stencils may be much bigger than the 5/7-point Laplacian, requiring us to fill more ghost cells across the coarse-fine interface near the EB. We achieve
  this by switching to least-squares based interpolation here, again using only the valid fine and coarse-level cells. The way of achieving this consists of building
  dual-level stencils for each ghost cell; these stencils represent the part of the interpolation stencil that reach into the fine and coarse grids. For efficiency 
  reasons, we copy the coarse grid data onto data layouts that may have more ghost cells than your regular data. Since the least squares stencils require a 
  specified number of cells to reach specified order, this permits us to use fine-grid stencils that only reach into valid data without dropping order. This optimization
  is useful for pruning extra communications due homogeneous multigrid interpolation. This means that our interpolation, which is in the form
  
     phiGhost = sum(vofCoar) w(vofCoar)*phi(vofCoar) + sum(vofFine)[w(vofFine)*phi(vofFine)]
  
  will have a coarse and and a fine grid stencil, but the fine stencil does not reach out of the ghosted fine-grid patch. The fine-grid part of the stencil is thus
  allowed to use the ghost cells in the patch that are filled from neighboring patches. 
*/
class EBLeastSquaresMultigridInterpolator : public EBMultigridInterpolator
{
public:
  /*!
    @brief Specification of cell-center or cell-centroid type of interpolation
  */
  using CellLocation = Location::Cell;

  /*!
    @brief Disallowed constructor. 
  */
  EBLeastSquaresMultigridInterpolator() = delete;

  /*!
    @brief Disallowed copy constructor.
  */
  EBLeastSquaresMultigridInterpolator(const EBLeastSquaresMultigridInterpolator& a_other) = delete;

  /*!
    @brief Full constructor. 
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgCoFi     Coarsened fine grids
    @param[in] a_eblgCoar     Coarse grids
    @param[in] a_dataLocation Interpration of data centering. Either on the cell center or the cell centroid. 
    @param[in] a_ghostVector  (Minimum) number of ghost cells in input/output data.
    @param[in] a_refRat       Refinement ratio
    @param[in] a_ghostCF      Width of the ghost region to be filled. Only relevant near the EBCF. 
    @param[in] a_order        Interpolation order
    @param[in] a_weighting    Weighting factor least squares.
    @note This class does not need the full grid resolution since we are only interested in least squares reconstruction of the 0th order term, so the
    other unknowns in the Taylor expansion that is computed are scaled by the grid resolution and then promptly discarded.
  */
  EBLeastSquaresMultigridInterpolator(const EBLevelGrid& a_eblgFine,
                                      const EBLevelGrid& a_eblgCoFi,
                                      const EBLevelGrid& a_eblgCoar,
                                      const CellLocation a_dataLocation,
                                      const IntVect&     a_ghostVector,
                                      const int          a_refRat,
                                      const int          a_ghostCF,
                                      const int          a_order,
                                      const int          a_weighting = 0) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBLeastSquaresMultigridInterpolator() noexcept;

  /*!
    @brief Return number of ghost cells across CF interface.
  */
  virtual int
  getGhostCF() const noexcept override;

  /*!
    @brief Get the interpolation stencil near the refinement boundary.
    @details This version can assume that all cells are regular, and that returned stencils may later be replaced
    by getInterpolationStencilEB. This function is not used in the interpolation calls, but is useful if the user
    wants to extract the stencils for filling matrices directly (e.g., if using PETSc). 
    @param[in] a_fineGhost Ghost cell index on the fine grid
    @param[in] a_din Grid index
    @param[in] a_dir Coordinate direction
    @param[in] a_side Lo or high side of the grid patch
  */
  virtual std::pair<DerivStencil, DerivStencil>
  getInterpolationStencilRegular(const IntVect&       a_fineGhost,
                                 const DataIndex&     a_din,
                                 const int            a_dir,
                                 const Side::LoHiSide a_side) const noexcept override;

  /*!
    @brief Get an explicit stencil for a ghost cell near the refinement boundary.
    @details This should return the fine + coarse stencils required for interpolating to the fine ghost cell. This
    function is not used in the interpolation calls, but is useful if the user wants to extract the stencils for
    filling matrices directly (e.g., if using PETSc). 
    @param[in] a_fineGhost Ghost cell on the fine side of the refinement boundary.
    @param[in] a_din Grid index.
    @return This first/second entry in the returned pair should be the fine/coarse stencil.    
  */
  virtual std::pair<VoFStencil, VoFStencil>
  getInterpolationStencilEB(const VolIndex& a_fineGhost, const DataIndex& a_din) const noexcept override;

  /*!
    @brief Do inhomogeneous interpolation
    @param[inout] a_phiFine   Fine phi
    @param[in]    a_phiCoar   Coarse phi
    @param[in]    a_variables Variables to interpolate
  */
  virtual void
  coarseFineInterp(LevelData<EBCellFAB>&       a_phiFine,
                   const LevelData<EBCellFAB>& a_phiCoar,
                   const Interval              a_variables) const noexcept override;

  /*!
    @brief Do homogeneous interpolation
    @param[in] a_phiFine   Fine phi
    @param[in] a_variables Variables to interpolate
  */
  virtual void
  coarseFineInterpH(LevelData<EBCellFAB>& a_phiFine, const Interval a_variables) const noexcept override;

  /*!
    @brief Do homogeneous coarse-fine interpolation on a grid patch
    @param[in] a_phiFine   Fine-level data
    @param[in] a_variables Variables to be interpolated
    @param[in] a_din       Data index
    @note This is the kernel version. The constructor guarantees that the fine-grid stencil will not out of the ghost region defined by a_ghostVector. However, this is
    the direct kernel version which is public because multigrid operators might need it for performance reasons. This routine does NOT warn against input data that do not
    contain enough ghost cells. 
  */
  virtual void
  coarseFineInterpH(EBCellFAB& a_phiFine, const Interval a_variables, const DataIndex& a_din) const noexcept override;

protected:
  /*!
    @brief Component that interpolation stencil is stored at
  */
  static constexpr int m_stenComp = 0;

  /*!
    @brief Number of variables in stencil data layouts. 
  */
  static constexpr int m_numStenComp = 1;

  /*!
    @brief Component number. Used for setting interval when data is copied into EBLeastSquaresMultigridInterpolator's buffers. 
  */
  static constexpr int m_comp = 0;

  /*!
    @brief How to interpret data locations. 
  */
  CellLocation m_dataLocation;

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Coarse grids
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Coarsened fine grids.
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Copier for making copying into m_grownCoarData go faster
  */
  Copier m_copier;

  /*!
    @brief Iterator over ghost cells
  */
  mutable LayoutData<VoFIterator> m_ghostIterFine;

  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Minimum number of ghost cells in input data. 
  */
  IntVect m_ghostVectorFine;

  /*!
    @brief Number of ghost cells on the coarsened-fine grids. Needed because the fine grids need some of
    the coarse grid cells.
  */
  IntVect m_ghostVectorCoFi;

  /*!
    @brief Refinement factor between fine and coarse level. 
  */
  int m_refRat;

  /*!
    @brief Number of ghost cells to fill across coarse-fine interface. 
    @note Only applies the the cut-cells. 
  */
  int m_ghostCF;

  /*!
    @brief Interpolation order. 
  */
  int m_order;

  /*!
    @brief Least squares weighting factor
  */
  int m_weight;

  /*!
    @brief Ghost cells that will be interpolated near the cut-cells. 
    @note This does not include the ghost cells on the refinement boundary away from the cut-cells. 
  */
  LayoutData<IntVectSet> m_ghostCells;

  /*!
    @brief Fine-grid interpolation stencils.
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_fineStencils;

  /*!
    @brief Coarse-grid interpolation stencils.
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_coarStencils;

  /*!
    @brief Fine grid interpolation stencils. 
    @details This is used for performance optimization. Note that m_aggFineStencils
    will reach into whatever data is being interpolated, so the input data must have
    the same number of ghost cells as specified in the constructor. 
  */
  LayoutData<RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB>>> m_aggFineStencils;

  /*!
    @brief Coarse-grid interpolation stencils. 
    @details This is used for performance optimization. Note that m_aggCoarStencils
    will only reach into m_grownCoarData, which gives well-defined cache offsets. 
  */
  LayoutData<RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB>>> m_aggCoarStencils;

  /*!
    @brief Regular coarse-fine interface boxes. 
    @note This is a map to the "regular" ghost cells that will be interpolated using a standard
    multigrid interpolation stencil (away from the EB). 
  */
  LayoutData<std::map<std::pair<int, Side::LoHiSide>, Box>> m_cfivs;

  /*!
    @brief For doing the coarse-side interpolation near the AMR interface
  */
  LayoutData<CoarseInterpQuadCF> m_loCoarseInterpCF[SpaceDim];

  /*!
    @brief For doing the coarse-side interpolation near the AMR interface
  */
  LayoutData<CoarseInterpQuadCF> m_hiCoarseInterpCF[SpaceDim];

  /*!
    @brief Define ghost cells to be filled in range m_ghostCF from each cut-cell. 
  */
  virtual void
  defineGhostRegions() noexcept;

  /*!
    @brief Define buffer data holder. 
  */
  virtual void
  defineBuffers() noexcept;

  /*!
    @brief Define the coarse-sided interpolation stencils for doing the orthogonal derivatives. 
  */
  virtual void
  defineCoarseInterp() noexcept;

  /*!
    @brief Define stencils for interpolation across the EBCF refinement boundary. 
  */
  virtual void
  defineStencilsEBCF() noexcept;

  /*!
    @brief Build the AggStencil objects for faster stencil aggregation
  */
  virtual void
  makeAggStencils() noexcept;

  /*!
    @brief Build an interpolation stencil
    @param[out] a_stencilFine    Fine-grid stencil
    @param[out] a_stencilCoar    Coarse-grid stencil
    @param[in]  a_dataLocation   Data location
    @param[in]  a_ghostVofFine   Ghost cell which is to be filled on the coarse grid. 
    @param[in]  a_ghostVofCoar   Vof corresponding to a_ghostVofFine on the coarse grid
    @param[in]  a_ebisBoxFine    Fine grid EBISBox
    @param[in]  a_ebisBxoCoar    Coarse grid EBISBox
    @param[in]  a_validFineCells Fine-grid region which we are allowed to use in the interpolation. 
    @param[in]  a_validCoarCells Coarse-grid region which we are allowed to use in the interpolation. 
    @param[in]  a_dxFine         Grid resolution on the fine level
    @param[in]  a_dxCoar         Grid resolution on the coarse level
    @param[in]  a_order          Interpolation order
    @param[in]  a_weight         Equation weighting for least squares reconstruction. 
    @return Returns true if a stencil was found, in which case it is defined in a_stencilFine and a_stencilCoar. 
  */
  virtual bool
  getStencil(VoFStencil&            a_stencilFine,
             VoFStencil&            a_stencilCoar,
             const CellLocation&    a_dataLocation,
             const VolIndex&        a_ghostVofFine,
             const VolIndex&        a_ghostVofCoar,
             const EBISBox&         a_ebisboxFine,
             const EBISBox&         a_ebisboxCoar,
             const DenseIntVectSet& a_validFineCells,
             const DenseIntVectSet& a_validCoarCells,
             const Real&            a_dxFine,
             const Real&            a_dxCoar,
             const int&             a_order,
             const int&             a_weight) const noexcept;

  /*!
    @brief Do regular coarse-fine interpolation as if the EB is not there.
    @param[inout] a_finePhi Fine-level data.
    @param[in] a_coarPhi Coarse-grid data.
    @param[in] a_fineVar Fine variable
    @param[in] a_coarVar Coarse variable
  */
  virtual void
  regularCoarseFineInterp(LevelData<EBCellFAB>&       a_finePhi,
                          const LevelData<EBCellFAB>& a_coarPhi,
                          const int                   a_fineVar,
                          const int                   a_coarVar) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
