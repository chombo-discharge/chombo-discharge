/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMultigridInterpolator.H
  @brief  Declaration of a factory class for interpolating ghost cells in a multigrid context. 
  @author Robert Marskar
  @todo  m_zeroCoar exists because I want EBQuadCFInterp to run under the hood right now. Later, we will use the stencils directly and there
  will be no need for m_zeroCoar
  @todo   Dropping in lo/hi CFIVS so we can move performance penalty to this class rather than EBHelmholtzOp
*/

#ifndef CD_EBMultigridInterpolator_H
#define CD_EBMultigridInterpolator_H

// Chombo includes
#include <EBLevelGrid.H>
#include <EBQuadCFInterp.H>
#include <CFIVS.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Multigrid interpolator class. 
  @note Currently just a shell of EBQuadCFInterp (so we have something to work with). 
*/
class EBMultigridInterpolator : public EBQuadCFInterp {
public:

  /*!
    @brief Specification of cell-center or cell-centroid type of interpolation
  */
  using CellLocation = Location::Cell;

  /*!
    @brief Default constructor. Must subsequently call define
  */
  EBMultigridInterpolator();

  /*!
    @brief Disallowed copy constructor.
  */
  EBMultigridInterpolator(const EBMultigridInterpolator& a_other) = delete;

  /*!
    @brief Full constructor
    @param[in] a_eblgFine   Fine grids
    @param[in] a_eblgCoar   Coarse grids
    @param[in] a_refRat     Refinement ratio
    @param[in] a_nVar       Number of variables
    @param[in] a_ghostCells Number of ghost cells to interpolate. 
  */
  EBMultigridInterpolator(const EBLevelGrid& a_eblgFine,
			  const EBLevelGrid& a_eblgCoar,
			  const int          a_nRef,
			  const int          a_nVar,
			  const int          a_ghostCF);

  /*!
    @brief Full constructor
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgCoar     Coarse grids
    @param[in] a_cellLocation Interpration of data centering. Either on the cell center or the cell centroid. 
    @param[in] a_refRat       Refinement ratio
    @param[in] a_nVar         Number of variables
    @param[in] a_ghostCF      Width of the ghost region to be filled
    @param[in] a_order        Interpolation order
    @param[in] a_order        Interpolation order
    @param[in] a_weighting    Weighting factor least squares. 
  */
  EBMultigridInterpolator(const EBLevelGrid& a_eblgFine,
			  const EBLevelGrid& a_eblgCoar,
			  const CellLocation a_cellLocation,
			  const int          a_refRat,
			  const int          a_nVar,
			  const int          a_ghostCF,
			  const int          a_order,
			  const int          a_weighting);
			  

  /*!
    @brief Dtor. Releases stencils. 
  */
  ~EBMultigridInterpolator();

  /*!
    @brief Return number of ghost cells across CF interface.
  */
  int getGhostCF() const;

  /*!
    @brief Do inhomogeneous interpolation
    @param[inout] a_phiFine   Fine phi
    @param[in]    a_phiCoar   Coarse phi
    @param[in]    a_variables Variables to interpolate
  */
  void coarseFineInterp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar, const Interval a_variables);

  /*!
    @brief Do homogeneous interpolation
    @param[in] a_phiFine   Fine phi
    @param[in] a_variables Variables to interpolate
  */
  void coarseFineInterpH(LevelData<EBCellFAB>& a_phiFine, const Interval a_variables);

  /*!
    @brief Do coarse-fine interpolation on a grid patch
  */
  void coarseFineInterpH(EBCellFAB& a_phiFine, const Interval a_variables, const DataIndex& a_dit);

  /*!
    @brief Inhomogeneous interpolation
    @param[inout] a_phiFine   Fine-level data to be interpolated. 
    @param[in]    a_phiCoar   Coarse-level data. 
    @param[in]    a_variables Variables to be interpolated. 
  */
  void interp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar, const Interval a_variables) const;

  /*!
    @brief Homogeneous interpolation
    @param[in] a_phiFine Fine-level data. We interpolate the ghost cells assuming that coarse-level data is zero everywhere. 
  */
  void interpH(LevelData<EBCellFAB>& a_phiFine, const Interval a_variables) const;

protected:
  
  /*!
    @brief Component that interpolation stencil is stored at
  */
  static constexpr int m_stenComp = 0;

  /*!
    @brief Number of stencil components
  */
  static constexpr int m_nStenComp = 1;
  
  /*!
    @brief Defines the CFIVS
  */
  void defineCFIVS();

    /*!
    @brief Coarse-fine interface. 
  */
  LayoutData<CFIVS> m_loCFIVS[SpaceDim];

  /*!
    @brief Coarse-fine interface. 
  */
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

// ==================================================
// NEW STUFF BELOW HERE
// ==================================================

  CellLocation m_cellLocation;

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Coarse grids
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Refined coarse grids
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Grown box layout on the fine level
  */
  BoxLayout m_fineBoxes;

  /*!
    @brief Grown box layout on the coarse level
  */
  BoxLayout m_coarBoxes;

  /*!
    @brief Fine data. Contains a lot of ghost storage
  */
  mutable BoxLayoutData<EBCellFAB> m_fineData;

  /*!
    @brief Coarsened data
  */
  mutable BoxLayoutData<EBCellFAB> m_coarData;

  /*!
    @brief Number of ghost cells to fill across coarse-fine interface. 
  */
  int m_ghostCF;

  /*!
    @brief Order of the interpolation
  */
  int m_order;

  /*!
    @brief Least squares weighting factor
  */
  int m_weight;

  /*!
    @brief Ghost cells that may require interpolation. 
  */
  LayoutData<IntVectSet> m_ghostCells;

  /*!
    @brief Fine stencils
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_fineStencils;

  /*!
    @brief Fine stencils
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_coarStencils;

  /*!
    @brief Define ghost cells to be filled in range m_ghostCF from each cut-cell
  */
  void defineGhosts();

  /*!
    @brief Define grids
  */
  void defineGrids(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar);

  /*!
    @brief Define data holders. 
  */
  void defineData();

  /*!
    @brief Define stencils
  */
  void defineStencils();

  /*!
    @brief Build an interpolation stencil
  */
  bool getStencil(VoFStencil&         a_stencilFine,
		  VoFStencil&         a_stencilCoar,
		  const CellLocation& a_cellLocation,
		  const VolIndex&     a_ghostVof,
		  const EBISBox&      a_ebisboxFine,
		  const EBISBox&      a_ebisboxCoar,
		  const IntVectSet&   a_validFineCells,
		  const IntVectSet&   a_validCoarCells,
		  const Real&         a_dxFine,
		  const Real&         a_dxCoar,
		  const int&          a_order,
		  const int&          a_weight);
};

#include <CD_NamespaceFooter.H>

#endif
