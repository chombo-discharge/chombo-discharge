/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMultigridInterpolator.H
  @brief  Declaration of a factory class for interpolating ghost cells in a multigrid context. 
  @author Robert Marskar
  @todo   Dropping in lo/hi CFIVS so we can move performance penalty to this class rather than EBHelmholtzOp
*/

#ifndef CD_EBMultigridInterpolator_H
#define CD_EBMultigridInterpolator_H

// Std includes
#include <map>

// Chombo includes
#include <EBLevelGrid.H>
#include <QuadCFInterp.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Multigrid interpolator class. 
  @note Currently just a shell of EBQuadCFInterp (so we have something to work with). 
*/
class EBMultigridInterpolator : public QuadCFInterp {
public:

  /*!
    @brief Specification of cell-center or cell-centroid type of interpolation
  */
  using CellLocation = Location::Cell;

  /*!
    @brief Disallowed constructor. 
  */
  EBMultigridInterpolator() = delete;

  /*!
    @brief Disallowed copy constructor.
  */
  EBMultigridInterpolator(const EBMultigridInterpolator& a_other) = delete;

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
			  const IntVect&     a_minGhost,
			  const int          a_refRat,
			  const int          a_nVar,
			  const int          a_ghostCF,
			  const int          a_order,
			  const int          a_weighting = 0);
			  

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
    @param[in] a_phiFine   Fine-level data
    @param[in] a_variables Variables to be interpolated
    @param[in] a_dit       Data index
  */
  void coarseFineInterpH(EBCellFAB& a_phiFine, const Interval a_variables, const DataIndex& a_dit);

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
    @brief Component number
  */
  static constexpr int m_comp = 0;

  /*!
    @brief How to interpret cell location
  */
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
    @brief Minimum number of ghost cells in input
  */
  IntVect m_minGhost;

  /*!
    @brief Refinement factor
  */
  int m_refRat;

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
    @brief Regular coarse-fine interface boxes
  */
  LayoutData<std::map<std::pair<int, Side::LoHiSide>, Box> > m_cfivs;

  /*!
    @brief Defines the CFIVS
  */
  void defineCFIVS();

  /*!
    @brief Define ghost cells to be filled in range m_ghostCF from each cut-cell
  */
  void defineGhostRegions();

  /*!
    @brief Define grids
  */
  void defineGrids(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar);

  /*!
    @brief Define data holders. 
  */
  void defineBuffers();

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
