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
*/

#ifndef CD_EBMultigridInterpolator_H
#define CD_EBMultigridInterpolator_H

// Chombo includes
#include <EBLevelGrid.H>
#include <EBQuadCFInterp.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Multigrid interpolator class. 
  @note Currently just a shell of EBQuadCFInterp (so we have something to work with). 
*/
class EBMultigridInterpolator : public EBQuadCFInterp {
public:

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
    @param[in] a_ghostCells Ghost cells to interpolate
  */
  EBMultigridInterpolator(const EBLevelGrid&            a_eblgFine,
			  const EBLevelGrid&            a_eblgCoar,
			  const int                     a_nRef,
			  const int                     a_nVar,
			  const LayoutData<IntVectSet>& a_ghostCells);


  /*!
    @brief Dtor. Releases stencils. 
  */
  ~EBMultigridInterpolator();

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

protected:

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Coarse grids
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Zero data on coarse grid

  */
  LevelData<EBCellFAB> m_zeroCoar;
};

#include <CD_NamespaceFooter.H>

#endif
