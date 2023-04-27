/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBFineInterp.H
  @brief  Declaration of an aggregated class for regrid operations.
  @author Robert Marskar
*/

#ifndef CD_EBFineInterp_H
#define CD_EBFineInterp_H

// Chombo includes
#include <EBLevelGrid.H>
#include <EBPWLFineInterp.H>

// Our includes
#include <CD_Regrid.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for interpolating data to fine grids. Can use constant interpolation or include limiters.
*/
class EBFineInterp : public EBPWLFineInterp
{
public:
  /*!
    @brief Default constructor. Must subsequently call the define method. 
  */
  EBFineInterp() noexcept;

  /*!
    @brief Constructor. Constructs a valid object which can be used for regridding.
    @param[in] a_eblgFine The fine grid. 
    @param[in] a_eblgCoFi Coarsened fine grid. 
    @param[in] a_eblgCoar The coarse grid. 
    @param[in] a_refRat   Refinement ratio between fine and coarse. 
    @param[in] a_nComp    Number of components. 
    @param[in] a_ebisPtr  EB index space. 
  */
  EBFineInterp(const EBLevelGrid&        a_eblgFine,
               const EBLevelGrid&        a_eblgCoFi,
               const EBLevelGrid&        a_eblgCoar,
               const int&                a_refRat,
               const int&                a_nComp,
               const EBIndexSpace* const a_ebisPtr) noexcept;

  /*!
    @brief Disallowed copy constructor
  */
  EBFineInterp(const EBFineInterp&) = delete;

  /*!
    @brief Disallowed move constructor
  */
  EBFineInterp(const EBFineInterp&&) = delete;

  /*!
    @brief Destructor
  */
  virtual ~EBFineInterp() noexcept;

  /*!
    @brief Disallowed copy operator
  */
  EBFineInterp&
  operator=(const EBFineInterp&) = delete;

  /*!
    @brief Disallowed move operator
  */
  EBFineInterp&
  operator=(const EBFineInterp&&) = delete;

  /*!
    @brief Define operation. Puts object in usable state.
    @param[in] a_eblgFine The fine grid. 
    @param[in] a_eblgCoFi Coarsened fine grid. 
    @param[in] a_eblgCoar The coarse grid. 
    @param[in] a_refRat   Refinement ratio between fine and coarse. 
    @param[in] a_nComp    Number of components. 
    @param[in] a_ebisPtr  EB index space. 
  */
  virtual void
  define(const EBLevelGrid&        a_eblgFine,
         const EBLevelGrid&        a_eblgCoFi,
         const EBLevelGrid&        a_eblgCoar,
         const int&                a_refRat,
         const int&                a_nComp,
         const EBIndexSpace* const a_ebisPtr) noexcept;

  /*!
    @brief Do a piecewise constant interpolation to the new grids.
    @param[out] a_fineData  Piecewise constant interpolation of the coarse data.
    @param[in]  a_coarData  Coarse data
    @param[in]  a_variables Variables to interpolate. 
  */
  virtual void
  regridConservativeNoSlopes(LevelData<EBCellFAB>&       a_fineData,
                             const LevelData<EBCellFAB>& a_coarData,
                             const Interval&             a_variables) const noexcept;

  /*!
    @brief Use min-mod slopes when interpolating to the new grids. 
    @param[out] a_fineData  Piecewise constant interpolation of the coarse data.
    @param[in]  a_coarData  Coarse data
    @param[in]  a_variables Variables to interpolate. 
  */
  virtual void
  regridMinMod(LevelData<EBCellFAB>& a_fineData, const LevelData<EBCellFAB>& a_coarData, const Interval& a_variables);

  /*!
    @brief Regrid irregular data using conservative interpolation
    @param[out] a_fineData  Piecewise interpolation of coarse data.
    @param[in]  a_coarData  Coarse data
    @param[in]  a_variables Variables to interpolate. 
  */
  virtual void
  regridConservative(LevelData<BaseIVFAB<Real>>&       a_fineData,
                     const LevelData<BaseIVFAB<Real>>& a_coarData,
                     const Interval&                   a_variables);

  /*!
    @brief Regrid irregular data using arithmetic interpolation
    @param[out] a_fineData  Piecewise interpolation of coarse data.
    @param[in]  a_coarData  Coarse data
    @param[in]  a_variables Variables to interpolate. 
  */
  virtual void
  regridArithmetic(LevelData<BaseIVFAB<Real>>&       a_fineData,
                   const LevelData<BaseIVFAB<Real>>& a_coarData,
                   const Interval&                   a_variables);

protected:
  /*!
    @brief Fine grid.
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Coarse grid
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Coarsened fine grid (i.e., has same ProblemDomain as m_eblgCoar);
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Data holder for holding data on the coarsened fine grids. 
  */
  LevelData<BaseIVFAB<Real>> m_irregCoFi;

  /*!
    @brief Iterator for valid cut-cells. 
    @details This iterator runs over the valid fine cells.
  */
  LayoutData<VoFIterator> m_fineVoFs;

  /*!
    @brief Weights in cut-cells when regridding conservatively
    @details Defined over the fine layout. 
  */
  LevelData<BaseIVFAB<Real>> m_conservativeWeights;

  /*!
    @brief Define weights for cut-cell conservative regridding
  */
  virtual void
  defineWeights() noexcept;

  /*!
    @brief Piecewise constant interpolation of an EBCellFAB 
    @param[out] a_fineData Piecewise constant interpolation of the coarse data.
    @param[in]  a_coarData Coarse data
    @param[in]  a_dit Grid index. 
    @param[in]  a_fineVar Fine-grid variable
    @param[in]  a_coarVar Coarse-grid variable
  */
  virtual void
  regridConservativeNoSlopes(EBCellFAB&       a_fineData,
                             const EBCellFAB& a_coarData,
                             const DataIndex& a_dit,
                             const int&       a_fineVar,
                             const int&       a_coarVar) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
