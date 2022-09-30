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
  EBFineInterp();

  /*!
    @brief Constructor. Constructs a valid object which can be used for regridding.
    @param[in] a_eblgFine The fine grid. 
    @param[in] a_eblgCoar The coarse grid. 
    @param[in] a_refRat   Refinement ratio between fine and coarse. 
    @param[in] a_nComp    Number of components. 
    @param[in] a_ebisPtr  EB index space. 
  */
  EBFineInterp(const EBLevelGrid&        a_eblgFine,
               const EBLevelGrid&        a_eblgCoar,
               const int&                a_refRat,
               const int&                a_nComp,
               const EBIndexSpace* const a_ebisPtr);

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
  virtual ~EBFineInterp();

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
    @param[in] a_eblgCoar The coarse grid. 
    @param[in] a_refRat   Refinement ratio between fine and coarse. 
    @param[in] a_nComp    Number of components. 
    @param[in] a_ebisPtr  EB index space. 
  */
  void
  define(const EBLevelGrid&        a_eblgFine,
         const EBLevelGrid&        a_eblgCoar,
         const int&                a_refRat,
         const int&                a_nComp,
         const EBIndexSpace* const a_ebisPtr);

  /*!
    @brief Do a piecewise constant interpolation to the new grids.
    @param[out] a_fineData  Piecewise constant interpolation of the coarse data.
    @param[in]  a_coarData  Coarse data
    @param[in]  a_variables Variables to interpolate. 
  */
  void
  regridNoSlopes(LevelData<EBCellFAB>& a_fineData, const LevelData<EBCellFAB>& a_coarData, const Interval& a_variables);

  /*!
    @brief Use min-mod slopes when interpolating to the new grids. 
    @param[out] a_fineData  Piecewise constant interpolation of the coarse data.
    @param[in]  a_coarData  Coarse data
    @param[in]  a_variables Variables to interpolate. 
  */
  void
  regridMinMod(LevelData<EBCellFAB>& a_fineData, const LevelData<EBCellFAB>& a_coarData, const Interval& a_variables);

  /*!
    @brief Regrid irregular data using conservative interpolation
    @param[out] a_fineData  Piecewise interpolation of coarse data.
    @param[in]  a_coarData  Coarse data
    @param[in]  a_variables Variables to interpolate. 
  */
  void
  regridConservative(LevelData<BaseIVFAB<Real>>&       a_fineData,
                     const LevelData<BaseIVFAB<Real>>& a_coarData,
                     const Interval&                   a_variables);

  /*!
    @brief Regrid irregular data using arithmetic interpolation
    @param[out] a_fineData  Piecewise interpolation of coarse data.
    @param[in]  a_coarData  Coarse data
    @param[in]  a_variables Variables to interpolate. 
  */
  void
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
    @brief Data holder for holding data on the coarsened fine grids. 
  */
  LevelData<BaseIVFAB<Real>> m_coarsenedFineIrregData;

  /*!
    @brief Iterator for valid cut-cells. 
    @details This iterator runs over the coarse cells on the coarsened fine layout. 
  */
  LayoutData<VoFIterator> m_irregVoFs;

  /*!
    @brief Piecewise constant interpolation of an EBCellFAB 
    @param[out] a_fineData  Piecewise constant interpolation of the coarse data.
    @param[in]  a_coarData  Coarse data
    @param[in]  a_dit       Grid index. 
    @param[in]  a_variables Variables to interpolate. 
  */
  void
  regridNoSlopes(EBCellFAB&       a_fineData,
                 const EBCellFAB& a_coarData,
                 const DataIndex& a_dit,
                 const Interval&  a_variables);
};

#include <CD_NamespaceFooter.H>

#endif
