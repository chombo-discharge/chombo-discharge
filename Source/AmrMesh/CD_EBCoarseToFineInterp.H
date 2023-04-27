/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBCoarseToFineInterp.H
  @brief  Declaration of an aggregated class for regrid operations.
  @author Robert Marskar
*/

#ifndef CD_EBCoarseToFineInterp_H
#define CD_EBCoarseToFineInterp_H

// Chombo includes
#include <RefCountedPtr.H>
#include <EBLevelGrid.H>

#include <EBPWLFineInterp.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for interpolating data to fine grids. Can use constant interpolation or include limiters.
*/
class EBCoarseToFineInterp
{
public:
  /*!
    @brief Type of interpolation methods supported.
    PWC = Piecewise constant, ignoring the embedded boundary. This is not conservative. 
    ConservativePWC = Piecewise constant, including the EB volume/area fractions boundary. 
    ConservativeMinMod = Min-mod slope for interpolation. 
    ConservativeMonotonizedCentral = Monotonized central slope for interpolation. 
  */
  enum Type
  {
    PWC,
    ConservativePWC,
    ConservativeMinMod,
    ConservativeMonotizedCentral,
  };

  /*!
    @brief Default constructor. Must subsequently call the define method. 
  */
  EBCoarseToFineInterp() noexcept;

  /*!
    @brief Constructor. Constructs a valid object which can be used for regridding.
    @param[in] a_eblgFine The fine grid. 
    @param[in] a_eblgCoFi Coarsened fine grid. 
    @param[in] a_eblgCoar The coarse grid. 
    @param[in] a_refRat   Refinement ratio between fine and coarse. 
    @param[in] a_nComp    Number of components. 
    @param[in] a_ebisPtr  EB index space. 
  */
  EBCoarseToFineInterp(const EBLevelGrid&        a_eblgFine,
                       const EBLevelGrid&        a_eblgCoFi,
                       const EBLevelGrid&        a_eblgCoar,
                       const int&                a_refRat,
                       const int&                a_nComp,
                       const EBIndexSpace* const a_ebisPtr) noexcept;

  /*!
    @brief Disallowed copy constructor
  */
  EBCoarseToFineInterp(const EBCoarseToFineInterp&) = delete;

  /*!
    @brief Disallowed move constructor
  */
  EBCoarseToFineInterp(const EBCoarseToFineInterp&&) = delete;

  /*!
    @brief Destructor
  */
  virtual ~EBCoarseToFineInterp() noexcept;

  /*!
    @brief Disallowed copy operator
  */
  EBCoarseToFineInterp&
  operator=(const EBCoarseToFineInterp&) = delete;

  /*!
    @brief Disallowed move operator
  */
  EBCoarseToFineInterp&
  operator=(const EBCoarseToFineInterp&&) = delete;

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
    @brief Interpolate to the fine grid using the specified method. 
    @param[out] a_fineData  Fine-grid data.
    @param[in]  a_coarData  Coarse data
    @param[in]  a_variables Variables to interpolate to the fine grid.
    @param[in]  a_interpType Interpolation type
  */
  virtual void
  interpolate(LevelData<EBCellFAB>&             a_fineData,
              const LevelData<EBCellFAB>&       a_coarData,
              const Interval&                   a_variables,
              const EBCoarseToFineInterp::Type& a_interpType) const noexcept;

  /*!
    @brief Interpolate to the fine grid using the specified method. 
    @param[out] a_fineData  Fine-grid data
    @param[in]  a_coarData  Coarse data
    @param[in]  a_variables Variables to interpolate to the fine grid.
    @param[in]  a_interpType Interpolation type
  */
  virtual void
  interpolate(LevelData<BaseIVFAB<Real>>&       a_fineData,
              const LevelData<BaseIVFAB<Real>>& a_coarData,
              const Interval&                   a_variables,
              const EBCoarseToFineInterp::Type& a_interpType) const noexcept;

protected:
  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Refinement factor
  */
  int m_refRat;

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
    @brief Irregular grid region.
  */
  LayoutData<IntVectSet> m_irregRegions;

  /*!
    @brief Weights in cut-cells when regridding conservatively
    @details Defined over the fine layout. 
  */
  LevelData<BaseIVFAB<Real>> m_conservativeWeights;

  /*!
    @brief Chombo object. Will be taken out of the design.
  */
  mutable RefCountedPtr<EBPWLFineInterp> m_fineInterp;

  /*!
    @brief Data holder for holding data on the coarsened fine grids. 
    @note I'm letting this live as a permanent buffer, simply because its memory footprint is so small. 
  */
  mutable LevelData<BaseIVFAB<Real>> m_irregCoFi;

  /*!
    @brief Iterator for valid cut-cells. 
    @details This iterator runs over the valid fine cells.
  */
  mutable LayoutData<VoFIterator> m_fineVoFs;

  /*!
    @brief Define weights for cut-cell conservative regridding
  */
  virtual void
  defineWeights() noexcept;

  /*!
    @brief Piecewise constant interpolation of an EBCellFAB. This ignores the volume fraction. This sets
    the fine-grid data to the coarse-grid data, so it is not conservative. 
    @param[out] a_fineData Fine-grid data. 
    @param[in]  a_coarData Coarse data
    @param[in]  a_dit Grid index. 
    @param[in]  a_fineVar Fine-grid variable
    @param[in]  a_coarVar Coarse-grid variable
  */
  virtual void
  interpolatePWC(EBCellFAB&       a_fineData,
                 const EBCellFAB& a_coarData,
                 const DataIndex& a_dit,
                 const int&       a_fineVar,
                 const int&       a_coarVar) const noexcept;

  /*!
    @brief Piecewise constant interpolation of an EBCellFAB. Does not use slopes
    @param[out] a_fineData Piecewise constant interpolation of the coarse data.
    @param[in]  a_coarData Coarse data
    @param[in]  a_dit Grid index. 
    @param[in]  a_fineVar Fine-grid variable
    @param[in]  a_coarVar Coarse-grid variable
  */
  virtual void
  interpolateConservativePWC(EBCellFAB&       a_fineData,
                             const EBCellFAB& a_coarData,
                             const DataIndex& a_dit,
                             const int&       a_fineVar,
                             const int&       a_coarVar) const noexcept;

  /*!
    @brief Piecewise constant interpolation of an EB-data. This ignores the area fraction. This sets
    the fine-grid data to the coarse-grid data. 
    @param[out] a_fineData Fine-grid data. 
    @param[in]  a_coarData Coarse data
    @param[in]  a_dit Grid index. 
    @param[in]  a_fineVar Fine-grid variable
    @param[in]  a_coarVar Coarse-grid variable
  */
  virtual void
  interpolatePWC(BaseIVFAB<Real>&       a_fineData,
                 const BaseIVFAB<Real>& a_coarData,
                 const DataIndex&       a_dit,
                 const int&             a_fineVar,
                 const int&             a_coarVar) const noexcept;

  /*!
    @brief Piecewise constant interpolation of an EB-data. This includes the area fraction and is therefore
    conservative. 
    @param[out] a_fineData Fine-grid data. 
    @param[in]  a_coarData Coarse data
    @param[in]  a_dit Grid index. 
    @param[in]  a_fineVar Fine-grid variable
    @param[in]  a_coarVar Coarse-grid variable
  */
  virtual void
  interpolateConservativePWC(BaseIVFAB<Real>&       a_fineData,
                             const BaseIVFAB<Real>& a_coarData,
                             const DataIndex&       a_dit,
                             const int&             a_fineVar,
                             const int&             a_coarVar) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
