/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBReflux.H
  @brief  Declaration of a class which can reflux over the coarse-fine interface. 
  @author Robert Marskar
*/

#ifndef CD_EBReflux_H
#define CD_EBReflux_H

// Std includes
#include <map>

// Chombo includes
#include <Copier.H>
#include <EBLevelGrid.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which can do refluxing across a coarse-fine interface. 
  @details This is like Chombos EBReflux, but trimmed for being faster & cheaper. This class will only reflux on the
  coarse-fine interface. If the user wants to accumulate fluxes, he/she should do it outside of this routine and only call
  the reflux method at the very end (only relevant for subcycling).

  This reflux method operates slightly different from Chombo in that the coarse-grid fluxes are first constructed by 
  conservatively coarsening the fine-grid fluxes. Then we do as we always do, subtracting the originan coarse-grid fluxes
  and adding in the corresponding fine-grid fluxes. Standard scaling factor for scaling the coarse/fine flux contributions
  are present, as always. 
*/
class EBReflux
{
public:
  /*!
    @brief Disallowed constructor. 
  */
  EBReflux() noexcept;

  /*!
    @brief Disallowed copy constructor.
  */
  EBReflux(const EBReflux& a_other) = delete;

  /*!
    @brief Full constructor. Calls the define function. 
    @param[in] a_eblg         Grids
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgCoFi     Coarsened fine grids
    @param[in] a_refRat       Refinement ratio between grids. 
  */
  EBReflux(const EBLevelGrid& a_eblg,
           const EBLevelGrid& a_eblgFine,
           const EBLevelGrid& a_eblgCoFi,
           const int          a_refRat) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBReflux() noexcept;

  /*!
    @brief Define method. Puts object in usable state. 
    @param[in] a_eblg     Grids on this level
    @param[in] a_eblgFine Fine grids
    @param[in] a_eblgCoFi Coarsened fine grids
    @param[in] a_refRat    Refinement ratio between coarse and fine grids.
  */
  virtual void
  define(const EBLevelGrid& a_eblg,
         const EBLevelGrid& a_eblgFine,
         const EBLevelGrid& a_eblgCoFi,
         const int          a_refRat) noexcept;

  /*!
    @brief Reflux into the coarse data. 
    @details This subtracts the coarse-fluxes from the EB and adds in the sum offine-grid fluxes. Multiplication by area 
    fractions and grid resolutions are done in this routine. 
    @param[inout] a_data      Data
    @param[in]    a_flux      Fluxes
    @param[in]    a_fineFlux  Fine-grid flux. 
    @param[in]    a_variables Variables to reflux.
    @param[in]    a_scaleCoarFlux Scaling constant for fine-grid fluxes
    @param[in]    a_scaleFineFlux Scaling constant for fine-grid fluxes
  */
  virtual void
  reflux(LevelData<EBCellFAB>&       a_Lphi,
         const LevelData<EBFluxFAB>& a_flux,
         const LevelData<EBFluxFAB>& a_fineFlux,
         const Interval              a_variables,
         const Real                  a_scaleCoarFlux,
         const Real                  a_scaleFineFlux) const noexcept;

protected:
  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Refinement factor between fine and coarse level. 
  */
  int m_refRat;

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Grid on this level. 
  */
  EBLevelGrid m_eblg;

  /*!
    @brief Coarsened fine grids. 
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Copier for copying from the coarsened fine grids to the coarse grids
  */
  Copier m_copier;

  /*!
    @brief Stencils for conservatively coarsening fine-grid fluxes onto the coarse grid
  */
  LayoutData<std::map<std::pair<int, Side::LoHiSide>, BaseIFFAB<FaceStencil>>> m_fluxCoarseningStencils;

  /*!
    @brief Regular coarse-fine regions for coarse-grid patches
  */
  LayoutData<std::map<std::pair<int, Side::LoHiSide>, DenseIntVectSet>> m_regularCoarseFineRegions;

  /*!
    @brief Stencils for conservatively coarsening fine-grid fluxes onto the coarse grid
  */
  mutable LayoutData<std::map<std::pair<int, Side::LoHiSide>, FaceIterator>> m_fluxCoarseningRegions;

  /*!
    @brief Cut-cells in each patch that lie on the coarse-fine interface
    @note Defined over the coarse grids
  */
  mutable LayoutData<std::map<std::pair<int, Side::LoHiSide>, VoFIterator>> m_irregularCoarseFineRegions;

  /*!
    @brief Define coarse-fine interface regions so that we can compute things later. 
  */
  virtual void
  defineRegionsCF() noexcept;

  /*!
    @brief Define the stencils required for coarsening the fine-grid fluxes onto the coarse grid
  */
  virtual void
  defineStencils() noexcept;

  /*!
    @brief Define buffers
  */
  virtual void
  defineBuffers() noexcept;

  /*!
    @brief Coarsen the input fine fluxes onto the coarse fluxes. 
    @details This replaces the coarse fluxes along the coarse-fine interface by the conservative average of the fine fluxes.
    @param[inout] a_coarFluxes Coarse-grid fluxes
    @param[in]    a_fineFluxes Fine-grid fluxes
    @param[in]    a_coarVar    Coarse-flux variable
    @param[in]    a_fineVar    Fine-flux variable
  */
  virtual void
  coarsenFluxesCF(LevelData<EBFluxFAB>&       a_coarFluxes,
                  const LevelData<EBFluxFAB>& a_fineFluxes,
                  const int                   a_coarVar,
                  const int                   a_fineVar) const noexcept;

  /*!
    @brief Reflux data into a_Lphi
    @details This subtracts the old fluxes (scaled by a_scaleOldFlux) from a_Lphi and adds in the new fluxes (scaled by a_scaleFineFlux).
    @param[inout] a_Lphi
    @param[in]    a_oldFluxes     Previous fluxes from which a_Lphi was computed
    @param[in]    a_newFluxes     New fluxes to be put into a_Lphi
    @param[in]    a_phiVar        Variable in a_Lphi to reflux
    @param[in]    a_oldFluxVar    Variable in a_oldFluxes to subtract
    @param[in]    a_newFluxVar    Variable in a_newFluxes to add
    @param[in]    a_scaleCoarFlux Scaling constant for fine-grid fluxes
    @param[in]    a_scaleFineFlux Scaling constant for fine-grid fluxes
  */
  virtual void
  refluxIntoCoarse(LevelData<EBCellFAB>&       a_Lphi,
                   const LevelData<EBFluxFAB>& a_oldFluxes,
                   const LevelData<EBFluxFAB>& a_newFluxes,
                   const int                   a_phiVar,
                   const int                   a_oldFluxVar,
                   const int                   a_newFluxVar,
                   const Real                  a_scaleCoarFlux,
                   const Real                  a_scaleFineFlux) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
