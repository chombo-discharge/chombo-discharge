/* chombo-discharge
 * Copyright © 2025 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBNonConservativeDivergence.H
  @brief  Declaration of a class for computing non-conservative divergences.
  @author Robert Marskar
*/

#ifndef CD_EBNonConservativeDivergence_H
#define CD_EBNonConservativeDivergence_H

// Chombo includes
#include <EBLevelGrid.H>
#include <LevelData.H>
#include <LayoutData.H>
#include <VoFIterator.H>
#include <Stencils.H>
#include <BaseIVFAB.H>
#include <EBCellFAB.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for computing a non-conservative by smooshing a conservative divergence.
  @details This class will compute a stencil s(phi) = sum(kappa * phi)/sum(kappa) where
  the stencil points consist of cells that can be reached through a monotone path.

  This is a per-level class without AMR functionality. On coarse-fine interfaces the user
  is responsible for filling the ghost cells and redistributing the mass correctly (see EBFluxRedistribution)
*/
class EBNonConservativeDivergence
{
public:
  /*!
    @brief Default constructor. Must call the define function afterwards
  */
  EBNonConservativeDivergence() noexcept;

  /*!
    @brief Copy constructor. Not allowed.
  */
  EBNonConservativeDivergence(const EBNonConservativeDivergence& a_other) = delete;

  /*!
    @brief Defining constructor. This calls the define function.
    @param[in] a_eblg Grids
    @param[in] a_radius Stencil radius
  */
  EBNonConservativeDivergence(const EBLevelGrid& a_eblg, const int a_radius) noexcept;

  /*!
    @brief Destructor. Does nothing.
  */
  virtual ~EBNonConservativeDivergence() noexcept;

  /*!
    @brief Define function. Puts object in usable state. 
    @param[in] a_eblg Grids
    @param[in] a_radius Stencil radius
  */
  virtual void
  define(const EBLevelGrid& a_eblg, const int a_radius) noexcept;

  /*!
    @brief Function for computing the non-conservative divergence.
    @param[out] a_nonConsDivF The non-conservative divergence.
    @param[in] a_kappaDivF The conservative divergence, but not weighted by kappa.
  */
  virtual void
  nonConservativeDivergence(LevelData<BaseIVFAB<Real>>& a_nonConsDivF,
                            const LevelData<EBCellFAB>& a_kappaDivF) const noexcept;

protected:
  /*!
    @brief Defined or not
  */
  bool m_isDefined;

  /*!
    @brief Grids
  */
  EBLevelGrid m_eblg;

  /*!
    @brief VoF iterator for iterating through the stencils
  */
  mutable LayoutData<VoFIterator> m_vofIterator;  

  /*!
    @brief Stencils.
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_stencils;
};

#include <CD_NamespaceFooter.H>

#endif
