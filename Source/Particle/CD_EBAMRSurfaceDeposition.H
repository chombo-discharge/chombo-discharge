/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRSurfaceDeposition.H
  @brief  Declaration of a class for handling surface deposition of particles with EB and AMR
  @author Robert Marskar
*/

#ifndef CD_EBAMRSurfaceDeposition_H
#define CD_EBAMRSurfaceDeposition_H

// Chombo includes
#include <RefCountedPtr.H>
#include <Copier.H>
#include <EBLevelGrid.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief class for handling surface deposition of particles with EB and AMR.
*/
class EBAMRSurfaceDeposition
{
public:
  /*!
    @brief Default constructor. Leaves the object in undefined state and requires the user to call the define function
  */
  EBAMRSurfaceDeposition() noexcept;

  /*!
    @brief Full constructor. Calls the define function.
    @param[in] a_ebGrids              Grids
    @param[in] a_ebGridsCoarsenedFine Coarsened fine grids
    @param[in] a_ebGridsRefinedCoar   Refined coarse grids
    @param[in] a_refRat               Refinement factory between levels
    @param[in] a_dx                   Grid resolutions.
    @πaram[in] a_probLo               Lower-left corner of physical domain.
    @param[in] a_finestLevel          Finest grid level.
    @param[in] a_radius               Deposition radius
  */
  EBAMRSurfaceDeposition(const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGrids,
                         const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsCoarsenedFine,
                         const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsRefinedCoar,
                         const Vector<int>&                        a_refRat,
                         const Vector<Real>&                       a_dx,
                         const RealVect&                           a_probLo,
                         const int                                 a_finestLevel,
                         const int                                 a_radius) noexcept;

  /*!
    @brief Destructor
  */
  virtual ~EBAMRSurfaceDeposition() noexcept;

  /*!
    @brief Define function. Puts objects in usable state. 
    @param[in] a_ebGrids              Grids
    @param[in] a_ebGridsCoarsenedFine Coarsened fine grids
    @param[in] a_ebGridsRefinedCoar   Refined coarse grids
    @param[in] a_refRat               Refinement factory between levels
    @param[in] a_dx                   Grid resolutions.
    @πaram[in] a_probLo               Lower-left corner of physical domain.
    @param[in] a_finestLevel          Finest grid level.
    @param[in] a_radius               Deposition radius
  */
  virtual void
  define(const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGrids,
         const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsCoarsenedFine,
         const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsRefinedCoar,
         const Vector<int>&                        a_refRat,
         const Vector<Real>&                       a_dx,
         const RealVect&                           a_probLo,
         const int                                 a_finestLevel,
         const int                                 a_radius) noexcept;

  /*!
    @brief Deposit function. Deposits particle on surface.
    @details P is the particle type
    @note This will throw an error if a particle is located on a grid cell that is NOT an irregular cell
  */
  template <class P, const Real& (P::*particleScalarField)() const>
  void
  deposit(EBAMRIVData& a_meshData, const ParticleContainer<P>& a_particles) const noexcept;

  /*!
    @brief Deposit function. Deposits particle on surface.
    @details P is the particle type
    @note This will throw an error if a particle is located on a grid cell that is NOT an irregular cell
  */
  template <class P, Real (P::*particleScalarField)()>
  void
  deposit(EBAMRIVData& a_meshData, const ParticleContainer<P>& a_particles) const noexcept;

protected:
  /*!
    @brief Debug or not
  */
  bool m_debug;

  /*!
    @brief Is defined or not.
  */
  bool m_isDefined;

  /*!
    @brief Turn on/off chattiness
  */
  bool m_verbose;

  /*!
    @brief Lower-left corner of computational domain
  */
  RealVect m_probLo;

  /*!
    @brief Finest grid level
  */
  int m_finestLevel;

  /*!
    @brief Deposition radius
  */
  int m_radius;

  /*!
    @brief Grids on each level
  */
  Vector<RefCountedPtr<EBLevelGrid>> m_ebGrids;

  /*!
    @brief Coarsened grids. 
    @note m_ebGridsCoarsenedFine[i] contains the coarsening of m_ebGrids[i] and m_ebGridsCoarsenedFine[0] contains the nullpointer
  */
  Vector<RefCountedPtr<EBLevelGrid>> m_ebGridsCoarsenedFine;

  /*!
    @brief Refined grids. 
    @note m_ebGridsCoarsenedFine[i] contains the refinement of m_ebGrids[i] and m_ebGridsCoarsenedFine[m_finestLevel-1] contains the nullpointer
  */
  Vector<RefCountedPtr<EBLevelGrid>> m_ebGridsRefinedCoar;

  /*!
    @brief Refinement ratios between levels
  */
  Vector<int> m_refRat;

  /*!
    @brief Grid resolutions.
  */
  Vector<Real> m_dx;

  /*!
    @brief Deposition weights for each irregular cell. This contains 1.0/sum(alpha * dx^(D-1)) and the cells 
  */
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil>>>> m_depositionStencils;

  /*!
    @brief Stencils for interpolating coarse-grid data to the fine grid. 
    @details This is defined on the coarse grid but reaches into the refined coarse grid. 
  */
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil>>>> m_coarseToFineStencils;

  /*!
    @brief Stencils for conservatively coarsening fine-grid data (in ghost cells) to the coarse grid.
    @details Defined on the fine grid but reaches into the refined coarse grid.
  */
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil>>>> m_fineToCoarseStencils;

  /*!
    @brief Data on each level. 
    @details Allocating this separately because we define Copiers that need a specific amount of ghost cells.
  */
  mutable EBAMRIVData m_data;

  /*!
    @brief Data defined on a refinement of the coarse grid. 
  */
  mutable EBAMRIVData m_refinedCoarData;

  /*!
    @brief Coarsening of the fine-grid data
  */
  mutable EBAMRIVData m_coarsenedFineData;

  /*!
    @brief Valid to valid copiers
  */
  Vector<Copier> m_validToValidCopiers;

  /*!
    @brief Level copier for adding data from valid+ghost -> valid. 
    @details Needed when adding ghost-deposited data back into valid cells
  */
  Vector<Copier> m_copierLevel;

  /*!
    @brief Copier for copying data from the refined coarse level to the fine level.
  */
  Vector<Copier> m_copierRefinedCoarToFineNoGhosts;

  /*!
    @brief Copier for copying data from the coarsened fine level to the coarse level.
  */
  Vector<Copier> m_copierCoarsenedFineToCoar;

  /*!
    @brief Define buffer data.
  */
  virtual void
  defineBuffers() noexcept;

  /*!
    @brief Define copiers
  */
  virtual void
  defineDataMotion() noexcept;

  /*!
   @brief Define deposition weights
  */
  virtual void
  defineDepositionStencils() noexcept;

  /*!
    @brief Define stencils for interpolating coarse-grid data to the fine grid
  */
  virtual void
  defineCoarseToFineStencils() noexcept;

  /*!
    @brief Define stencils for coarsening data in fine-grid ghost cells to the coarse grid.
  */
  virtual void
  defineFineToCoarseStencils() noexcept;

  /*!
    @brief Add the invalid coarse data (i.e., underneath the fine grid) to the fine level. 
    @details This does conservative interpolation to the fine grid. This all takes place on our buffer storage.
    @param[inout] a_meshData EB surface data
  */
  virtual void
  addInvalidCoarseDataToFineData() const noexcept;

  /*!
    @brief Add the ghosted fine-level data to the coarse data. 
    @details This does conservative coarsening of the fine-grid ghosted data to the coarse grid. This all takes place on our buffer storage.
  */
  virtual void
  addFineGhostDataToValidCoarData() const noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_EBAMRSurfaceDepositionImplem.H>

#endif
