/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_AmrMesh.H
  @brief  Declaration of core class for handling AMR-related operations (with embedded boundaries)
  @author Robert Marskar
  @todo   a_lmax is a dead parameter in buildGrids. Fix this and the dependencies. 
  @todo   a_regsize is a dead parameter which should be removed. 
*/

#ifndef CD_AmrMesh_H
#define CD_AmrMesh_H

// Chombo includes
#include <DisjointBoxLayout.H>
#include <EBFluxRegister.H>
#include <EBLevelRedist.H>
#include <EBCoarToFineRedist.H>
#include <EBMGInterp.H>
#include <EBCoarToCoarRedist.H>
#include <EBFineToCoarRedist.H>
#include <EBPWLFineInterp.H>
#include <ProblemDomain.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_EbCoarAve.H>
#include <CD_ComputationalGeometry.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_IrregAmrStencil.H>
#include <CD_LoadBalancing.H>
#include <CD_MFLevelGrid.H>
#include <CD_DomainFluxIFFAB.H>
#include <CD_DomainFluxIFFABFactory.H>
#include <CD_ParticleContainer.H>
#include <CD_NonConservativeDivergenceStencil.H>
#include <CD_Realm.H>
#include <CD_LoadBalancing.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for handling spatial operations. 
  @details This class handles most (if not all) spatial operations. Examples are grid generation, coarsening of data, interpolation of data, memory allocations etc. 
*/
class AmrMesh
{
public:
  /*!
    @brief Default constructor. 
  */
  AmrMesh();

  /*!
    @brief Disallowed copy constructor
    @param[in] a_other Other AmrMesh
  */
  AmrMesh(const AmrMesh& a_other) = delete;

  /*!
    @brief Disallowed move constructor
    @param[in] a_other Other AmrMesh
  */
  AmrMesh(const AmrMesh&& a_other) = delete;

  /*!
    @brief Destructor
  */
  ~AmrMesh();

  /*!
    @brief Disallowed copy assignment
    @param[in] a_other Other AmrMesh
  */
  AmrMesh&
  operator=(const AmrMesh& a_other) = delete;

  /*!
    @brief Disallowed move assignment
    @param[in] a_other Other AmrMesh
  */
  AmrMesh&
  operator=(const AmrMesh&& a_other) = delete;

  /*!
    @brief Deallocate data. 
    @param[inout] a_data Data to be deleted. 
  */
  template <typename T>
  void
  deallocate(Vector<T*>& a_data) const;

  /*!
    @brief Deallocate smart pointer data.
    @param[inout] a_data Data to be deallocated. 
  */
  template <typename T>
  void
  deallocate(Vector<RefCountedPtr<T>>& a_data) const;

  /*!
    @brief Deallocate data.
    @param[inout] a_data Data to be deallocated.
  */
  template <typename T>
  void
  deallocate(EBAMRData<T>& a_data) const;

  /*!
    @brief Turn smart-pointer data structure into regular-pointer data structure
    @param[out] a_alias Raw pointer aliased 
    @param[in]  a_data  Smart-pointer data
    @details This fetches the bare pointer from a_data and sets it in a_alias. 
  */
  template <typename T>
  void
  alias(Vector<T*>& a_alias, const Vector<RefCountedPtr<T>>& a_data) const;

  /*!
    @brief Turn smart-pointer data structure into regular-pointer data structure
    @param[out] a_alias Raw pointer aliased 
    @param[in]  a_data  Smart-pointer data
    @details This fetches the bare pointer from a_data and sets it in a_alias. 
  */
  template <typename T, typename S>
  void
  alias(Vector<T*>& a_alias, const EBAMRData<S>& a_data) const;

  /*!
    @brief Template class for generic allocation of particle data. 
    @param[out] a_particles Particle data holder. 
    @param[out] a_realm     Realm on which the particles will be allocated. 
  */
  template <typename T>
  void
  allocate(Vector<RefCountedPtr<ParticleData<T>>>& a_particles, const std::string a_realm) const;

  /*!
    @brief Template class for generic allocation of particle containers. 
    @param[out] a_container Particle container to be allocated
    @param[out] a_realm     Realm on which the particles will be allocated. 
    @note Particles that live within a_pvrBuffer on the fine side of the refinement boundary will be mapped to the coarse level. 
  */
  template <typename T>
  void
  allocate(ParticleContainer<T>& a_container, const std::string a_realm) const;

  /*!
    @brief Allocate pointer but not any memory blocks. 
    @param[out] a_data Data structure for which pointers will be constructed. 
    @details This sets all the elements in the vector to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void
  allocatePointer(Vector<RefCountedPtr<T>>& a_data) const;

  /*!
    @brief Allocate pointer but not any memory blocks. 
    @param[out] a_data       Data structure for which pointers will be constructed. 
    @param[in] a_finestLevel Finest grid level on which a_data will be defined. 
    @details This sets all the elements in the vector to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void
  allocatePointer(Vector<RefCountedPtr<T>>& a_data, const int a_finestLevel) const;

  /*!
    @brief Allocate pointer but not any memory blocks. 
    @param[out] a_data Data structure for which pointers will be constructed. 
    @details This sets all the elements in a_data to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void
  allocatePointer(EBAMRData<T>& a_data) const;

  /*!
    @brief Allocate pointer but not any memory.
    @param[out] a_data       Data structure for which pointers will be constructed. 
    @param[in] a_finestLevel Finest grid level on which a_data will be defined. 
    @details This sets all the elements in the vector to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void
  allocatePointer(EBAMRData<T>& a_data, const int a_finestLevel) const;

  /*!
    @brief Parse options. Called during the constructor. 
  */
  void
  parseOptions();

  /*!
    @brief Parse runtime options
    @note Called by Driver as simulations progress. 
  */
  void
  parseRuntimeOptions();

  /*!
    @brief Register an operator over a realm and a phase. 
    @param[in] a_operator Operator to be registered. 
    @param[in] a_realm    Realm name
    @param[in] a_phase    Phase (gas or solid). 
    @note Run-time errors will occur if the operator or realm do not exist. 
  */
  void
  registerOperator(const std::string a_operator, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Register a boolean mask over a realm. 
    @param[in] a_mask     Mask name to be registered. 
    @param[in] a_phase    Phase (gas or solid). 
    @param[in] a_buffer   Relevant buffer for the mask.
    @note The meaning of a_buffer differs for different masks.
    @note Run-time errors will occur if the mask or realm do not exist.
  */
  void
  registerMask(const std::string a_mask, const int a_buffer, const std::string a_realm);

  /*!
    @brief Build domains
  */
  void
  buildDomains();

  /*!
    @brief Set the finest level
    @param[in] a_finestLevel
    @note This is dangerous routine, but it's needed because Driver reads checkpoint files and needs to set the finest level from those. 
  */
  void
  setFinestLevel(const int a_finestLevel);

  /*!
    @brief Set the coarsest grid cells
  */
  void
  setCoarsestGrid(const IntVect& a_nCells);

  /*!
    @brief Query if a realm exists
    @param[in] a_realm Name of the realm. 
    @return True if the realm exists and false otherwise. 
  */
  bool
  queryRealm(const std::string a_realm) const;

  /*!
    @brief Check if mesh has an EBCF
    @note Always returns true, non-EBCF cases are a thing of the past. 
  */
  bool
  getEbCf() const;

  /*!
    @brief Register a new realm
    @param[in] a_realm Realm name
    @note It is safe to re-register a realm. 
  */
  void
  registerRealm(const std::string a_realm);

  /*!
    @brief Alias multifluid data up to a specified level. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @param[in]  a_finestLevel Maximum grid level for aliasing operation. 
  */
  void
  alias(EBAMRCellData&           a_data,
        const phase::which_phase a_phase,
        const MFAMRCellData&     a_mfdata,
        const int                a_finestLevel) const;

  /*!
    @brief Alias multifluid data up to a specified level. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @param[in]  a_finestLevel Maximum grid level for aliasing operation. 
  */
  void
  alias(EBAMRFluxData&           a_data,
        const phase::which_phase a_phase,
        const MFAMRFluxData&     a_mfdata,
        const int                a_finestLevel) const;

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
  */
  void
  alias(EBAMRCellData& a_data, const phase::which_phase a_phase, const MFAMRCellData& a_mfdata) const;

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
  */
  void
  alias(EBAMRFluxData& a_data, const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata) const;

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
  */
  void
  alias(EBAMRIVData& a_data, const phase::which_phase a_phase, const MFAMRIVData& a_mfdata) const;

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @return A single-phase alias of the multifluid data, having the same size as the input data. 
  */
  EBAMRCellData
  alias(const phase::which_phase a_phase, const MFAMRCellData& a_mfdata) const;

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @return A single-phase alias of the multifluid data, having the same size as the input data. 
  */
  EBAMRFluxData
  alias(const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata) const;

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @return A single-phase alias of the multifluid data, having the same size as the input data. 
  */
  EBAMRIVData
  alias(const phase::which_phase a_phase, const MFAMRIVData& a_mfdata) const;

  /*!
    @brief Compute cell-centered gradient over an AMR hierarchy. 
    @param[out] a_gradient Cell centered gradient. 
    @param[in]  a_phi      The scalar for which the gradient is computed. 
    @param[in]  a_realm    Name of the realm where the data lives. 
    @param[in]  a_phase    Phase on which the data lives. 
    @note This routine will reach into ghost cells and across refinement boundaries. The user must make sure that 
    ghost cells are updated before using this routine. 
  */
  void
  computeGradient(EBAMRCellData&           a_gradient,
                  const EBAMRCellData&     a_phi,
                  const std::string        a_realm,
                  const phase::which_phase a_phase) const;

  /*!
    @brief Compute face-centered gradient over an AMR hierarchy. 
    @param[out] a_gradient Face centered gradient. 
    @param[in]  a_phi      The scalar for which the gradient is computed. 
    @param[in]  a_realm    Name of the realm where the data lives. 
    @param[in]  a_phase    Phase on which the data lives. 
    @note This routine will reach into ghost cells and across refinement boundaries. The user must make sure that 
    ghost cells are updated before using this routine. 
  */
  void
  computeGradient(EBAMRFluxData&           a_gradient,
                  const EBAMRCellData&     a_phi,
                  const std::string        a_realm,
                  const phase::which_phase a_phase) const;

  /*!
    @brief Compute cell-centered gradient over an AMR hierarchy. 
    @param[out] a_gradient Cell centered gradient. 
    @param[in]  a_phi      The scalar for which the gradient is computed. 
    @param[in]  a_realm    Name of the Realm where the data lives. 
    @note This routine will reach into ghost cells and across refinement boundaries. The user must make sure that 
    ghost cells are updated before using this routine. 
  */
  void
  computeGradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi, const std::string a_realm) const;

  /*!
    @brief Compute face-centered gradient over an AMR hierarchy. 
    @param[out] a_gradient Face centered gradient. 
    @param[in]  a_phi      The scalar for which the gradient is computed. 
    @param[in]  a_realm    Name of the Realm where the data lives. 
    @note This routine will reach into ghost cells and across refinement boundaries. The user must make sure that 
    ghost cells are updated before using this routine. 
  */
  void
  computeGradient(MFAMRFluxData& a_gradient, const MFAMRCellData& a_phi, const std::string a_realm) const;

  /*!
    @brief Allocate a data holder over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_phase  Phase (gas or solid)
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void
  allocate(EBAMRCellData&           a_data,
           const std::string        a_realm,
           const phase::which_phase a_phase,
           const int                a_nComp,
           const int                a_nGhost = -1) const;

  /*!
    @brief Allocate a data holder over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_phase  Phase (gas or solid)
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void
  allocate(EBAMRFluxData&           a_data,
           const std::string        a_realm,
           const phase::which_phase a_phase,
           const int                a_nComp,
           const int                a_nGhost = -1) const;

  /*!
    @brief Allocate a data holder over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_phase  Phase (gas or solid)
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void
  allocate(EBAMRIVData&             a_data,
           const std::string        a_realm,
           const phase::which_phase a_phase,
           const int                a_nComp,
           const int                a_nGhost = 0) const;

  /*!
    @brief Allocate a data holder over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void
  allocate(EBAMRIFData&             a_data,
           const std::string        a_realm,
           const phase::which_phase a_phase,
           const int                a_nComp,
           const int                a_nGhost = 0) const;

  /*!
    @brief Allocate Boolean data over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void
  allocate(EBAMRBool& a_data, const std::string a_realm, const int a_nComp, const int a_nGhost = -1) const;

  /*!
    @brief Allocate a data holder over a specific realm.
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void
  allocate(MFAMRCellData& a_data, const std::string a_realm, const int a_nComp, const int a_nGhost = -1) const;

  /*!
    @brief Allocate a data holder over a specific realm.
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void
  allocate(MFAMRFluxData& a_data, const std::string a_realm, const int a_nComp, const int a_ghost = -1) const;

  /*!
    @brief Allocate a data holder over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void
  allocate(MFAMRIVData& a_data, const std::string a_realm, const int a_nComp, const int a_ghost = 0) const;

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_lmin) const;

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_lmin) const;

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_lmin) const;

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_lmin) const;

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(EBAMRBool& a_data, const int a_lmin) const;

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(MFAMRCellData& a_data, const int a_lmin) const;

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(MFAMRFluxData& a_data, const int a_lmin) const;

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(MFAMRIVData& a_data, const int a_lmin) const;
  /*!
    @brief Average down multifluid data over a realm. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
  */
  void
  averageDown(MFAMRCellData& a_data, const std::string a_realm) const;

  /*!
    @brief Average down multifluid data over a realm. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
  */
  void
  averageDown(MFAMRFluxData& a_data, const std::string a_realm) const;

  /*!
    @brief Average multifluid data over a realm. 
    @details This computes an arithmetic average of the face data. Data on coarse faces is replaced by the arithmetic average of the fine face data. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
  */
  void
  averageFaces(MFAMRFluxData& a_data, const std::string a_realm) const;

  /*!
    @brief Average down on specific realm and phase. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  averageDown(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Average down on specific realm and phase. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
    @param[in]    a_level The coarse level.
    @details This coarsens data from level a_lvl+1 to level a_lvl
  */
  void
  averageDown(EBAMRCellData&           a_data,
              const std::string        a_realm,
              const phase::which_phase a_phase,
              const int                a_lvl) const;

  /*!
    @brief Coarsen cells over a realm. This does a specific level. 
    @param[inout] a_coarData Coarse grid data
    @param[inout] a_fineData Fine grid data
    @param[in]    a_level    The grid level corresponding to a_fineData
    @param[in]    a_realm    Realm name
    @param[in]    a_phase    Phase (gas or solid)
  */
  void
  averageDown(LevelData<EBCellFAB>&       a_coarData,
              const LevelData<EBCellFAB>& a_fineData,
              const int                   a_level,
              const std::string           a_realm,
              const phase::which_phase    a_phase) const;

  /*!
    @brief Average down on specific realm and phase. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  averageDown(EBAMRFluxData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Average down on specific realm and phase. 
    @details This computes an arithmetic average of the face data. Data on coarse faces is replaced by the arithmetic average of the fine face data. 
    @param[inout] a_data  Data to be average.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  averageFaces(EBAMRFluxData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Average down on specific realm and phase. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  averageDown(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Do conservative coarsening for data defined on cut-cells. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  conservativeAverage(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Deposit scalar particle quantities on the mesh. 
    @details This routine will switch between various algorithms depending on the user input. If the user chooses to deposit with "halos" near the refinement boundaries
    then this function requires the "halo" particles to be filled before entering this routine. Thus, if depositing with halos (which requires kernels with 2x or 4x 
    the particle width), the user must first call ParticleContainer::copyMaskParticles(AMRMask). If the user does not do this, he will lose mass near the refinement
    boundaries. Other than that, this function is just like regular particle-mesh deposition. The mesh field should have exactly one component and the quantity to be
    deposited is fetched from the specified template parameter 'a_particleScalarField'. This should be a pointer to a member function with signature 
    const Real& P::function() const. E.g. if we are depositing mass through a function const Real& P::mass() const we will specifiy
    
    interpolate<P, &P::mass> (...)
	
    To deposit a different quantity, simply change the function pointer in the input parameter. 

    @param[out] a_meshData             Mesh data. Must have exactly one compnent. 
    @param[in]  a_realm                Realm where data is registered.
    @param[in]  a_phase                Phase where data is registered.
    @param[in]  a_particles            Particle container. Must be in "usable state" for deposition.
    @param[in]  a_depositionType       Specification of deposition kernel (e.g., CIC)
    @param[in]  a_coarseFineDeposition Specification of handling of coarse-fine boundaries.
    @param[in]  a_forceIrregNGP        Force NGP deposition in irregular cells or not. 
  */
  template <class P, const Real& (P::*particleScalarField)() const>
  void
  depositParticles(EBAMRCellData&              a_meshData,
                   const std::string&          a_realm,
                   const phase::which_phase&   a_phase,
                   const ParticleContainer<P>& a_particles,
                   const DepositionType        a_depositionType,
                   const CoarseFineDeposition  a_coarseFineDeposition,
                   const bool                  a_forceIrregNGP = false);

  /*!
    @brief Deposit scalar particle quantities on the mesh. 
    @details Precisely like the above function, except that the function signature is Real P::particleScalarField() const (C++ is not too smart about these things)
    @param[out] a_meshData             Mesh data. Must have exactly one compnent. 
    @param[in]  a_realm                Realm where data is registered.
    @param[in]  a_phase                Phase where data is registered.
    @param[in]  a_particles            Particle container. Must be in "usable state" for deposition.
    @param[in]  a_depositionType       Specification of deposition kernel (e.g., CIC)
    @param[in]  a_coarseFineDeposition Specification of handling of coarse-fine boundaries.
    @param[in]  a_forceIrregNGP        Force NGP deposition in irregular cells or not. 
  */
  template <class P, Real (P::*particleScalarField)() const>
  void
  depositParticles(EBAMRCellData&              a_meshData,
                   const std::string&          a_realm,
                   const phase::which_phase&   a_phase,
                   const ParticleContainer<P>& a_particles,
                   const DepositionType        a_depositionType,
                   const CoarseFineDeposition  a_coarseFineDeposition,
                   const bool                  a_forceIrregNGP = false);

  /*!
    @brief Deposit vector particle quantities on the mesh. 
    @details This routine will switch between various algorithms depending on the user input. If the user chooses to deposit with "halos" near the refinement boundaries
    then this function requires the "halo" particles to be filled before entering this routine. Thus, if depositing with halos (which requires kernels with 2x or 4x 
    the particle width), the user must first call ParticleContainer::copyMaskParticles(AMRMask). If the user does not do this, he will lose mass near the refinement
    boundaries. Other than that, this function is just like regular particle-mesh deposition. The mesh field should have exactly SpaceDim components and the quantity to be
    deposited is fetched from the specified template parameter 'a_particleVectorField'. This should be a pointer to a member function with signature 
    const RealVect& P::function() const. E.g. if we are depositing current through a function const RealVect& P::current() const we will specifiy
    
    interpolate<P, &P::current> (...)
	
    To deposit a different quantity, simply change the function pointer in the input parameter. 

    @param[out] a_meshData             Mesh data. Must have exactly one compnent. 
    @param[in]  a_realm                Realm where data is registered.
    @param[in]  a_phase                Phase where data is registered.
    @param[in]  a_particles            Particle container. Must be in "usable state" for deposition.
    @param[in]  a_depositionType       Specification of deposition kernel (e.g., CIC)
    @param[in]  a_coarseFineDeposition Specification of handling of coarse-fine boundaries.
    @param[in]  a_forceIrregNGP        Force NGP deposition in irregular cells or not. 
  */
  template <class P, const RealVect& (P::*particleVectorField)() const>
  void
  depositParticles(EBAMRCellData&              a_meshData,
                   const std::string&          a_realm,
                   const phase::which_phase&   a_phase,
                   const ParticleContainer<P>& a_particles,
                   const DepositionType        a_depositionType,
                   const CoarseFineDeposition  a_coarseFineDeposition,
                   const bool                  a_forceIrregNGP = false);

  /*!
    @brief Deposit vector particle quantities on the mesh. 
    @details Precisely like the above function, except that the function signature is RealVect P::particleVectorField() const (C++ is not too smart about these things)
    @param[out] a_meshData             Mesh data. Must have exactly one compnent. 
    @param[in]  a_realm                Realm where data is registered.
    @param[in]  a_phase                Phase where data is registered.
    @param[in]  a_particles            Particle container. Must be in "usable state" for deposition.
    @param[in]  a_depositionType       Specification of deposition kernel (e.g., CIC)
    @param[in]  a_coarseFineDeposition Specification of handling of coarse-fine boundaries.
    @param[in]  a_forceIrregNGP        Force NGP deposition in irregular cells or not. 
  */
  template <class P, RealVect (P::*particleVectorField)() const>
  void
  depositParticles(EBAMRCellData&              a_meshData,
                   const std::string&          a_realm,
                   const phase::which_phase&   a_phase,
                   const ParticleContainer<P>& a_particles,
                   const DepositionType        a_depositionType,
                   const CoarseFineDeposition  a_coarseFineDeposition,
                   const bool                  a_forceIrregNGP = false);

  /*!
    @brief Interpolate a scalar field onto the particle position. 
    @details This is just like regular particle-mesh interpolation. The input field should have exactly one component and the
    the field will be interpolated onto the template parameter's input field. The template parameter 'particleScalarField' should be a pointer
    to a member function which will set the particle field. E.g. the function must have a signature Real& P::particleScalarField(). A valid
    expression is e.g.
    
        interpolate<P, &P::mass> (...)
	
    To interpolate onto a different field, replace ::mass by the other scalar field.

    @param[inout] a_particles            Particles to be interpolated. 
    @param[in]    a_realm                Realm where data is registered.
    @param[in]    a_phase                Phase where data is registered.
    @param[in]    a_meshScalarField      Scalar field on the mesh 
    @param[in]    a_interpType           Interpolation type. 
    @param[in]    a_forceIrregNGP        Force NGP interpolation in cut-cells. 
  */
  template <class P, Real& (P::*particleScalarField)()>
  void
  interpolateParticles(ParticleContainer<P>&     a_particles,
                       const std::string&        a_realm,
                       const phase::which_phase& a_phase,
                       const EBAMRCellData&      a_meshScalarField,
                       const DepositionType      a_interpType,
                       const bool                a_forceIrregNGP = false) const;

  /*!
    @brief Interpolate a vector field onto the particle position. 
    @details This is just like regular particle-mesh interpolation. The input field should have exactly SpaceDim components and the
    the field will be interpolated onto the template parameter's input field. The template parameter 'particleVectorField' should be a pointer
    to a member function which will set the particle field. E.g. the function must have a signature RealVect& P::particleVectorField(). A valid
    expression is e.g.
    
        interpolate<P, &P::velocity> (...)
	
    To interpolate onto a different field, replace ::velocity by another vector field. 

    @param[inout] a_particles            Particles to be interpolated. 
    @param[in]    a_realm                Realm where data is registered.
    @param[in]    a_phase                Phase where data is registered.
    @param[in]    a_meshScalarField      Scalar field on the mesh 
    @param[in]    a_interpType           Interpolation type. 
    @param[in]    a_forceIrregNGP        Force NGP interpolation in cut-cells. 
  */
  template <class P, RealVect& (P::*particleVectorField)()>
  void
  interpolateParticles(ParticleContainer<P>&     a_particles,
                       const std::string&        a_realm,
                       const phase::which_phase& a_phase,
                       const EBAMRCellData&      a_meshVectorField,
                       const DepositionType      a_interpType,
                       const bool                a_forceIrregNGP = false) const;

  /*!
    @brief Function which removes particles from the domain if they fall inside the EB.
    @details The template parameter indicate the particle type -- it MUST have a const RealVect& P::position() const function. 
    @details This version uses the implicit function to evaluate whether or not the particles are inside the EB. The particles will be removed from
    the container if f(x) > a_tolerance*dx. 
    @param[inout] a_particles Particle data
    @param[in]    a_phase     Phase where the particles live. 
    @param[in]    a_tolerance Tolerance. 
    @note Because the implicit functions are global, this function will work even if the particles aren't mapped to the correct EBISBox.
  */
  template <class P>
  void
  removeCoveredParticlesIF(ParticleContainer<P>&     a_particles,
                           const phase::which_phase& a_phase,
                           const Real                a_tolerance = 0.0) const;

  /*!
    @brief Function which removes particles from the domain if they fall inside the EB.
    @details The template parameter indicate the particle type -- it MUST have a const RealVect& P::position() const function. 
    @details This version uses the discrete information to evaluate whether or not the particles are inside the EB. The particles will be removed from
    the container if they lie within a covered cell OR if they lie in an irregular cell but on the "wrong" side of the EB face (to within precision tolerance*dx)
    @param[inout] a_particles Particle data
    @param[in]    a_phase     Phase where the particles live. 
    @param[in]    a_tolerance Tolerance when to remove particles
    @note This function will only work if the particles are mapped to the correct EBISBox. 
  */
  template <class P>
  void
  removeCoveredParticlesDiscrete(ParticleContainer<P>&     a_particles,
                                 const phase::which_phase& a_phase,
                                 const Real                a_tolerance = 0.0) const;

  /*!
    @brief Function which removes particles from the domain if they fall inside the EB.
    @details The template parameter indicate the particle type -- it must have a const RealVect& P::position() function. 
    @details This version removes all particles that live in covered cells.
    @param[inout] a_particles Particle data
    @param[in]    a_phase     Phase where the particles live. 
    @note This function will only work if the particles are mapped to the correct EBISBox. 
  */
  template <class P>
  void
  removeCoveredParticlesVoxels(ParticleContainer<P>& a_particles, const phase::which_phase& a_phase) const;

  /*!
    @brief Function which transferse particles from one particle container to another if they fall inside the EB.
    @details The template parameter indicate the particle type -- it MUST have a const RealVect& P::position() const function. 
    @details This version uses the implicit function to evaluate whether or not the particles are inside the EB. The particles will be transferred from
    the container if f(x) > a_tolerance*dx.
    @param[inout] a_particlesFrom Container to transfer from
    @param[inout] a_particlesTo   Container to transfer to
    @param[in]    a_phase         Phase where the particles live. 
    @param[in]    a_tolerance     Tolerance (Note: relative to grid resolution)
    @note Because the implicit functions are global, this function will work even if the particles aren't mapped to the correct boxes (but remapping will be necessary)
  */
  template <class P>
  void
  transferCoveredParticlesIF(ParticleContainer<P>&     a_particlesFrom,
                             ParticleContainer<P>&     a_particlesTo,
                             const phase::which_phase& a_phase,
                             const Real                a_tolerance = 0.0) const;

  /*!
    @brief Function which transferse particles from one particle container to another if they fall inside the EB.
    @details The template parameter indicate the particle type -- it MUST have a const RealVect& P::position() const function. 
    @details This version uses the discrete information to evaluate whether or not the particles are inside the EB. The particles will be transferred from
    the container if they lie within a covered cell OR if they lie in an irregular cell but on the "wrong" side of the EB face (to within precision tolerance*dx)
    @param[inout] a_particlesFrom Container to transfer from
    @param[inout] a_particlesTo   Container to transfer to
    @param[in]    a_phase         Phase where the particles live. 
    @param[in]    a_tolerance     Tolerance. 
    @note This function will only work if the particles are mapped to the correct EBISBox. 
  */
  template <class P>
  void
  transferCoveredParticlesDiscrete(ParticleContainer<P>&     a_particlesFrom,
                                   ParticleContainer<P>&     a_particlesTo,
                                   const phase::which_phase& a_phase,
                                   const Real                a_tolerance = 0.0) const;

  /*!
    @brief Function which transferse particles from one particle container to another if they fall inside the EB.
    @details The template parameter indicate the particle type -- it MUST have a const RealVect& P::position() const function. 
    @details This version transfers all particles that live in covered cells.
    @param[inout] a_particlesFrom Container to transfer from
    @param[inout] a_particlesTo   Container to transfer to
    @param[in]    a_phase         Phase where the particles live. 
    @param[in]    a_tolerance     Tolerance. 
    @note This function will only work if the particles are mapped to the correct EBISBox. 
  */
  template <class P>
  void
  transferCoveredParticlesVoxels(ParticleContainer<P>&     a_particlesFrom,
                                 ParticleContainer<P>&     a_particlesTo,
                                 const phase::which_phase& a_phase) const;

  /*!
    @brief Particle intersection algorithm based on ray-casting.
    @details This routine will iterate through all the particles and check if they intersect the geometry. The template
    parameter indicates the particle type -- it MUST have const RealVec& position() const and  const RealVect& oldPosition() const 
    functions that determine the start and stop position of the particle trajectory. This routine uses a ray-casting method
    to check for intersections with the EB (the domain side is much easier). If the particles are closer to the EB than
    a_tolerance, they are absorbed and placed on the EB. Their position are updated and they are placed in the 
    a_ebParticles argument. This routine uses a ray-casting method where it computes the distance from the EB
    (assuming that the implicit function is a signed distance function). Particles are then moved that
    distance along their trajectory and we then update the new distance to the EB. This is done recursively until the particles
    have either moved the entire length or been absorbed by the EB or domain side. 
    @param[inout] a_activeParticles Particles to be intersected with geometry
    @param[out]   a_ebParticles     Particles that intersected with the EB
    @param[out]   a_domainParticles Particles that intersected with the domain faces
    @param[in]    a_phase           Phase where the input particles live
    @param[in]    a_tolerance       Tolerance
    @param[in]    a_deleteParticles If true, particles will be removed from a_activeParticles if they intersect the geometry.
  */
  template <class P>
  void
  intersectParticlesRaycastIF(ParticleContainer<P>&    a_activeParticles,
                              ParticleContainer<P>&    a_ebParticles,
                              ParticleContainer<P>&    a_domainParticles,
                              const phase::which_phase a_phase,
                              const Real               a_tolerance,
                              const bool               a_deleteParticles) const;

  /*!
    @brief Particle intersection algorithm based on bisection. 
    @details This routine will iterate through all the particles and check if they intersect the geometry. The template
    parameter indicates the particle type -- it MUST have const RealVec& position() const and  const RealVect& oldPosition() const 
    functions that determine the start and stop position of the particle trajectory. This routine uses a bisection method
    to check for intersections with the EB (the domain side is much easier). Their position are updated and they are placed
    in the a_ebParticles argument. 
    @param[inout] a_activeParticles Particles to be intersected with geometry
    @param[out]   a_ebParticles     Particles that intersected with the EB
    @param[out]   a_domainParticles Particles that intersected with the domain faces
    @param[in]    a_phase           Phase where the input particles live
    @param[in]    a_bisectionStep   Length of the bisection step
    @param[in]    a_deleteParticles If true, particles will be removed from a_activeParticles if they intersect the geometry.
  */
  template <class P>
  void
  intersectParticlesBisectIF(ParticleContainer<P>&    a_activeParticles,
                             ParticleContainer<P>&    a_ebParticles,
                             ParticleContainer<P>&    a_domainParticles,
                             const phase::which_phase a_phase,
                             const Real               a_bisectionStep,
                             const bool               a_deleteParticles) const;
  /*!
    @brief Interpolate ghost vectors over a realm, using the default ghost cell interpolation method. 
    @param[inout] a_data  Data to be interpolated.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  interpGhost(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Interpolate ghost cells over a realm, using the default ghost cell interpolation method on a specific level. 
    @param[inout] a_fineData Fine grid data
    @param[inout] a_coarData Coarse grid data
    @param[in]    a_level    The grid level corresponding to a_fineData
    @param[in]    a_realm    Realm name
    @param[in]    a_phase    Phase (gas or solid)
  */
  void
  interpGhost(LevelData<EBCellFAB>&       a_fineData,
              const LevelData<EBCellFAB>& a_coarData,
              const int                   a_level,
              const std::string           a_realm,
              const phase::which_phase    a_phase) const;

  /*!
    @brief Interpolate ghost cells over a realm. Calls the default ghost cell interpolation method. 
    @param[inout] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
  */
  void
  interpGhost(MFAMRCellData& a_data, const std::string a_realm) const;

  /*!
    @brief Interpolate ghost cells over a realm and phase. This uses piecewise linear interpolation (with limiters)
    @param[inout] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  interpGhostPwl(MFAMRCellData& a_data, const std::string a_realm) const;

  /*!
    @brief Interpolate ghost cells over a realm and phase. This uses piecewise linear interpolation (with limiters)
    @param[inout] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  interpGhostPwl(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Interpolate ghost cells over a realm and phase. 
    @details This routine uses the multigrid interpolator, which does not reach into invalid regions. I.e., all
    ghost cells are interpolated using valid data only. 
    @param[inout] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  interpGhostMG(MFAMRCellData& a_data, const std::string a_realm) const;

  /*!
    @brief Interpolate ghost cells over a realm and phase. 
    @details This routine uses the multigrid interpolator, which does not reach into invalid regions. I.e., all
    ghost cells are interpolated using valid data only. 
    @param[inout] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  interpGhostMG(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Interpolate data to new grids
    @details This is called when requiring data to be interpolated to new grids. Takes old data as argument
    and fills the new grid data with an interpolation of the old grid data. 
    @param[out] a_newData        New grid data.
    @param[in]  a_oldData        Old grid data.
    @param[in]  a_lmin           Coarsest level that did not change (but distribution may have changed). 
    @param[in]  a_oldFinestLevel Previous finest level.
    @param[in]  a_newFinestLevel New finest level. 
    @param[in]  a_useSlopes      If true, use slopes when regridding the data. 
  */
  void
  interpToNewGrids(MFAMRCellData&       a_newData,
                   const MFAMRCellData& a_oldData,
                   const int            a_lmin,
                   const int            a_oldFinestLevel,
                   const int            a_newFinestLevel,
                   const bool           a_useSlopes);

  /*!
    @brief Interpolate data to new grids
    @details This is called when requiring data to be interpolated to new grids. Takes old data as argument
    and fills the new grid data with an interpolation of the old grid data. 
    @param[out] a_newData        New grid data.
    @param[in]  a_oldData        Old grid data.
    @param[in]  a_phase          Phase on which we regrid. 
    @param[in]  a_lmin           Coarsest level that did not change (but distribution may have changed). 
    @param[in]  a_oldFinestLevel Previous finest level.
    @param[in]  a_newFinestLevel New finest level. 
    @param[in]  a_useSlopes      If true, use slopes when regridding the data. 
  */
  void
  interpToNewGrids(EBAMRCellData&           a_newData,
                   const EBAMRCellData&     a_oldData,
                   const phase::which_phase a_phase,
                   const int                a_lmin,
                   const int                a_oldFinestLevel,
                   const int                a_newFinestLevel,
                   const bool               a_useSlopes);

  /*!
    @brief Interpolate data to centroids on realm and phase
    @param[inout] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
    @details This leaves regular data intact, while the values in irregular are replaced by the corresponding interpolation to centroids. 
  */
  void
  interpToCentroids(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Sets multifluid index space. 
    @param[in] a_multiFluidIndexSpace Multifluid index space wrapper. 
  */
  void
  setMultifluidIndexSpace(const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace);

  /*!
    @brief Set implicit function for a specific phase. Need e.g. for level-sets. 
    @param[in] a_phase  Phase
    @param[in] a_baseIF The level-set function describing the phase. 
  */
  void
  setBaseImplicitFunction(const phase::which_phase a_phase, const RefCountedPtr<BaseIF>& a_baseIF);

  /*!
    @brief Regrid AMR. This versions generates the grids and Realms, but not the operator. 
    @param[in] a_tags    Cell tags which will generate the grid. 
    @param[in] a_lmin    Coarsest grid level allowed to change. 
    @param[in] a_hardcap Grid generation hardcap. If < 0 there are no limitations to grid depth. 
  */
  void
  regridAmr(const Vector<IntVectSet>& a_tags, const int a_lmin, const int a_hardcap = -1);

  /*!
    @brief Regrid a realm. This generates the grids for the realm, but does not do the operators on the realm. 
    @param[in] a_realm Realm name
    @param[in] a_procs Processor IDs
    @param[in] a_boxes Grid boxes
    @param[in] a_lmin  Coarsest level that changes. 
  */
  void
  regridRealm(const std::string          a_realm,
              const Vector<Vector<int>>& a_procs,
              const Vector<Vector<Box>>& a_boxes,
              const int                  a_lmin);

  /*!
    @brief Set grids from boxes and computational loads.
    @param[in] a_boxes          Grid boxes
    @param[in] a_realmsAndLoads Realm names and computational loads for realms
    @details This will set the grids for AmrMesh and realms, load balancing the various realms with the provided loads. 
  */
  void
  setGrids(const Vector<Vector<Box>>& a_boxes, const std::map<std::string, Vector<Vector<long int>>>& a_realmsAndLoads);

  /*!
    @brief Regrid AMR operators. This is done for all realms. 
    @param[in] a_lmin Coarsest grid level that changes. 
  */
  void
  regridOperators(const int a_lmin);

  /*!
    @brief Do a sanity check to make sure everything is set up correctly
  */
  void
  sanityCheck() const;

  /*!
    @brief Get finest grid level. 
  */
  int
  getFinestLevel() const;

  /*!
    @brief Get maximum permitted amr depth.
  */
  int
  getMaxAmrDepth() const;

  /*!
    @brief Get maximum permitted simulation depth.
  */
  int
  getMaxSimulationDepth() const;

  /*!
    @brief Get blocking factor (i.e. the smallest possible box).
  */
  int
  getBlockingFactor() const;

  /*!
    @brief Get maximum permitted box size
  */
  int
  getMaxBoxSize() const;

  /*!
    @brief Return buffer for B-R mesh refinement algorithm. 
  */
  int
  getBrBuffer() const;

  /*!
    @brief Get maximum box size for EBIS generation. 
  */
  int
  getMaxEbisBoxSize() const;

  /*!
    @brief Get the default number of ghost cells
  */
  int
  getNumberOfGhostCells() const;

  /*!
    @brief Get number of ghost cells used for EB grid generation. 
  */
  int
  getNumberOfEbGhostCells() const;

  /*!
    @brief Get default redistribution radius
  */
  int
  getRedistributionRadius() const;

  /*!
    @brief Get resolution on the finest grid level. 
  */
  Real
  getFinestDx() const;

  /*!
    @brief Get lower-left corner of computational domain. 
  */
  RealVect
  getProbLo() const;

  /*!
    @brief Get upper-left corner of computational domain. 
  */
  RealVect
  getProbHi() const;

  /*!
    @brief Get finest domain
  */
  ProblemDomain
  getFinestDomain() const;

  /*!
    @brief Get spatial resolutions
  */
  const Vector<Real>&
  getDx() const;

  /*!
    @brief Get refinement ratios
  */
  const Vector<int>&
  getRefinementRatios() const;

  /*!
    @brief Get refinement factor between two levels
    @param[in] a_level1 First level
    @param[in] a_level2 Second level
    @details This will work correctly in reverse as well.
  */
  int
  getRefinementRatio(const int a_level1, const int a_level2) const;

  /*!
    @brief Get implicit function for a specific phase
    @param[in] a_phase Phase 8gas or solid)
  */
  const RefCountedPtr<BaseIF>&
  getBaseImplicitFunction(const phase::which_phase a_phase) const;

  /*!
    @brief Get the "proxy" grids in AmrMesh.
    @details Returns m_grids
  */
  const Vector<DisjointBoxLayout>&
  getProxyGrids() const;

  /*!
    @brief Get the grids
  */
  const Vector<DisjointBoxLayout>&
  getGrids(const std::string a_realm) const;

  /*!
    @brief Get domains
  */
  const Vector<ProblemDomain>&
  getDomains() const;

  /*!
    @brief Get a registered mask.
    @param[in] a_mask   Mask name
    @param[in] a_buffer Buffer (used when registering mask)
    @param[in] a_realm  Realm name
  */
  const AMRMask&
  getMask(const std::string a_mask, const int a_buffer, const std::string a_realm) const;

  /*!
    @brief Get a map of all valid cells on a specified realm
    @param[in] a_realm Realm name
  */
  const AMRMask&
  getValidCells(const std::string a_realm) const;

  /*!
    @brief Get the EBLevelGrid for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const Vector<RefCountedPtr<EBLevelGrid>>&
  getEBLevelGrid(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get EBISLayouts for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const Vector<EBISLayout>&
  getEBISLayout(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get EBISLayouts for a Realm
    @param[in] a_realm Realm name
  */
  const Vector<RefCountedPtr<MFLevelGrid>>&
  getMFLevelGrid(const std::string a_realm) const;

  /*!
    @brief Get vof iterators for a Realm and phase. This has the capability of iterating through cut-cells. 
    @details Not const because we need to reset the vofiterator
  */
  Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
  getVofIterator(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get neighboring boxes
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex>>>>&
  getNeighbors(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get levelset function, allocated over a grid for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const EBAMRFAB&
  getLevelset(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get EBAMRParticleMesh operator
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)    
  */
  EBAMRParticleMesh&
  getParticleMesh(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get the coarsening utility
  */
  Vector<RefCountedPtr<EbCoarAve>>&
  getCoarseAverage(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get multigrid interpolation utility
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBMultigridInterpolator>>&
  getMultigridInterpolator(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get piecewise linear ghost cell interpolator
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<AggEBPWLFillPatch>>&
  getFillPatch(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get interpolator 
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBFineInterp>>&
  getFineInterp(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get flux register
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBFluxRegister>>&
  getFluxRegister(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get redistribution
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBLevelRedist>>&
  getLevelRedist(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get coar-to-fine redistribution for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBCoarToFineRedist>>&
  getCoarToFineRedist(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get coar-to-coar redistribution objects
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBCoarToCoarRedist>>&
  getCoarToCoarRedist(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get coar-to-coar redistribution redistribution objects
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBFineToCoarRedist>>&
  getFineToCoarRedist(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get interpolation stencils for going from cell center to cell centroid
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const IrregAmrStencil<CentroidInterpolationStencil>&
  getCentroidInterpolationStencils(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get interpolation stencils for going from cell center to cell boundary centroid
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const IrregAmrStencil<EbCentroidInterpolationStencil>&
  getEbCentroidInterpolationStencils(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get the nonconservative divergence stencils
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const IrregAmrStencil<NonConservativeDivergenceStencil>&
  getNonConservativeDivergenceStencils(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get the name of all Realms
    @return Names of all Realms
  */
  std::vector<std::string>
  getRealms() const;

  /*!
    @brief Get box sorting method
    @return Box sorting method.
  */
  BoxSorting
  getBoxSorting() const;

  /*!
    @brief Get EBIndexSpace corresponding to a particular phase
    @param[in] a_phase Phase (gas/solid)
  */
  RefCountedPtr<EBIndexSpace>&
  getEBIndexSpace(const phase::which_phase a_phase)
  {
    return m_multifluidIndexSpace->getEBIndexSpace(a_phase);
  }

protected:
  /*!
    @brief Enum for having understandable notation for grid generation
  */
  enum class GridGenerationMethod
  {
    BergerRigoutsous,
    Tiled,
  };

  /*!
    @brief These are all the Realms
  */
  mutable std::map<std::string, RefCountedPtr<Realm>> m_realms;

  /*!
    @brief Implicit functions
  */
  std::map<phase::which_phase, RefCountedPtr<BaseIF>> m_baseif;

  /*!
    @brief Grid generation method
  */
  GridGenerationMethod m_gridGenerationMethod;

  /*!
    @brief Box sorting
  */
  BoxSorting m_boxSort;

  /*!
    @brief MultiFluidIndexSpace
  */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /*!
    @brief Stencil type
  */
  IrregStencil::StencilType m_stencilType;

  /*!
    @brief Stencil type for centroid interpolation
  */
  IrregStencil::StencilType m_centroidStencilType;

  /*!
    @brief Stencil type for EB interpolation
  */
  IrregStencil::StencilType m_ebCentroidStencilType;

  /*!
    @brief Coarsest box where we compute
  */
  IntVect m_numCells;

  /*!
    @brief Fill ratio
  */
  Real m_fillRatioBR;

  /*!
    @brief Domain simulation corner
  */
  RealVect m_probLo;

  /*!
    @brief Domain simulation corner
  */
  RealVect m_probHi;

  /*!
    @brief Refinement ratio
  */
  int m_refRatio;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Finest level
  */
  int m_finestLevel;

  /*!
    @brief Maximum amr depth
  */
  int m_maxAmrDepth;

  /*!
    @brief Maximum allowed depth for simulation
  */
  int m_maxSimulationDepth;

  /*!
    @brief Max box size
  */
  int m_maxBoxSize;

  /*!
    @brief Maximum box size for EBIS generation
  */
  int m_maxEbisBoxSize;

  /*!
    @brief Set buffer size
  */
  int m_bufferSizeBR;

  /*!
    @brief Blocking factor
  */
  int m_blockingFactor;

  /*!
    @brief Number of ghost cells to use for eb stuff
  */
  int m_numEbGhostsCells;

  /*!
    @brief Number of ghost cells
  */
  int m_numGhostCells;

  /*!
    @brief Number of ghost cells to use when writing level-set to grid
  */
  int m_numLsfGhostCells;

  /*!
    @brief Multigrid interpolation order
  */
  int m_multigridInterpOrder;

  /*!
    @brief Multigrid interpolation radius
  */
  int m_multigridInterpRadius;

  /*!
    @brief Multigrid interpolation weights
  */
  int m_multigridInterpWeight;

  /*!
    @brief Redistribution radius
  */
  int m_redistributionRadius;

  /*!
    @brief Order for interpolation stencils to centroids
  */
  int m_centroidStencilOrder;

  /*!
    @brief Radius for centroid interpolation stencil
  */
  int m_centroidStencilRadius;

  /*!
    @brief Order for interpolation stencilst to EB
  */
  int m_ebCentroidStencilOrder;

  /*!
    @brief Radius for EB interpolation stencil
  */
  int m_ebCentroidStencilRadius;

  /*!
    @brief Stencil order for irregular stencils
  */
  int m_irregStenOrder;

  /*!
    @brief radius for irregular stencils
  */
  int m_irregStenRadius;

  /*!
    @brief Has grids or not
  */
  bool m_hasGrids;

  /*!
    @brief Grids
  */
  Vector<DisjointBoxLayout> m_grids;

  /*!
    @brief Problem domains
  */
  Vector<ProblemDomain> m_domains;

  /*!
    @brief AMR resolutions
  */
  Vector<int> m_refinementRatios;

  /*!
    @brief Level resolutions
  */
  Vector<Real> m_dx;

  /*!
    @brief Define Realms
  */
  void
  defineRealms();

  /*!
    @brief Build new internal AMR grids.
    @param[inout] a_tags Sets of cell tags used for the grid generation.
    @param[in]    a_lmin The finest grid level which changes. 
    @param[in]    a_hardcap Hardcap for the maximum grid level which can be generated. If a_hardcap < 0 there is no restriction beyond the AmrMesh restrictions.
    @details This will call the specified grid generation method and fill m_grids with the new grids. 
  */
  void
  buildGrids(Vector<IntVectSet>& a_tags, const int a_lmin, const int a_hardcap = -1);

  /*!
    @brief Parse the low/high corners of the computational domain. 
  */
  void
  parseProbLoHiCorners();

  /*!
    @brief Parse the grid generation algorithm. 
  */
  void
  parseGridGeneration();

  /*!
    @brief Parse the verbosity for AmrMesh. 
  */
  void
  parseVerbosity();

  /*!
    @brief Parse the coarsest domain grid.
  */
  void
  parseCoarsestLevelNumCells();

  /*!
    @brief Parse the maximum permitted AMR depth
  */
  void
  parseMaxAmrDepth();

  /*!
    @brief Set maximum simulation depth. 
  */
  void
  parseMaxSimulationDepth();

  /*!
    @brief Parse refinement ratios
  */
  void
  parseRefinementRatios();

  /*!
    @brief Parse buffer size for Berger-Rigoutsous grid algorithm. 
  */
  void
  parseBrBufferSize();

  /*!
    @brief Parse the Berger-Rigoutsos fill ratio. 
  */
  void
  parseBrFillRatio();

  /*!
    @brief Parse the maximum permitted box size. 
  */
  void
  parseMaxBoxSize();

  /*!
    @brief Parse the maximum permitted box size. 
    @details Only relevant for the grid generation step. 
  */
  void
  parseMaxEbisBoxSize();

  /*!
    @brief Parse the B-R blocking factor. For tiled mesh refinement this parses the tile size. 
  */
  void
  parseBlockingFactor();

  /*!
    @brief Parse number of ghost cells for eb stuff
    @details Only matters for the grid generation step. 
  */
  void
  parseEbGhostCells();

  /*!
    @brief Parse the number of ghost cells 
  */
  void
  parseNumGhostCells();

  /*!
    @brief Parse settings for the multigrid interpolator
  */
  void
  parseMultigridInterpolator();

  /*!
    @brief Parse the default redistribution radius
  */
  void
  parseRedistributionRadius();

  /*!
    @brief Parse centroid interpolation stencils
  */
  void
  parseCentroidStencils();

  /*!
    @brief Parse EB interpolation (or extrapolation) stencils
  */
  void
  parseEbCentroidStencils();

  /*!
    @brief Compute cell-centered gradient for a grid level.
    @param[out] a_gradient Cell centered gradient. 
    @param[in]  a_phi      The scalar for which the gradient is computed. 
    @param[in]  a_realm    Name of the realm where the data lives. 
    @param[in]  a_phase    Phase on which the data lives. 
    @param[in]  a_lvl      Grid level where the data lives. 
    @note This routine will reach into ghost cells and across refinement boundaries. The user must make sure that 
    ghost cells updated before using this routine. 
  */
  void
  computeGradient(LevelData<EBCellFAB>&       a_gradient,
                  const LevelData<EBCellFAB>& a_phi,
                  const std::string           a_realm,
                  const phase::which_phase    a_phase,
                  const int                   a_lvl) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_AmrMeshImplem.H>

#endif
