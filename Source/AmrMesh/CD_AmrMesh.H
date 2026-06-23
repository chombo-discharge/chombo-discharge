/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_AmrMesh.H
  @brief  Declaration of core class for handling AMR-related operations (with embedded boundaries)
  @author Robert Marskar
  @todo   a_lmax is a dead parameter in buildGrids. Fix this and the dependencies. 
  @todo   a_regsize is a dead parameter which should be removed. 
*/

#ifndef CD_AMRMESH_H
#define CD_AMRMESH_H

// Chombo includes
#include <DisjointBoxLayout.H>
#include <ProblemDomain.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_EBCoarAve.H>
#include <CD_ComputationalGeometry.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_LoadBalancing.H>
#include <CD_MFLevelGrid.H>
#include <CD_EBReflux.H>
#include <CD_DomainFluxIFFAB.H>
#include <CD_DomainFluxIFFABFactory.H>
#include <CD_ParticleContainer.H>
#include <CD_ParticleContainerSoA.H>
#include <CD_Realm.H>
#include <CD_CopyStrategy.H>
#include <CD_LoadBalancing.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Class for handling spatial operations. 
  @details This class handles most (if not all) spatial operations. Examples are grid generation, coarsening of data, interpolation of data, memory allocations etc. 
*/
class AmrMesh
{
public:
  /**
    @brief Default constructor. 
  */
  AmrMesh();

  /**
    @brief Disallowed copy constructor
    @param[in] a_other Other AmrMesh
  */
  AmrMesh(const AmrMesh& a_other) = delete;

  /**
    @brief Disallowed move constructor
    @param[in] a_other Other AmrMesh
  */
  AmrMesh(const AmrMesh&& a_other) = delete;

  /**
    @brief Destructor
  */
  virtual ~AmrMesh();

  /**
    @brief Disallowed copy assignment
    @param[in] a_other Other AmrMesh
  */
  AmrMesh&
  operator=(const AmrMesh& a_other) = delete;

  /**
    @brief Disallowed move assignment
    @param[in] a_other Other AmrMesh
  */
  AmrMesh&
  operator=(const AmrMesh&& a_other) = delete;

  /**
    @brief Method for copying from a source container to a destination container. User supplies information
    about from/to regions. 
    @param[in,out] a_dst        Destination data 
    @param[in]    a_src        Source data
    @param[in]    a_toRegion   Region we copy into
    @param[in]    a_fromRegion Region we copy from
    @note If the user requests copying into ghosted regions, the Copiers MUST USE THE CORRECT NUMBER OF GHOST CELLS
    @note This routine will not work if copying between grids before/after regrids. 
  */
  template <typename T>
  void
  copyData(EBAMRData<T>&       a_dst,
           const EBAMRData<T>& a_src,
           const CopyStrategy& a_toRegion   = CopyStrategy::Valid,
           const CopyStrategy& a_fromRegion = CopyStrategy::Valid) const noexcept;

  /**
    @brief Method for copying from a source container to a destination container. User supplies information
    @param[in,out] a_dst        Destination data 
    @param[in]    a_src        Source data
    @param[in]    a_dstComps   Destination components
    @param[in]    a_srcComps   Source components
    @param[in]    a_toRegion   Region we copy into
    @param[in]    a_fromRegion Region we copy from
    @note If the user requests copying into ghosted regions, the Copiers MUST USE THE CORRECT NUMBER OF GHOST CELLS
    @note This routine will not work if copying between grids before/after regrids. 
  */
  template <typename T>
  void
  copyData(EBAMRData<T>&       a_dst,
           const EBAMRData<T>& a_src,
           const Interval&     a_dstComps,
           const Interval&     a_srcComps,
           const CopyStrategy& a_toRegion   = CopyStrategy::Valid,
           const CopyStrategy& a_fromRegion = CopyStrategy::Valid) const noexcept;

  /**
    @brief Method for copying from a source container to a destination container. User supplies information
    about from/to regions. 
    @param[in,out] a_dst        Destination data 
    @param[in]    a_src        Source data
    @param[in]    a_level      Grid level
    @param[in]    a_toRealm    Destination realm
    @param[in]    a_fromRealm  Source realm
    @param[in]    a_fromRegion Region we copy from
    @param[in]    a_toRegion   Region we copy into
    @note If the user requests copying into ghosted regions, the Copiers MUST USE THE CORRECT NUMBER OF GHOST CELLS
    @note This routine will not work if copying between grids before/after regrids. 
  */
  template <typename T>
  void
  copyData(LevelData<T>&       a_dst,
           const LevelData<T>& a_src,
           const int           a_level,
           const std::string   a_toRealm,
           const std::string   a_fromRealm,
           const CopyStrategy& a_toRegion   = CopyStrategy::Valid,
           const CopyStrategy& a_fromRegion = CopyStrategy::Valid) const noexcept;

  /**
    @brief Method for copying from a source container to a destination container. User supplies information
    about from/to regions. 
    @param[in,out] a_dst        Destination data 
    @param[in]    a_src        Source data
    @param[in]    a_level      Grid level
    @param[in]    a_toRealm    Destination realm
    @param[in]    a_fromRealm  Source realm
    @param[in]    a_dstComps   Destination components
    @param[in]    a_srcComps   Source components
    @param[in]    a_fromRegion Region we copy from
    @param[in]    a_toRegion   Region we copy into
    @note If the user requests copying into ghosted regions, the Copiers MUST USE THE CORRECT NUMBER OF GHOST CELLS
    @note This routine will not work if copying between grids before/after regrids. 
  */
  template <typename T>
  void
  copyData(LevelData<T>&       a_dst,
           const LevelData<T>& a_src,
           const int           a_level,
           const std::string   a_toRealm,
           const std::string   a_fromRealm,
           const Interval&     a_dstComps,
           const Interval&     a_srcComps,
           const CopyStrategy& a_toRegion   = CopyStrategy::Valid,
           const CopyStrategy& a_fromRegion = CopyStrategy::Valid) const noexcept;

  /**
    @brief Deallocate data. 
    @param[in,out] a_data Data to be deleted. 
  */
  template <typename T>
  void
  deallocate(Vector<T*>& a_data) const;

  /**
    @brief Deallocate smart pointer data.
    @param[in,out] a_data Data to be deallocated. 
  */
  template <typename T>
  void
  deallocate(Vector<RefCountedPtr<T>>& a_data) const;

  /**
    @brief Deallocate data.
    @param[in,out] a_data Data to be deallocated.
  */
  template <typename T>
  void
  deallocate(EBAMRData<T>& a_data) const;

  /**
    @brief Turn smart-pointer data structure into regular-pointer data structure
    @param[out] a_alias Raw pointer aliased 
    @param[in]  a_data  Smart-pointer data
    @details This fetches the bare pointer from a_data and sets it in a_alias. 
  */
  template <typename T>
  void
  alias(Vector<T*>& a_alias, const Vector<RefCountedPtr<T>>& a_data) const;

  /**
    @brief Turn smart-pointer data structure into regular-pointer data structure
    @param[out] a_alias Raw pointer aliased 
    @param[in]  a_data  Smart-pointer data
    @details This fetches the bare pointer from a_data and sets it in a_alias. 
  */
  template <typename T, typename S>
  void
  alias(Vector<T*>& a_alias, const EBAMRData<S>& a_data) const;

  /**
    @brief Template class for generic allocation of particle data. 
    @param[out] a_particles Particle data holder. 
    @param[out] a_realm     Realm on which the particles will be allocated. 
  */
  template <typename T>
  void
  allocate(Vector<RefCountedPtr<ParticleData<T>>>& a_particles, const std::string& a_realm) const;

  /**
    @brief Template class for generic allocation of particle containers. 
    @param[out] a_container Particle container to be allocated
    @param[out] a_realm     Realm on which the particles will be allocated. 
  */
  template <typename T>
  void
  allocate(ParticleContainer<T>& a_container, const std::string& a_realm) const;

  /**
    @brief Allocate a struct-of-arrays particle container on the given realm.
    @param[out] a_container SoA particle container to be allocated.
    @param[in]  a_realm     Realm on which the particles will be allocated.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  allocate(ParticleContainerSoA<P, Traits>& a_container, const std::string& a_realm) const;

  /**
    @brief Allocate pointer but not any memory blocks. 
    @param[out] a_data  Data structure for which pointers will be constructed. 
    @details This sets all the elements in the vector to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void
  allocatePointer(Vector<RefCountedPtr<T>>& a_data) const;

  /**
    @brief Allocate pointer but not any memory blocks. 
    @param[out] a_data       Data structure for which pointers will be constructed. 
    @param[in]  a_finestLevel Finest grid level on which a_data will be defined. 
    @details This sets all the elements in the vector to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void
  allocatePointer(Vector<RefCountedPtr<T>>& a_data, const int a_finestLevel) const;

  /**
    @brief Allocate pointer but not any memory blocks. 
    @param[out] a_data Data structure for which pointers will be constructed. 
    @param[in]  a_realm Realm associated with the data
    @details This sets all the elements in a_data to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void
  allocatePointer(EBAMRData<T>& a_data, const std::string& a_realm) const;

  /**
    @brief Allocate pointer but not any memory.
    @param[out] a_data        Data structure for which pointers will be constructed. 
    @param[in]  a_realm       Realm associated with the data
    @param[in]  a_finestLevel Finest grid level on which a_data will be defined. 
    @details This sets all the elements in the vector to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void
  allocatePointer(EBAMRData<T>& a_data, const std::string& a_realm, const int a_finestLevel) const;

  /**
    @brief Slice cell-centered data in order to fetch a subset of components
    @param[in] a_variables Variables
    @return Subset of variables. The first component in the returned data holder corresponds to a_variables.begin(), the last to a_variables.end()
    @param[in] a_original Original
  */
  EBAMRCellData
  slice(EBAMRCellData& a_original, const Interval a_variables) const noexcept;

  /**
    @brief Slice cell-centered data in order to fetch a subset of components
    @param[in] a_variables Variables
    @return Subset of variables. The first component in the returned data holder corresponds to a_variables.begin(), the last to a_variables.end()
    @param[in] a_original Original
  */
  EBAMRCellData
  slice(const EBAMRCellData& a_original, const Interval a_variables) const noexcept;

  /**
    @brief Parse options. Called during the constructor. 
  */
  void
  parseOptions();

  /**
    @brief Parse runtime options
    @note Called by Driver as simulations progress. 
  */
  void
  parseRuntimeOptions();

  /**
    @brief Register an operator over a realm and a phase. 
    @param[in] a_operator Operator to be registered. 
    @param[in] a_realm    Realm name
    @param[in] a_phase    Phase (gas or solid). 
    @note Run-time errors will occur if the operator or realm do not exist. 
  */
  void
  registerOperator(const std::string& a_operator, const std::string& a_realm, const phase::which_phase a_phase);

  /**
    @brief Register a boolean mask over a realm. 
    @param[in] a_mask     Mask name to be registered. 
    @param[in] a_buffer   Relevant buffer for the mask.
    @note The meaning of a_buffer differs for different masks.
    @note Run-time errors will occur if the mask or realm do not exist.
    @param[in] a_realm Realm
  */
  void
  registerMask(const std::string& a_mask, const int a_buffer, const std::string& a_realm);

  /**
    @brief Build domains
  */
  void
  buildDomains();

  /**
    @brief Set the finest level
    @param[in] a_finestLevel
    @note This is dangerous routine, but it's needed because Driver reads checkpoint files and needs to set the finest level from those. 
  */
  void
  setFinestLevel(const int a_finestLevel);

  /**
    @brief Set the coarsest grid cells
    @param[in] a_nCells Number of cells on the coarsest grid level
  */
  void
  setCoarsestGrid(const IntVect& a_nCells);

  /**
    @brief Query if a realm exists
    @param[in] a_realm Name of the realm. 
    @return True if the realm exists and false otherwise. 
  */
  bool
  queryRealm(const std::string& a_realm) const;

  /**
    @brief Check if mesh has an EBCF
    @note Always returns true, non-EBCF cases are a thing of the past.
    @return Always true
  */
  bool
  getEbCf() const;

  /**
    @brief Register a new realm
    @param[in] a_realm Realm name
    @note It is safe to re-register a realm. 
  */
  void
  registerRealm(const std::string& a_realm);

  /**
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

  /**
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

  /**
    @brief Alias an entire hierarchy of multifluid data. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
  */
  void
  alias(EBAMRCellData& a_data, const phase::which_phase a_phase, const MFAMRCellData& a_mfdata) const;

  /**
    @brief Alias an entire hierarchy of multifluid data. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
  */
  void
  alias(EBAMRFluxData& a_data, const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata) const;

  /**
    @brief Alias an entire hierarchy of multifluid data. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
  */
  void
  alias(EBAMRIVData& a_data, const phase::which_phase a_phase, const MFAMRIVData& a_mfdata) const;

  /**
    @brief Alias an entire hierarchy of multifluid data. 
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @return A single-phase alias of the multifluid data, having the same size as the input data. 
  */
  EBAMRCellData
  alias(const phase::which_phase a_phase, const MFAMRCellData& a_mfdata) const;

  /**
    @brief Alias an entire hierarchy of multifluid data. 
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @return A single-phase alias of the multifluid data, having the same size as the input data. 
  */
  EBAMRFluxData
  alias(const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata) const;

  /**
    @brief Alias an entire hierarchy of multifluid data. 
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @return A single-phase alias of the multifluid data, having the same size as the input data. 
  */
  EBAMRIVData
  alias(const phase::which_phase a_phase, const MFAMRIVData& a_mfdata) const;

  /**
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
                  const std::string&       a_realm,
                  const phase::which_phase a_phase) const;

  /**
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
                  const std::string&       a_realm,
                  const phase::which_phase a_phase) const;

  /**
    @brief Compute cell-centered gradient over an AMR hierarchy. 
    @param[out] a_gradient Cell centered gradient. 
    @param[in]  a_phi      The scalar for which the gradient is computed. 
    @param[in]  a_realm    Name of the Realm where the data lives. 
    @note This routine will reach into ghost cells and across refinement boundaries. The user must make sure that 
    ghost cells are updated before using this routine. 
  */
  void
  computeGradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi, const std::string& a_realm) const;

  /**
    @brief Compute face-centered gradient over an AMR hierarchy. 
    @param[out] a_gradient Face centered gradient. 
    @param[in]  a_phi      The scalar for which the gradient is computed. 
    @param[in]  a_realm    Name of the Realm where the data lives. 
    @note This routine will reach into ghost cells and across refinement boundaries. The user must make sure that 
    ghost cells are updated before using this routine. 
  */
  void
  computeGradient(MFAMRFluxData& a_gradient, const MFAMRCellData& a_phi, const std::string& a_realm) const;

  /**
    @brief Allocate a data holder over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_phase  Phase (gas or solid)
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_ghost  Number of ghost cells for a_data
    @note If a_ghost < 0, this routine will use the default number of ghost cells.
  */
  void
  allocate(EBAMRCellData&           a_data,
           const std::string&       a_realm,
           const phase::which_phase a_phase,
           const int                a_nComp,
           const int                a_ghost = -1) const;

  /**
    @brief Allocate data on specific level
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_phase  Phase (gas or solid)
    @param[in]  a_level  AMR level
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_ghost  Number of ghost cells for a_data
    @note If a_ghost < 0, this routine will use the default number of ghost cells.
  */
  void
  allocate(LevelData<EBCellFAB>&    a_data,
           const std::string&       a_realm,
           const phase::which_phase a_phase,
           const int                a_level,
           const int                a_nComp,
           const int                a_ghost = -1) const;

  /**
    @brief Allocate a data holder over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_phase  Phase (gas or solid)
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_ghost  Number of ghost cells for a_data
    @note If a_ghost < 0, this routine will use the default number of ghost cells.
  */
  void
  allocate(EBAMRFluxData&           a_data,
           const std::string&       a_realm,
           const phase::which_phase a_phase,
           const int                a_nComp,
           const int                a_ghost = -1) const;

  /**
    @brief Allocate a data holder over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_phase  Phase (gas or solid)
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_ghost  Number of ghost cells for a_data
    @note If a_ghost < 0, this routine will use the default number of ghost cells.
  */
  void
  allocate(EBAMRIVData&             a_data,
           const std::string&       a_realm,
           const phase::which_phase a_phase,
           const int                a_nComp,
           const int                a_ghost = -1) const;

  /**
    @brief Allocate a data holder over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated.
    @param[in]  a_phase  Phase (gas or solid).
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_ghost  Number of ghost cells for a_data
    @note If a_ghost < 0, this routine will use the default number of ghost cells.
  */
  void
  allocate(EBAMRIFData&             a_data,
           const std::string&       a_realm,
           const phase::which_phase a_phase,
           const int                a_nComp,
           const int                a_ghost = -1) const;

  /**
    @brief Allocate Boolean data over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated.
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_ghost  Number of ghost cells for a_data
    @note If a_ghost < 0, this routine will use the default number of ghost cells.
  */
  void
  allocate(EBAMRBool& a_data, const std::string& a_realm, const int a_nComp, const int a_ghost = -1) const;

  /**
    @brief Allocate a data holder over a specific realm.
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated.
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_ghost  Number of ghost cells for a_data
    @note If a_ghost < 0, this routine will use the default number of ghost cells.
  */
  void
  allocate(MFAMRCellData& a_data, const std::string& a_realm, const int a_nComp, const int a_ghost = -1) const;

  /**
    @brief Allocate a data holder over a specific realm.
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated.
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_ghost  Number of ghost cells for a_data
    @note If a_ghost < 0, this routine will use the default number of ghost cells.
  */
  void
  allocate(MFAMRFluxData& a_data, const std::string& a_realm, const int a_nComp, const int a_ghost = -1) const;

  /**
    @brief Allocate a data holder over a specific realm
    @param[out] a_data   Data holder to be allocated
    @param[in]  a_realm  Realm of the name where the data will be allocated.
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_ghost  Number of ghost cells for a_data
    @note If a_ghost < 0, this routine will use the default number of ghost cells.
  */
  void
  allocate(MFAMRIVData& a_data, const std::string& a_realm, const int a_nComp, const int a_ghost = 0) const;

  /**
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_lmin) const;

  /**
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_lmin) const;

  /**
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_lmin) const;

  /**
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_lmin) const;

  /**
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(EBAMRBool& a_data, const int a_lmin) const;

  /**
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(MFAMRCellData& a_data, const int a_lmin) const;

  /**
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(MFAMRFluxData& a_data, const int a_lmin) const;

  /**
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void
  reallocate(MFAMRIVData& a_data, const int a_lmin) const;

  /**
    @brief Average multifluid data over a specified realm
    @param[in,out] a_data    Data to be coarsened.
    @param[in]    a_realm   Realm name
    @param[in]    a_average Averaging method
  */
  void
  average(MFAMRCellData& a_data, const std::string& a_realm, const Average& a_average) const;

  /**
    @brief Arithmetic coarsening of multifluid data.
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
  */
  void
  arithmeticAverage(MFAMRCellData& a_data, const std::string& a_realm) const;

  /**
    @brief Harmonic coarsening of multifluid data.
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
  */
  void
  harmonicAverage(MFAMRCellData& a_data, const std::string& a_realm) const;

  /**
    @brief Conservative coarsening of multifluid data.
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
  */
  void
  conservativeAverage(MFAMRCellData& a_data, const std::string& a_realm) const;

  /**
    @brief Average multifluid data over a specified realm
    @param[in,out] a_data    Data to be coarsened.
    @param[in]    a_realm   Realm name
    @param[in]    a_average Averaging method
  */
  void
  average(MFAMRFluxData& a_data, const std::string& a_realm, const Average& a_average) const;

  /**
    @brief Average multifluid data over a realm. 
    @details This computes an arithmetic average of the face data. Data on coarse faces is replaced by the arithmetic average of the fine face data. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
  */
  void
  arithmeticAverage(MFAMRFluxData& a_data, const std::string& a_realm) const;

  /**
    @brief Average multifluid data over a realm. 
    @details This computes an arithmetic average of the face data. Data on coarse faces is replaced by the arithmetic average of the fine face data. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
  */
  void
  harmonicAverage(MFAMRFluxData& a_data, const std::string& a_realm) const;

  /**
    @brief Average down multifluid data over a realm. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
  */
  void
  conservativeAverage(MFAMRFluxData& a_data, const std::string& a_realm) const;

  /**
    @brief Average down on specific realm and phase. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
    @param[in]    a_average Averaging method
  */
  void
  average(EBAMRCellData&           a_data,
          const std::string&       a_realm,
          const phase::which_phase a_phase,
          const Average&           a_average) const;

  /**
    @brief Average down on specific realm and phase. 
    @param[in,out] a_data      Data to be coarsened.
    @param[in]    a_realm     Realm
    @param[in]    a_phase     Phase (gas or solid)
    @param[in]    a_variables Variables to average
    @param[in]    a_average   Averaging method
  */
  void
  average(EBAMRCellData&           a_data,
          const std::string&       a_realm,
          const phase::which_phase a_phase,
          const Average&           a_average,
          const Interval&          a_variables) const;

  /**
    @brief Arithmetic average of data. Does all components. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  arithmeticAverage(EBAMRCellData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Arithmetic average of data. 
    @param[in,out] a_data      Data to be coarsened.
    @param[in]    a_realm     Realm name
    @param[in]    a_phase     Phase (gas or solid)
    @param[in]    a_variables Variables
  */
  void
  arithmeticAverage(EBAMRCellData&           a_data,
                    const std::string&       a_realm,
                    const phase::which_phase a_phase,
                    const Interval&          a_variables) const;

  /**
    @brief Harmonic average of data. Does all components. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  harmonicAverage(EBAMRCellData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Harmonic average of data. Does all components. 
    @param[in,out] a_data      Data to be coarsened.
    @param[in]    a_realm     Realm name
    @param[in]    a_phase     Phase (gas or solid)
    @param[in]    a_variables Variables
  */
  void
  harmonicAverage(EBAMRCellData&           a_data,
                  const std::string&       a_realm,
                  const phase::which_phase a_phase,
                  const Interval&          a_variables) const;

  /**
    @brief Conseratively average data. Does all components. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  conservativeAverage(EBAMRCellData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Conseratively average data. Does all components. 
    @param[in,out] a_data      Data to be coarsened.
    @param[in]    a_realm     Realm name
    @param[in]    a_phase     Phase (gas or solid)
    @param[in]    a_variables Variables
  */
  void
  conservativeAverage(EBAMRCellData&           a_data,
                      const std::string&       a_realm,
                      const phase::which_phase a_phase,
                      const Interval&          a_variables) const;

  /**
    @brief Average down on specific realm and phase. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
    @param[in]    a_average Averaging method
  */
  void
  average(EBAMRFluxData&           a_data,
          const std::string&       a_realm,
          const phase::which_phase a_phase,
          const Average&           a_average) const;

  /**
    @brief Average down on specific realm and phase. 
    @details This computes an arithmetic average of the face data. Data on coarse faces is replaced by the arithmetic average of the fine face data. 
    @param[in,out] a_data  Data to be average.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  arithmeticAverage(EBAMRFluxData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Average down on specific realm and phase. 
    @details This computes an arithmetic average of the face data. Data on coarse faces is replaced by the arithmetic average of the fine face data. 
    @param[in,out] a_data  Data to be average.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  harmonicAverage(EBAMRFluxData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Average down on specific realm and phase. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  conservativeAverage(EBAMRFluxData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Average down on specific realm and phase. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
    @param[in]    a_average Averaging method
  */
  void
  average(EBAMRIVData&             a_data,
          const std::string&       a_realm,
          const phase::which_phase a_phase,
          const Average&           a_average) const;

  /**
    @brief Arithmetic average on specific realm and phase. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  arithmeticAverage(EBAMRIVData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Harmonic average on specific realm and phase. 
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  harmonicAverage(EBAMRIVData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Conservative averaging on specific realm and phase.
    @param[in,out] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  conservativeAverage(EBAMRIVData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Deposit scalar particle quantities on the mesh.
    @details P is the particle type, Ret is the returned value of the particle member function, and must be a Real or a
    RealVect. MemFunc is a pointer to a member function of P. E.g., depositParticles<P, const RealVect&, &P::position>.    
    @param[out] a_meshData             Mesh data. Must have exactly one component. 
    @param[in]  a_realm                Realm where data is registered.
    @param[in]  a_phase                Phase where data is registered.
    @param[in]  a_particles            Particle container. Must be in "usable state" for deposition.
    @param[in]  a_depositionType       Specification of deposition kernel (e.g., CIC)
    @param[in]  a_coarseFineDeposition Specification of handling of coarse-fine boundaries.
    @param[in]  a_forceIrregNGP        Force NGP deposition in irregular cells or not. 
  */
  template <class P, typename Ret, Ret (P::*MemberFunc)() const>
  void
  depositParticles(EBAMRCellData&              a_meshData,
                   const std::string&          a_realm,
                   const phase::which_phase&   a_phase,
                   const ParticleContainer<P>& a_particles,
                   const DepositionType        a_depositionType,
                   const CoarseFineDeposition  a_coarseFineDeposition,
                   const bool                  a_forceIrregNGP = false);

  /**
    @brief Deposit scalar particle quantities on the mesh.
    @details P is the particle type, Ret is the returned value of the particle member function, and must be a Real or a
    const Real&. MemFunc is a pointer to a member function of P. E.g., depositParticles<P, const Real&, &P::weight>.        
    @param[out] a_meshData   Mesh data. MUST have exactly one component. 
    @param[in]  a_realm      Realm where data is registered.
    @param[in]  a_phase      Phase where data is registered.
    @param[in]  a_particles  Particle container. Must be in "usable state" for deposition.
  */
  template <class P, typename Ret, Ret (P::*MemberFunc)() const>
  void
  depositParticles(EBAMRIVData&                a_meshData,
                   const std::string&          a_realm,
                   const phase::which_phase&   a_phase,
                   const ParticleContainer<P>& a_particles) const noexcept;

  /**
    @brief Interpolate mesh data onto a particle position.
    @details P is the particle type, Ret is the returned value of the particle member function, and must be a Real& or a
    RealVect&. MemFunc is a pointer to a member function of P. E.g., interpolate<P, Real&, &P::weight>.        
    @param[in,out] a_particles            Particles to be interpolated. 
    @param[in]    a_realm                Realm where data is registered.
    @param[in]    a_phase                Phase where data is registered.
    @param[in]    a_meshScalarField      Scalar field on the mesh 
    @param[in]    a_interpType           Interpolation type. 
    @param[in]    a_forceIrregNGP        Force NGP interpolation in cut-cells. 
  */
  template <class P, class Ret, Ret (P::*MemberFunc)()>
  void
  interpolateParticles(ParticleContainer<P>&     a_particles,
                       const std::string&        a_realm,
                       const phase::which_phase& a_phase,
                       const EBAMRCellData&      a_meshScalarField,
                       const DepositionType      a_interpType,
                       const bool                a_forceIrregNGP = false) const;

  /**
    @brief Deposit the SoA container's weight column on the mesh (all coarse-fine strategies).
    @param[out] a_meshData             Mesh data. Must have exactly one component.
    @param[in]  a_realm                Realm where data is registered.
    @param[in]  a_phase                Phase where data is registered.
    @param[in]  a_particles            SoA particle container.
    @param[in]  a_depositionType       Deposition kernel (e.g. CIC).
    @param[in]  a_coarseFineDeposition Coarse-fine handling.
    @param[in]  a_forceIrregNGP        Force NGP deposition in cut cells.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  depositWeight(EBAMRCellData&                         a_meshData,
                const std::string&                     a_realm,
                const phase::which_phase&              a_phase,
                const ParticleContainerSoA<P, Traits>& a_particles,
                const DepositionType                   a_depositionType,
                const CoarseFineDeposition             a_coarseFineDeposition,
                const bool                             a_forceIrregNGP);

  /**
    @brief Deposit a custom per-particle scalar (a_cellGather) of an SoA container on the mesh.
    @details Generalizes depositWeight to a derived per-particle quantity (e.g. weight*mobility), via the
    EBAMRParticleMeshSoA gatherer-based deposit. Deposits one component.
    @param[out] a_meshData             Mesh data (one component).
    @param[in]  a_realm                Realm where data is registered.
    @param[in]  a_phase                Phase where data is registered.
    @param[in]  a_particles            SoA particle container.
    @param[in]  a_depositionType       Deposition kernel (e.g. CIC).
    @param[in]  a_coarseFineDeposition Coarse-fine handling.
    @param[in]  a_forceIrregNGP        Force NGP deposition in cut cells.
    @param[in]  a_cellGather           Per-particle value gatherer (leaf, index) -> Real.
    @tparam P          Payload type of the container.
    @tparam Traits     Column descriptor for the payload.
    @tparam CellGather Callable (const ParticleSoA<P, Traits>&, std::size_t) -> Real.
  */
  template <typename P, typename Traits, typename CellGather>
  void
  depositGathered(EBAMRCellData&                         a_meshData,
                  const std::string&                     a_realm,
                  const phase::which_phase&              a_phase,
                  const ParticleContainerSoA<P, Traits>& a_particles,
                  const DepositionType                   a_depositionType,
                  const CoarseFineDeposition             a_coarseFineDeposition,
                  const bool                             a_forceIrregNGP,
                  CellGather                             a_cellGather);

  /**
    @brief Deposit one or more SoA payload columns on the mesh (all coarse-fine strategies).
    @details Members selects either one scalar column or SpaceDim vector-component columns, mapped to
    consecutive mesh components.
    @param[out] a_meshData             Mesh data.
    @param[in]  a_realm                Realm where data is registered.
    @param[in]  a_phase                Phase where data is registered.
    @param[in]  a_particles            SoA particle container.
    @param[in]  a_depositionType       Deposition kernel.
    @param[in]  a_coarseFineDeposition Coarse-fine handling.
    @param[in]  a_forceIrregNGP        Force NGP deposition in cut cells.
    @tparam Members Pack of payload data-member pointers selecting the column(s).
    @tparam P       Payload type (deduced).
    @tparam Traits  Column descriptor (deduced).
  */
  template <auto... Members, typename P, typename Traits>
  void
  depositParticles(EBAMRCellData&                         a_meshData,
                   const std::string&                     a_realm,
                   const phase::which_phase&              a_phase,
                   const ParticleContainerSoA<P, Traits>& a_particles,
                   const DepositionType                   a_depositionType,
                   const CoarseFineDeposition             a_coarseFineDeposition,
                   const bool                             a_forceIrregNGP);

  /**
    @brief Interpolate a mesh scalar field onto the SoA container's weight column.
    @param[in,out] a_particles     SoA particles whose weight column is overwritten.
    @param[in]     a_realm         Realm where data is registered.
    @param[in]     a_phase         Phase where data is registered.
    @param[in]     a_meshScalarField Scalar field on the mesh.
    @param[in]     a_interpType    Interpolation kernel.
    @param[in]     a_forceIrregNGP Force NGP interpolation in cut cells.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  interpolateWeight(ParticleContainerSoA<P, Traits>& a_particles,
                    const std::string&               a_realm,
                    const phase::which_phase&        a_phase,
                    const EBAMRCellData&             a_meshScalarField,
                    const DepositionType             a_interpType,
                    const bool                       a_forceIrregNGP) const;

  /**
    @brief Interpolate a mesh field onto one or more SoA payload columns.
    @details Members selects either one scalar column or SpaceDim vector-component columns, read from
    consecutive mesh components.
    @param[in,out] a_particles     SoA particles whose target columns are overwritten.
    @param[in]     a_realm         Realm where data is registered.
    @param[in]     a_phase         Phase where data is registered.
    @param[in]     a_meshField     Field on the mesh.
    @param[in]     a_interpType    Interpolation kernel.
    @param[in]     a_forceIrregNGP Force NGP interpolation in cut cells.
    @tparam Members Pack of payload data-member pointers selecting the target column(s).
    @tparam P       Payload type (deduced).
    @tparam Traits  Column descriptor (deduced).
  */
  template <auto... Members, typename P, typename Traits>
  void
  interpolateParticles(ParticleContainerSoA<P, Traits>& a_particles,
                       const std::string&               a_realm,
                       const phase::which_phase&        a_phase,
                       const EBAMRCellData&             a_meshField,
                       const DepositionType             a_interpType,
                       const bool                       a_forceIrregNGP) const;

  /**
    @brief Function which removes particles from the domain if they fall inside the EB.
    @details The template parameter indicate the particle type -- it MUST have a const RealVect& P::position() const function. 
    @details This version uses the implicit function to evaluate whether or not the particles are inside the EB. The particles will be removed from
    the container if f(x) > a_tolerance*dx. 
    @param[in,out] a_particles Particle data
    @param[in]    a_phase     Phase where the particles live. 
    @param[in]    a_tolerance Tolerance. 
    @note Because the implicit functions are global, this function will work even if the particles aren't mapped to the correct EBISBox.
  */
  template <class P>
  void
  removeCoveredParticlesIF(ParticleContainer<P>&     a_particles,
                           const phase::which_phase& a_phase,
                           const Real                a_tolerance = 0.0) const;

  /**
    @brief Remove SoA particles that fall inside the EB (implicit-function test): f(x) > a_tolerance*dx.
    @param[in,out] a_particles SoA particle data.
    @param[in]     a_phase     Phase where the particles live.
    @param[in]     a_tolerance Tolerance (multiplied by dx).
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  removeCoveredParticlesIF(ParticleContainerSoA<P, Traits>& a_particles,
                           const phase::which_phase&        a_phase,
                           const Real                       a_tolerance) const;

  /**
    @brief Function which removes particles from the domain if they fall inside the EB.
    @details The template parameter indicate the particle type -- it MUST have a const RealVect& P::position() const function. 
    @details This version uses the discrete information to evaluate whether or not the particles are inside the EB. The particles will be removed from
    the container if they lie within a covered cell OR if they lie in an irregular cell but on the "wrong" side of the EB face (to within precision tolerance*dx)
    @param[in,out] a_particles Particle data
    @param[in]    a_phase     Phase where the particles live. 
    @param[in]    a_tolerance Tolerance when to remove particles
    @note This function will only work if the particles are mapped to the correct EBISBox. 
  */
  template <class P>
  void
  removeCoveredParticlesDiscrete(ParticleContainer<P>&     a_particles,
                                 const phase::which_phase& a_phase,
                                 const Real                a_tolerance = 0.0) const;

  /**
    @brief Remove SoA particles inside the EB using discrete (EBISBox) information: covered cells, or
           irregular cells on the wrong side of every VoF's EB face (to within a_tolerance*dx).
    @param[in,out] a_particles SoA particle data.
    @param[in]     a_phase     Phase where the particles live.
    @param[in]     a_tolerance Tolerance (multiplied by dx).
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  removeCoveredParticlesDiscrete(ParticleContainerSoA<P, Traits>& a_particles,
                                 const phase::which_phase&        a_phase,
                                 const Real                       a_tolerance) const;

  /**
    @brief Function which removes particles from the domain if they fall inside the EB.
    @details The template parameter indicate the particle type -- it must have a const RealVect& P::position() function.
    @details This version removes all particles that live in covered cells.
    @param[in,out] a_particles Particle data
    @param[in]    a_phase     Phase where the particles live.
    @note This function will only work if the particles are mapped to the correct EBISBox.
  */
  template <class P>
  void
  removeCoveredParticlesVoxels(ParticleContainer<P>& a_particles, const phase::which_phase& a_phase) const;

  /**
    @brief Remove SoA particles that live in covered cells (voxel test).
    @param[in,out] a_particles SoA particle data.
    @param[in]     a_phase     Phase where the particles live.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  removeCoveredParticlesVoxels(ParticleContainerSoA<P, Traits>& a_particles, const phase::which_phase& a_phase) const;

  /**
    @brief Function which transferse particles from one particle container to another if they fall inside the EB.
    @details The template parameter indicate the particle type -- it MUST have a const RealVect& P::position() const function. 
    @details This version uses the implicit function to evaluate whether or not the particles are inside the EB. The particles will be transferred from
    the container if f(x) > a_tolerance*dx.
    @param[in,out] a_particlesFrom Container to transfer from
    @param[in,out] a_particlesTo   Container to transfer to
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

  /**
    @brief Function which transferse particles from one particle container to another if they fall inside the EB.
    @details The template parameter indicate the particle type -- it MUST have a const RealVect& P::position() const function. 
    @details This version uses the discrete information to evaluate whether or not the particles are inside the EB. The particles will be transferred from
    the container if they lie within a covered cell OR if they lie in an irregular cell but on the "wrong" side of the EB face (to within precision tolerance*dx)
    @param[in,out] a_particlesFrom Container to transfer from
    @param[in,out] a_particlesTo   Container to transfer to
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

  /**
    @brief Function which transferse particles from one particle container to another if they fall inside the EB.
    @details The template parameter indicate the particle type -- it MUST have a const RealVect& P::position() const function. 
    @details This version transfers all particles that live in covered cells.
    @param[in,out] a_particlesFrom Container to transfer from
    @param[in,out] a_particlesTo   Container to transfer to
    @param[in]    a_phase         Phase where the particles live. 
    @note This function will only work if the particles are mapped to the correct EBISBox. 
  */
  template <class P>
  void
  transferCoveredParticlesVoxels(ParticleContainer<P>&     a_particlesFrom,
                                 ParticleContainer<P>&     a_particlesTo,
                                 const phase::which_phase& a_phase) const;

  /**
    @brief Transfer SoA particles inside the EB (implicit-function test) from one container to another: f(x) > a_tolerance*dx.
    @param[in,out] a_particlesFrom Source SoA container (transferred particles are removed from it).
    @param[in,out] a_particlesTo   Destination SoA container.
    @param[in]     a_phase         Phase where the particles live.
    @param[in]     a_tolerance     Tolerance (multiplied by dx).
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
    @note Because the implicit functions are global, this works even if the particles aren't mapped to the correct boxes (remapping will be necessary afterwards).
  */
  template <typename P, typename Traits>
  void
  transferCoveredParticlesIF(ParticleContainerSoA<P, Traits>& a_particlesFrom,
                             ParticleContainerSoA<P, Traits>& a_particlesTo,
                             const phase::which_phase&        a_phase,
                             const Real                       a_tolerance) const;

  /**
    @brief Transfer SoA particles inside the EB using discrete (EBISBox) information from one container to another.
    @details Transfers particles in covered cells, or in irregular cells on the wrong side of every VoF's EB face (to within a_tolerance*dx).
    @param[in,out] a_particlesFrom Source SoA container (transferred particles are removed from it).
    @param[in,out] a_particlesTo   Destination SoA container.
    @param[in]     a_phase         Phase where the particles live.
    @param[in]     a_tolerance     Tolerance (multiplied by dx).
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  transferCoveredParticlesDiscrete(ParticleContainerSoA<P, Traits>& a_particlesFrom,
                                   ParticleContainerSoA<P, Traits>& a_particlesTo,
                                   const phase::which_phase&        a_phase,
                                   const Real                       a_tolerance) const;

  /**
    @brief Transfer SoA particles that live in covered cells (voxel test) from one container to another.
    @param[in,out] a_particlesFrom Source SoA container (transferred particles are removed from it).
    @param[in,out] a_particlesTo   Destination SoA container.
    @param[in]     a_phase         Phase where the particles live.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  transferCoveredParticlesVoxels(ParticleContainerSoA<P, Traits>& a_particlesFrom,
                                 ParticleContainerSoA<P, Traits>& a_particlesTo,
                                 const phase::which_phase&        a_phase) const;

  /**
    @brief Transfer particles that are on the wrong side of the EB to a different container.
    @param[in,out] a_dstParticles     Destination particle container
    @param[in,out] a_srcParticles     Source particle container
    @param[in]    a_phase            Phase where the particles should live
    @param[in]    a_transferModifier Optional input argument for letting the user manipulate particles that were transferred
  */
  template <class P>
  void
  transferIrregularParticles(
    ParticleContainer<P>&          a_dstParticles,
    ParticleContainer<P>&          a_srcParticles,
    const phase::which_phase       a_phase,
    const std::function<void(P&)>& a_transferModifier = [](P&) -> void {
      return;
    }) const noexcept;

  /**
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
    @param[in,out] a_activeParticles     Particles to be intersected with geometry
    @param[out]   a_ebParticles         Particles that intersected with the EB
    @param[out]   a_domainParticles     Particles that intersected with the domain faces
    @param[in]    a_phase               Phase where the input particles live
    @param[in]    a_deleteParticles     If true, particles will be removed from a_activeParticles if they intersect the geometry.
    @param[in]    a_nonDeletionModifier Optional input argument for letting the user manipulate particles that were intersected but not deleted
    @param[in] a_tolerance Tolerance
  */
  template <class P>
  void
  intersectParticlesRaycastIF(
    ParticleContainer<P>&          a_activeParticles,
    ParticleContainer<P>&          a_ebParticles,
    ParticleContainer<P>&          a_domainParticles,
    const phase::which_phase       a_phase,
    const Real                     a_tolerance,
    const bool                     a_deleteParticles,
    const std::function<void(P&)>& a_nonDeletionModifier = [](P&) -> void {
      return;
    }) const noexcept;

  /**
    @brief Particle intersection algorithm based on bisection. 
    @details This routine will iterate through all the particles and check if they intersect the geometry. The template
    parameter indicates the particle type -- it MUST have const RealVec& position() const and  const RealVect& oldPosition() const 
    functions that determine the start and stop position of the particle trajectory. This routine uses a bisection method
    to check for intersections with the EB (the domain side is much easier). Their position are updated and they are placed
    in the a_ebParticles argument. 
    @param[in,out] a_activeParticles     Particles to be intersected with geometry
    @param[out]   a_ebParticles         Particles that intersected with the EB
    @param[out]   a_domainParticles     Particles that intersected with the domain faces
    @param[in]    a_phase               Phase where the input particles live
    @param[in]    a_bisectionStep       Length of the bisection step
    @param[in]    a_deleteParticles     If true, particles will be removed from a_activeParticles if they intersect the geometry.
    @param[in]    a_nonDeletionModifier Optional input argument for letting the user manipulate particles that were intersected but not deleted
  */
  template <class P>
  void
  intersectParticlesBisectIF(
    ParticleContainer<P>&          a_activeParticles,
    ParticleContainer<P>&          a_ebParticles,
    ParticleContainer<P>&          a_domainParticles,
    const phase::which_phase       a_phase,
    const Real                     a_bisectionStep,
    const bool                     a_deleteParticles,
    const std::function<void(P&)>& a_nonDeletionModifier = [](P&) -> void {
      return;
    }) const noexcept;

  /**
    @brief SoA ray-casting particle intersection algorithm.
    @details Iterates through the active particles and checks if they intersect the EB or the domain faces along
    their trajectory (oldPosition -> position). Particles that intersect are moved onto the intersection point and
    placed in a_ebParticles or a_domainParticles. The EB test uses the ray-casting method (assuming the implicit
    function is a signed distance function). The container-owned position holds the trajectory end point; the start
    point is read from the payload columns selected by the OldPosition member-pointer pack.
    @param[in,out] a_activeParticles Particles to be intersected with geometry.
    @param[out]   a_ebParticles     Particles that intersected the EB (cleared and remapped here).
    @param[out]   a_domainParticles Particles that intersected the domain faces (cleared and remapped here).
    @param[in]    a_phase           Phase where the input particles live.
    @param[in]    a_tolerance       Tolerance.
    @param[in]    a_deleteParticles If true, intersected particles are removed from a_activeParticles.
    @param[in]    a_nonDeletionModifier Callback applied to each intersected-but-not-deleted particle, given its
    home leaf and slot index (used e.g. to flag the original particle via a payload column).
    @tparam OldPosition Pack of SpaceDim payload member pointers giving the particle's previous-position components.
    @tparam P           Payload type (deduced).
    @tparam Traits      Column descriptor (deduced).
  */
  template <auto... OldPosition, typename P, typename Traits>
  void
  intersectParticlesRaycastIF(
    ParticleContainerSoA<P, Traits>&                                 a_activeParticles,
    ParticleContainerSoA<P, Traits>&                                 a_ebParticles,
    ParticleContainerSoA<P, Traits>&                                 a_domainParticles,
    const phase::which_phase                                         a_phase,
    const Real                                                       a_tolerance,
    const bool                                                       a_deleteParticles,
    const std::function<void(ParticleSoA<P, Traits>&, std::size_t)>& a_nonDeletionModifier) const noexcept;

  /**
    @brief SoA bisection-based particle intersection algorithm.
    @details Identical to intersectParticlesRaycastIF but the EB intersection test uses bisection rather than
    ray-casting. The container-owned position holds the trajectory end point; the start point is read from the
    payload columns selected by the OldPosition member-pointer pack.
    @param[in,out] a_activeParticles Particles to be intersected with geometry.
    @param[out]   a_ebParticles     Particles that intersected the EB (cleared and remapped here).
    @param[out]   a_domainParticles Particles that intersected the domain faces (cleared and remapped here).
    @param[in]    a_phase           Phase where the input particles live.
    @param[in]    a_bisectionStep   Length of the bisection step.
    @param[in]    a_deleteParticles If true, intersected particles are removed from a_activeParticles.
    @param[in]    a_nonDeletionModifier Callback applied to each intersected-but-not-deleted particle, given its
    home leaf and slot index (used e.g. to flag the original particle via a payload column).
    @tparam OldPosition Pack of SpaceDim payload member pointers giving the particle's previous-position components.
    @tparam P           Payload type (deduced).
    @tparam Traits      Column descriptor (deduced).
  */
  template <auto... OldPosition, typename P, typename Traits>
  void
  intersectParticlesBisectIF(
    ParticleContainerSoA<P, Traits>&                                 a_activeParticles,
    ParticleContainerSoA<P, Traits>&                                 a_ebParticles,
    ParticleContainerSoA<P, Traits>&                                 a_domainParticles,
    const phase::which_phase                                         a_phase,
    const Real                                                       a_bisectionStep,
    const bool                                                       a_deleteParticles,
    const std::function<void(ParticleSoA<P, Traits>&, std::size_t)>& a_nonDeletionModifier) const noexcept;

  /**
    @brief Interpolate ghost vectors over a realm, using the default ghost cell interpolation method. 
    @param[in,out] a_data  Data to be interpolated.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  interpGhost(EBAMRCellData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Interpolate ghost cells over a realm, using the default ghost cell interpolation method on a specific level. 
    @param[in,out] a_fineData  Fine grid data
    @param[in,out] a_coarData  Coarse grid data
    @param[in]    a_fineLevel The grid level corresponding to a_fineData
    @param[in]    a_realm     Realm name
    @param[in]    a_phase     Phase (gas or solid)
  */
  void
  interpGhost(LevelData<EBCellFAB>&       a_fineData,
              const LevelData<EBCellFAB>& a_coarData,
              const int                   a_fineLevel,
              const std::string&          a_realm,
              const phase::which_phase    a_phase) const;

  /**
    @brief Interpolate ghost cells over a realm. Calls the default ghost cell interpolation method. 
    @param[in,out] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
  */
  void
  interpGhost(MFAMRCellData& a_data, const std::string& a_realm) const;

  /**
    @brief Interpolate ghost cells over a realm and phase. This uses piecewise linear interpolation (with limiters)
    @param[in,out] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
  */
  void
  interpGhostPwl(MFAMRCellData& a_data, const std::string& a_realm) const;

  /**
    @brief Interpolate ghost cells over a realm and phase. This uses piecewise linear interpolation (with limiters)
    @param[in,out] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  interpGhostPwl(EBAMRCellData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Interpolate ghost cells over a realm and phase. 
    @details This routine uses the multigrid interpolator, which does not reach into invalid regions. I.e., all
    ghost cells are interpolated using valid data only. 
    @param[in,out] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
  */
  void
  interpGhostMG(MFAMRCellData& a_data, const std::string& a_realm) const;

  /**
    @brief Interpolate ghost cells over a realm and phase. 
    @details This routine uses the multigrid interpolator, which does not reach into invalid regions. I.e., all
    ghost cells are interpolated using valid data only. 
    @param[in,out] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void
  interpGhostMG(EBAMRCellData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Interpolate data to new grids
    @details This is called when requiring data to be interpolated to new grids. Takes old data as argument
    and fills the new grid data with an interpolation of the old grid data. 
    @param[out] a_newData        New grid data.
    @param[in]  a_oldData        Old grid data.
    @param[in]  a_lmin           Coarsest level that did not change (but distribution may have changed). 
    @param[in]  a_oldFinestLevel Previous finest level.
    @param[in]  a_newFinestLevel New finest level. 
    @param[in]  a_type Interpolation type
  */
  void
  interpToNewGrids(MFAMRCellData&                   a_newData,
                   const MFAMRCellData&             a_oldData,
                   const int                        a_lmin,
                   const int                        a_oldFinestLevel,
                   const int                        a_newFinestLevel,
                   const EBCoarseToFineInterp::Type a_type);

  /**
    @brief Interpolate data to new grids
    @details This is called when requiring data to be interpolated to new grids. Takes old data as argument
    and fills the new grid data with an interpolation of the old grid data. 
    @param[out] a_newData        New grid data.
    @param[in]  a_oldData        Old grid data.
    @param[in]  a_phase          Phase on which we regrid. 
    @param[in]  a_lmin           Coarsest level that did not change (but distribution may have changed). 
    @param[in]  a_oldFinestLevel Previous finest level.
    @param[in]  a_newFinestLevel New finest level. 
    @param[in]  a_type Interpolation type
  */
  void
  interpToNewGrids(EBAMRCellData&                   a_newData,
                   const EBAMRCellData&             a_oldData,
                   const phase::which_phase         a_phase,
                   const int                        a_lmin,
                   const int                        a_oldFinestLevel,
                   const int                        a_newFinestLevel,
                   const EBCoarseToFineInterp::Type a_type);

  /**
    @brief Interpolate data to new grids
    @details This is called when requiring data to be interpolated to new grids. Takes old data as argument
    and fills the new grid data with an interpolation of the old grid data. 
    @param[out] a_newData        New grid data.
    @param[in]  a_oldData        Old grid data.
    @param[in]  a_phase          Phase on which we regrid. 
    @param[in]  a_lmin           Coarsest level that did not change (but distribution may have changed). 
    @param[in]  a_oldFinestLevel Previous finest level.
    @param[in]  a_newFinestLevel New finest level. 
    @param[in]  a_type           Interpolation type
  */
  void
  interpToNewGrids(EBAMRIVData&                     a_newData,
                   const EBAMRIVData&               a_oldData,
                   const phase::which_phase         a_phase,
                   const int                        a_lmin,
                   const int                        a_oldFinestLevel,
                   const int                        a_newFinestLevel,
                   const EBCoarseToFineInterp::Type a_type);

  /**
    @brief Regrid particle to new grids. 
    @note User needs to call ParticleContainer<P>::preRegrid prior to this.
    @param[in]  a_lmin           Coarsest level that did not change (but distribution may have changed). 
    @param[in]  a_newFinestLevel New finest level. 
    @param[in] a_particles Description
  */
  template <class P>
  void
  remapToNewGrids(ParticleContainer<P>& a_particles, const int a_lmin, const int a_newFinestLevel) const noexcept;

  /**
    @brief Regrid a struct-of-arrays particle container to new grids.
    @note User needs to call ParticleContainerSoA::preRegrid prior to this.
    @param[in,out] a_particles      SoA particle container to regrid.
    @param[in]     a_lmin           Coarsest level that did not change.
    @param[in]     a_newFinestLevel New finest level.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  remapToNewGrids(ParticleContainerSoA<P, Traits>& a_particles,
                  const int                        a_lmin,
                  const int                        a_newFinestLevel) const noexcept;

  /**
    @brief Interpolate data to centroids on realm and phase
    @param[in,out] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
    @details This leaves regular data intact, while the values in irregular are replaced by the corresponding interpolation to centroids. 
  */
  void
  interpToCentroids(EBAMRCellData& a_data, const std::string& a_realm, const phase::which_phase a_phase) const noexcept;

  /**
    @brief Interpolate data to centroids
    @param[in,out] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
    @param[in]    a_level Grid level
  */
  void
  interpToCentroids(LevelData<EBCellFAB>&    a_data,
                    const std::string&       a_realm,
                    const phase::which_phase a_phase,
                    const int                a_level) const noexcept;

  /**
    @brief Interpolate data to centroids
    @param[in,out] a_centroidData Centroid-centered data
    @param[in]    a_cellData Cell-centered data
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
    @param[in]    a_level Grid level
    @param[in]    a_din Grid index
  */
  void
  interpToCentroids(EBCellFAB&               a_centroidData,
                    const EBCellFAB&         a_cellData,
                    const std::string&       a_realm,
                    const phase::which_phase a_phase,
                    const int                a_level,
                    const DataIndex&         a_din) const noexcept;

  /**
    @brief Interpolate data to centroids on realm and phase
    @param[out] a_centroidData Centroid-centered data
    @param[in]  a_cellData Cell-centered data to be interpolated
    @param[in]  a_realm Realm name
    @param[in]  a_phase Phase (gas or solid)
    @details This leaves regular data intact, while the values in irregular are replaced by the corresponding interpolation to centroids. 
  */
  void
  interpToCentroids(EBAMRIVData&             a_centroidData,
                    const EBAMRCellData&     a_cellData,
                    const std::string&       a_realm,
                    const phase::which_phase a_phase) const noexcept;

  /**
    @brief Interpolate data to centroids
    @param[out] a_centroidData Centroid-centered data
    @param[in]  a_cellData Cell-centered data to be interpolated
    @param[in]  a_realm Realm name
    @param[in]  a_phase Phase (gas or solid)
    @param[in]  a_level Grid level
  */
  void
  interpToCentroids(LevelData<BaseIVFAB<Real>>& a_centroidData,
                    const LevelData<EBCellFAB>& a_cellData,
                    const std::string&          a_realm,
                    const phase::which_phase    a_phase,
                    const int                   a_level) const noexcept;

  /**
    @brief Interpolate data to EB centroids on realm and phase
    @param[out] a_centroidData EB-centered data
    @param[in]  a_cellData Cell-centered data to be interpolated
    @param[in]  a_realm Realm name
    @param[in]  a_phase Phase (gas or solid)
  */
  void
  interpToEB(EBAMRIVData&             a_centroidData,
             const EBAMRCellData&     a_cellData,
             const std::string&       a_realm,
             const phase::which_phase a_phase) const noexcept;

  /**
    @brief Interpolate data to EB centroids
    @param[out] a_centroidData EB-centered data
    @param[in]  a_cellData Cell-centered data to be interpolated
    @param[in]  a_realm Realm name
    @param[in]  a_phase Phase (gas or solid)
    @param[in]  a_level Grid level
  */
  void
  interpToEB(LevelData<BaseIVFAB<Real>>& a_centroidData,
             const LevelData<EBCellFAB>& a_cellData,
             const std::string&          a_realm,
             const phase::which_phase    a_phase,
             const int                   a_level) const noexcept;

  /**
    @brief Interpolate data to EB centroids
    @param[out] a_centroidData EB-centered data
    @param[in]  a_cellData Cell-centered data to be interpolated
    @param[in]  a_realm Realm name
    @param[in]  a_phase Phase (gas or solid)
    @param[in]  a_level Grid level
    @param[in]  a_din Grid index
  */
  void
  interpToEB(BaseIVFAB<Real>&         a_centroidData,
             const EBCellFAB&         a_cellData,
             const std::string&       a_realm,
             const phase::which_phase a_phase,
             const int                a_level,
             const DataIndex&         a_din) const noexcept;

  /**
    @brief Compute a non-conservative divergence.
    @details This computes divNC(F) = sum(kappa*divC(F))/sum(kappa), where kappa*divC(F) is the finite-volume approximation to
    kappa*div(F) = int_A F*dA (and dA is in grid units, not physical units). 
    @param[out] a_nonConsDivF The non-conservative divergence. 
    @param[in] a_kappaDivF Conservative divergence, not weighted by kappa.
    @param[in] a_realm Realm where the data lives.    
    @param[in] a_phase Phase where the data lives.
  */
  void
  nonConservativeDivergence(EBAMRIVData&              a_nonConsDivF,
                            const EBAMRCellData&      a_kappaDivF,
                            const std::string&        a_realm,
                            const phase::which_phase& a_phase) const noexcept;

  /**
    @brief Compute a non-conservative divergence.
    @details This computes divNC(F) = sum(kappa*divC(F))/sum(kappa), where kappa*divC(F) is the finite-volume approximation to
    kappa*div(F) = int_A F*dA (and dA is in grid units, not physical units). 
    @param[out] a_nonConsDivF The non-conservative divergence.
    @param[in] a_level Grid level
    @param[in] a_kappaDivF Conservative divergence, not weighted by kappa.
    @param[in] a_realm Realm where the data lives.    
    @param[in] a_phase Phase where the data lives.
  */
  void
  nonConservativeDivergence(LevelData<BaseIVFAB<Real>>& a_nonConsDivF,
                            const LevelData<EBCellFAB>& a_kappaDivF,
                            const int&                  a_level,
                            const std::string&          a_realm,
                            const phase::which_phase&   a_phase) const noexcept;

  /**
    @brief Sets multifluid index space. 
    @param[in] a_multiFluidIndexSpace Multifluid index space wrapper. 
  */
  void
  setMultifluidIndexSpace(const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace);

  /**
    @brief Set implicit function for a specific phase. Need e.g. for level-sets. 
    @param[in] a_phase  Phase
    @param[in] a_baseIF The level-set function describing the phase. 
  */
  void
  setBaseImplicitFunction(const phase::which_phase a_phase, const RefCountedPtr<BaseIF>& a_baseIF);

  /**
    @brief Run pre-regrid operations. 
  */
  void
  preRegrid();

  /**
    @brief Regrid AMR. This versions generates the grids and Realms, but not the operator. 
    @param[in] a_tags    Cell tags which will generate the grid. 
    @param[in] a_lmin    Coarsest grid level allowed to change. 
    @param[in] a_hardcap Grid generation hardcap. If < 0 there are no limitations to grid depth. 
  */
  void
  regridAmr(const Vector<IntVectSet>& a_tags, const int a_lmin, const int a_hardcap = -1);

  /**
    @brief Regrid a realm. This generates the grids for the realm, but does not do the operators on the realm. 
    @param[in] a_realm Realm name
    @param[in] a_procs Processor IDs
    @param[in] a_boxes Grid boxes
    @param[in] a_lmin  Coarsest level that changes. 
  */
  void
  regridRealm(const std::string&         a_realm,
              const Vector<Vector<int>>& a_procs,
              const Vector<Vector<Box>>& a_boxes,
              const int                  a_lmin);

  /**
    @brief Run post-regrid operations. 
  */
  void
  postRegrid();

  /**
    @brief Set grids from boxes and computational loads.
    @param[in] a_boxes          Grid boxes
    @param[in] a_realmsAndLoads Realm names and computational loads for realms
    @details This will set the grids for AmrMesh and realms, load balancing the various realms with the provided loads. 
  */
  void
  setGrids(const Vector<Vector<Box>>& a_boxes, const std::map<std::string, Vector<Vector<long int>>>& a_realmsAndLoads);

  /**
    @brief Regrid AMR operators. This is done for all realms. 
    @param[in] a_lmin Coarsest grid level that changes. 
  */
  void
  regridOperators(const int a_lmin);

  /**
    @brief Regrid AMR operators. This is done for a specific realm.
    @param[in] a_realm Realm
    @param[in] a_lmin  Coarsest grid level that changes. 
  */
  void
  regridOperators(const std::string& a_realm, const int a_lmin);

  /**
    @brief Do a sanity check to make sure everything is set up correctly
  */
  void
  sanityCheck() const;

  /**
    @brief Get finest grid level.
    @return Finest AMR level
  */
  int
  getFinestLevel() const;

  /**
    @brief Get maximum permitted amr depth.
    @return Maximum AMR depth
  */
  int
  getMaxAmrDepth() const;

  /**
    @brief Get maximum permitted simulation depth.
    @return Maximum simulation depth
  */
  int
  getMaxSimulationDepth() const;

  /**
    @brief Get blocking factor (i.e. the smallest possible box).
    @return Blocking factor
  */
  int
  getBlockingFactor() const;

  /**
    @brief Get maximum permitted box size
    @return Maximum box size
  */
  int
  getMaxBoxSize() const;

  /**
    @brief Return buffer for B-R mesh refinement algorithm.
    @return Buffer size
  */
  int
  getBrBuffer() const;

  /**
    @brief Get maximum box size for EBIS generation.
    @return Maximum EBIS box size
  */
  int
  getMaxEbisBoxSize() const;

  /**
    @brief Get the default number of ghost cells
    @return Number of ghost cells
  */
  int
  getNumberOfGhostCells() const;

  /**
    @brief Get number of ghost cells used for EB grid generation.
    @return Number of EB ghost cells
  */
  int
  getNumberOfEbGhostCells() const;

  /**
    @brief Get default redistribution radius
    @return Redistribution radius
  */
  int
  getRedistributionRadius() const;

  /**
    @brief Get resolution on the finest grid level.
    @return Grid spacing on finest level
  */
  Real
  getFinestDx() const;

  /**
    @brief Get lower-left corner of computational domain.
    @return Lower-left corner of the domain
  */
  RealVect
  getProbLo() const;

  /**
    @brief Get upper-right corner of computational domain.
    @return Upper-right corner of the domain
  */
  RealVect
  getProbHi() const;

  /**
    @brief Get finest domain
    @return Problem domain on the finest level
  */
  ProblemDomain
  getFinestDomain() const;

  /**
    @brief Get spatial resolutions
    @return Grid spacings on all levels
  */
  const Vector<Real>&
  getDx() const;

  /**
    @brief Get refinement ratios
    @return Refinement ratios between consecutive levels
  */
  const Vector<int>&
  getRefinementRatios() const;

  /**
    @brief Get refinement factor between two levels
    @param[in] a_level1 First level
    @param[in] a_level2 Second level
    @details This will work correctly in reverse as well.
    @return Refinement ratio
  */
  int
  getRefinementRatio(const int a_level1, const int a_level2) const;

  /**
    @brief Get implicit function for a specific phase
    @param[in] a_phase Phase (gas or solid)
    @return Implicit function for the specified phase
  */
  const RefCountedPtr<BaseIF>&
  getBaseImplicitFunction(const phase::which_phase a_phase) const;

  /**
    @brief Get the "proxy" grids in AmrMesh.
    @details Returns m_grids
    @return Proxy grids
  */
  const Vector<DisjointBoxLayout>&
  getProxyGrids() const;

  /**
    @brief Get the grids
    @param[in] a_realm Realm name
    @return Grids for the specified realm
  */
  const Vector<DisjointBoxLayout>&
  getGrids(const std::string& a_realm) const;

  /**
    @brief Get domains
    @return Problem domains on all AMR levels
  */
  const Vector<ProblemDomain>&
  getDomains() const;

  /**
    @brief Get a registered mask.
    @param[in] a_mask   Mask name
    @param[in] a_buffer Buffer (used when registering mask)
    @param[in] a_realm  Realm name
    @return Mask
  */
  const AMRMask&
  getMask(const std::string& a_mask, const int a_buffer, const std::string& a_realm) const;

  /**
    @brief Get a map of all valid cells on a specified realm
    @param[in] a_realm Realm name
    @return Valid cells
  */
  const AMRMask&
  getValidCells(const std::string& a_realm) const;

  /**
    @brief Get the tiled space representation
    @param[in] a_realm Realm name
    @return Level tiles
  */
  const Vector<RefCountedPtr<LevelTiles>>&
  getLevelTiles(const std::string& a_realm) const;

  /**
    @brief Get the EBLevelGrid for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return EBLevelGrids on all AMR levels
  */
  const Vector<RefCountedPtr<EBLevelGrid>>&
  getEBLevelGrid(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get the coarsened fine EBLevelGrid for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Coarsened-fine EBLevelGrids on all AMR levels
  */
  const Vector<RefCountedPtr<EBLevelGrid>>&
  getEBLevelGridCoFi(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get EBISLayouts for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return EBISLayouts on all AMR levels
  */
  const Vector<EBISLayout>&
  getEBISLayout(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get EBISLayouts for a Realm
    @param[in] a_realm Realm name
    @return M f level grid
  */
  const Vector<RefCountedPtr<MFLevelGrid>>&
  getMFLevelGrid(const std::string& a_realm) const;

  /**
    @brief Get VoF iterators covering all cut-cells (singly- and multiply-cut) for a Realm and phase.
    @details Not const because we need to reset the VoFIterator. For a pass covering only multiply-cut
             cells, use getMultiCutVofIterator instead.
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Per-level LayoutData of VoFIterators covering all irregular cells
  */
  Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
  getVofIterator(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get vof iterators restricted to multiply-cut cells for a Realm and phase.
    @details Multiply-cut cells are those with more than one VoF per cell. This iterator is used
             as a second pass after BoxLoops::loop over getSingleValuedFAB, which already covers
             VoF 0 of all cut-cells.
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Per-level LayoutData of VoFIterators restricted to multiply-cut cells
  */
  Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
  getMultiCutVofIterator(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get face iterators over cut-cell faces in each valid grid box (FaceStop::SurroundingWithBoundary).
    @details Each per-box entry is an array of SpaceDim FaceIterators. Domain boundary faces are included.
             For cell-to-face averaging that excludes domain boundary faces, use getFaceIteratorNoBoundary.
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Per-level LayoutData of FaceIterator arrays
  */
  Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
  getFaceIterator(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get face iterators over cut-cell faces in each valid grid box (FaceStop::SurroundingNoBoundary).
    @details Each per-box entry is an array of SpaceDim FaceIterators. Domain boundary faces are excluded.
             Intended for cell-to-face averaging where boundary faces are fixed up separately.
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Per-level LayoutData of FaceIterator arrays
  */
  Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
  getFaceIteratorNoBoundary(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get per-direction face iterators covering cut-cell faces that reach into the tangential ghost layer.
    @details Each per-box entry is an array of SpaceDim FaceIterators, one per coordinate direction.
             For the iterator in direction dir, the underlying IVS covers irregular cells in a box grown by 1
             in the two directions tangential to dir; the dir extent stays at the valid box. Domain boundary
             faces are excluded (FaceStop::SurroundingNoBoundary). Callers require at least 1 ghost layer in
             cell data.
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Per-level LayoutData of FaceIterator arrays
  */
  Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
  getFaceIteratorWithTangentialGhosts(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get face iterators restricted to multiply-cut faces in each valid grid box.
    @details Only faces abutting at least one multiply-cut cell (cells with more than one VoF)
             are visited. Faces are enumerated with FaceStop::SurroundingNoBoundary. Used as a
             second pass after a box loop over getSingleValuedFAB to avoid double-processing
             singly-cut faces.
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Per-level LayoutData of FaceIterator arrays
  */
  Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
  getMultiCutFaceIterator(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get levelset function, allocated over a grid for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Levelset
  */
  const EBAMRFAB&
  getLevelset(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get the regular-cell mask (value 1 in regular cells, 0 elsewhere) for a Realm and phase.
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Regular-cell mask
  */
  const EBAMRCellData&
  getRegularCells(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get the covered-cell mask (value 1 in covered cells, 0 elsewhere) for a Realm and phase.
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Covered-cell mask
  */
  const EBAMRCellData&
  getCoveredCells(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get the non-covered-cell mask (1 in regular/irregular cells, 0 in covered cells) for a Realm and phase.
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Non-covered-cell mask
  */
  const EBAMRCellData&
  getNotCoveredCells(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get the irregular-cell mask (value 1 in irregular cells, 0 elsewhere) for a Realm and phase.
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Irregular-cell mask
  */
  const EBAMRCellData&
  getIrregularCells(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get EBAMRParticleMesh operator
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)    
    @return Particle mesh
  */
  EBAMRParticleMesh&
  getParticleMesh(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get the SoA EBAMRParticleMesh operator
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return SoA particle mesh
  */
  EBAMRParticleMeshSoA&
  getParticleMeshSoA(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get EBAMRSurfaceDeposition surface deposition operator
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)    
    @return Surface deposition
  */
  EBAMRSurfaceDeposition&
  getSurfaceDeposition(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get the coarsening utility
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Coarse-averaging operators on all levels
  */
  Vector<RefCountedPtr<EBCoarAve>>&
  getCoarseAverage(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get multigrid interpolation utility
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Multigrid interpolator
  */
  Vector<RefCountedPtr<EBMultigridInterpolator>>&
  getMultigridInterpolator(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get interpolator
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Coarse-to-fine interpolation operators on all levels
  */
  Vector<RefCountedPtr<EBCoarseToFineInterp>>&
  getFineInterp(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get flux register
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Flux registers on all AMR levels
  */
  Vector<RefCountedPtr<EBReflux>>&
  getFluxRegister(const std::string& a_realm, const phase::which_phase a_phase) const;

  /**
    @brief Get the redistribution operators
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
    @return Redistribution op
  */
  Vector<RefCountedPtr<EBFluxRedistribution>>&
  getRedistributionOp(const std::string& a_realm, const phase::which_phase a_phase) const;

#ifdef CH_USE_PETSC
  /**
    @brief Get a PETSc grid object
    @param[in] a_realm Realm name.
  */
  const RefCountedPtr<PetscGrid>&
  getPetscGrid(const std::string& a_realm) const noexcept;
#endif

  /**
    @brief Get the name of all Realms
    @return Names of all Realms
  */
  std::vector<std::string>
  getRealms() const;

  /**
    @brief Get box sorting method
    @return Box sorting method.
  */
  BoxSorting
  getBoxSorting() const;

  /**
    @brief Get EBIndexSpace corresponding to a particular phase
    @param[in] a_phase Phase (gas/solid)
    @return EBIndexSpace for the specified phase
  */
  RefCountedPtr<EBIndexSpace>&
  getEBIndexSpace(const phase::which_phase a_phase)
  {
    return m_multifluidIndexSpace->getEBIndexSpace(a_phase);
  }

protected:
  /**
    @brief Enum for having understandable notation for grid generation
  */
  enum class GridGenerationMethod
  {
    BergerRigoutsous,
    Tiled,
  };

  /**
    @brief These are all the Realms
  */
  mutable std::map<std::string, RefCountedPtr<Realm>> m_realms;

  /**
    @brief Implicit functions
  */
  std::map<phase::which_phase, RefCountedPtr<BaseIF>> m_baseif;

  /**
    @brief Old grids
  */
  std::map<std::string, Vector<DisjointBoxLayout>> m_oldGrids;

  /**
    @brief Storage for copiers from the old grids to the new ones
    @note Defined in the realm regrid routine. Cleared in the postRegrid routine. 
  */
  std::map<std::string, Vector<Copier>> m_oldToNewCellCopiers;

  /**
    @brief Storage for copiers from the old grids to the new ones
    @note Defined in the realm regrid routine. Cleared in the postRegrid routine. 
  */
  std::map<std::string, Vector<Copier>> m_oldToNewEBCopiers;

  /**
    @brief Map for copying between various Realms. First index is the "from" realm and second index is the "to" realm
    @note This is important at scale when you have > 1M boxes since Copier definitions become costly.
  */
  std::map<std::pair<std::string, std::string>, Vector<Copier>> m_validToValidRealmCopiers;

  /**
    @brief Map for copying between various Realms. First index is the "from" realm and second index is the "to" realm
    @note This is important at scale when you have > 1M boxes since Copier definitions become costly.
  */
  std::map<std::pair<std::string, std::string>, Vector<Copier>> m_validToValidGhostRealmCopiers;

  /**
    @brief Map for copying between various Realms. First index is the "from" realm and second index is the "to" realm
    @note This is important at scale when you have > 1M boxes since Copier definitions become costly.
  */
  std::map<std::pair<std::string, std::string>, Vector<Copier>> m_validGhostToValidRealmCopiers;

  /**
    @brief Map for copying between various Realms. First index is the "from" realm and second index is the "to" realm
    @note This is important at scale when you have > 1M boxes since Copier definitions become costly.
  */
  std::map<std::pair<std::string, std::string>, Vector<Copier>> m_validGhostToValidGhostRealmCopiers;

  /**
    @brief Grid generation method
  */
  GridGenerationMethod m_gridGenerationMethod;

  /**
    @brief Box sorting
  */
  BoxSorting m_boxSort;

  /**
    @brief MultiFluidIndexSpace
  */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /**
    @brief Interpolation type for cell to cell centroid interpolation
  */
  CellCentroidInterpolation::Type m_cellCentroidInterpolationType;

  /**
    @brief Interpolation type for cell to EB centroid interpolation    
  */
  EBCentroidInterpolation::Type m_ebCentroidInterpolationType;

  /**
    @brief Coarsest box where we compute
  */
  IntVect m_numCells;

  /**
    @brief Fill ratio
  */
  Real m_fillRatioBR;

  /**
    @brief Domain simulation corner
  */
  RealVect m_probLo;

  /**
    @brief Domain simulation corner
  */
  RealVect m_probHi;

  /**
    @brief Refinement ratio
  */
  int m_refRatio;

  /**
    @brief Verbosity
  */
  int m_verbosity;

  /**
    @brief Finest level
  */
  int m_finestLevel;

  /**
    @brief Finest level before a regrid
  */
  int m_oldFinestLevel;

  /**
    @brief Maximum amr depth
  */
  int m_maxAmrDepth;

  /**
    @brief Maximum allowed depth for simulation
  */
  int m_maxSimulationDepth;

  /**
    @brief Max box size
  */
  int m_maxBoxSize;

  /**
    @brief Maximum box size for EBIS generation
  */
  int m_maxEbisBoxSize;

  /**
    @brief Set buffer size
  */
  int m_bufferSizeBR;

  /**
    @brief Blocking factor
  */
  int m_blockingFactor;

  /**
    @brief Number of ghost cells to use for eb stuff
  */
  int m_numEbGhostsCells;

  /**
    @brief Number of ghost cells
  */
  int m_numGhostCells;

  /**
    @brief Number of ghost cells to use when writing level-set to grid
  */
  int m_numLsfGhostCells;

  /**
    @brief Multigrid interpolation order
  */
  int m_multigridInterpOrder;

  /**
    @brief Multigrid interpolation radius
  */
  int m_multigridInterpRadius;

  /**
    @brief Multigrid interpolation weights
  */
  int m_multigridInterpWeight;

  /**
    @brief Redistribution radius
  */
  int m_redistributionRadius;

  /**
    @brief Has grids or not
  */
  bool m_hasGrids;

  /**
    @brief Has regrid copiers or not
  */
  bool m_hasRegridCopiers;

  /**
    @brief Grids
  */
  Vector<DisjointBoxLayout> m_grids;

  /**
    @brief Problem domains
  */
  Vector<ProblemDomain> m_domains;

  /**
    @brief AMR resolutions
  */
  Vector<int> m_refinementRatios;

  /**
    @brief Level resolutions
  */
  Vector<Real> m_dx;

  /**
    @brief Define Realms
  */
  void
  defineRealms();

  /**
    @brief Build new internal AMR grids.
    @param[in,out] a_tags Sets of cell tags used for the grid generation.
    @param[in]    a_lmin The finest grid level which changes. 
    @param[in]    a_hardcap Hardcap for the maximum grid level which can be generated. If a_hardcap < 0 there is no restriction beyond the AmrMesh restrictions.
    @details This will call the specified grid generation method and fill m_grids with the new grids. 
  */
  void
  buildGrids(const Vector<IntVectSet>& a_tags, const int a_lmin, const int a_hardcap = -1);

  /**
    @brief Build copiers for copying between realms
  */
  void
  buildCopiers();

  /**
    @brief Parse the low/high corners of the computational domain. 
  */
  void
  parseProbLoHiCorners();

  /**
    @brief Parse the grid generation algorithm. 
  */
  void
  parseGridGeneration();

  /**
    @brief Parse the verbosity for AmrMesh. 
  */
  void
  parseVerbosity();

  /**
    @brief Parse the coarsest domain grid.
  */
  void
  parseCoarsestLevelNumCells();

  /**
    @brief Parse the maximum permitted AMR depth
  */
  void
  parseMaxAmrDepth();

  /**
    @brief Set maximum simulation depth. 
  */
  void
  parseMaxSimulationDepth();

  /**
    @brief Parse refinement ratios
  */
  void
  parseRefinementRatios();

  /**
    @brief Parse buffer size for Berger-Rigoutsous grid algorithm. 
  */
  void
  parseBrBufferSize();

  /**
    @brief Parse the Berger-Rigoutsos fill ratio. 
  */
  void
  parseBrFillRatio();

  /**
    @brief Parse the maximum permitted box size. 
  */
  void
  parseMaxBoxSize();

  /**
    @brief Parse the maximum permitted box size. 
    @details Only relevant for the grid generation step. 
  */
  void
  parseMaxEbisBoxSize();

  /**
    @brief Parse the B-R blocking factor. For tiled mesh refinement this parses the tile size. 
  */
  void
  parseBlockingFactor();

  /**
    @brief Parse number of ghost cells for eb stuff
    @details Only matters for the grid generation step. 
  */
  void
  parseEbGhostCells();

  /**
    @brief Parse the number of ghost cells 
  */
  void
  parseNumGhostCells();

  /**
    @brief Parse settings for the multigrid interpolator
  */
  void
  parseMultigridInterpolator();

  /**
    @brief Parse the default redistribution radius
  */
  void
  parseRedistributionRadius();

  /**
    @brief Parse centroid interpolation stencils
  */
  void
  parseCellCentroidInterpolation();

  /**
    @brief Parse EB interpolation (or extrapolation) stencils
  */
  void
  parseEBCentroidInterpolation();

  /**
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
                  const std::string&          a_realm,
                  const phase::which_phase    a_phase,
                  const int                   a_lvl) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_AmrMeshImplem.H>

#endif
