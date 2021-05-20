/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   AmrMesh.H
  @brief  Declaration of core class for handling AMR-related operations (with embedded boundaries)
  @author Robert Marskar
  @todo   a_lmax is a dead parameter in buildGrids. Fix this and the dependencies. 
  @todo   a_regsize is a dead parameter which should be removed. 
*/

#ifndef CD_AmrMesh_H
#define CD_AmrMesh_H

// Chombo includes
#include <DisjointBoxLayout.H>
#include <ebcoarseaverage.H>
#include <EBFluxRegister.H>
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
#include <CD_ComputationalGeometry.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_IrregAmrStencil.H>
#include <CD_LoadBalancing.H>
#include <CD_NwoEbQuadCfInterp.H>
#include <MFLevelGrid.H>
#include <ebcoarseaverage.H>
#include <CD_DomainFluxIFFAB.H>
#include <CD_DomainFluxIFFABFactory.H>
#include <CD_ParticleContainer.H>
#include <CD_EbGhostCloud.H>
#include <CD_NonConservativeDivergenceStencil.H>
#include <CD_Realm.H>
#include <CD_LoadBalancing.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Class for handling spatial operations for streamer equations.
  @details This class handles most (if not all) spatial operations. Examples are grid generation, coarsening of data, interpolation of data, memory allocations etc. 
*/
class AmrMesh {
protected:
  
  /*!
    @brief Enum for specifying how to interpolate ghost cells. 
  */
  enum class GhostInterpolation {
    PiecewiseLinear,
    Quadratic
  };
  
public:

  /*!
    @brief Constructor
  */
  AmrMesh();

  /*!
    @brief Destructor
  */
  ~AmrMesh();

  /*!
    @brief Parse options. Called during the constructor. 
  */
  void parseOptions();

  /*!
    @brief Parse runtime options. 
  */
  void parseRuntimeOptions();

  /*!
    @brief Deallocate data. 
    @param[inout] a_data Data to be deleted. 
    @details You should never have to do this, as smart pointer are our friends. 
  */
  template<typename T>
  void deallocate(Vector<T*>& a_data);

  /*!
    @brief Deallocate smart pointer data.
    @param[inout] a_data Data to be deallocated. 
  */
  template<typename T>
  void deallocate(Vector<RefCountedPtr<T> >& a_data);

  /*!
    @brief Deallocate data.
    @param[inout] a_data Data to be deallocated.
  */
  template<typename T>
  void deallocate(EBAMRData<T>& a_data);

  /*!
    @brief Turn smart-pointer data structure into regular-pointer data structure
    @param[out] a_alias Aliased data
    @param[in] a_data  Smart-pointer data
    @details This fetches the bare pointer from a_data and sets it in a_alias. 
  */
  template<typename T>
  void alias(Vector<T*>& a_alias, const Vector<RefCountedPtr<T> >& a_data);

  /*!
    @brief Turn smart-pointer data structure into regular-pointer data structure
    @param[out] a_alias Aliased data
    @param[in]  a_data  Smart-pointer data
    @details This fetches the bare pointer from a_data and sets it in a_alias. 
  */
  template<typename T, typename S>
  void alias(Vector<T*>& a_alias, const EBAMRData<S>& a_data);

  /*!
    @brief Template class for generic allocation of particles
    @param[out] a_particles Particle data holder. 
    @param[out] a_realm     Realm on which the particles will be allocated. 
  */
  template <typename T >
  void allocate(Vector<RefCountedPtr<ParticleData<T> > >& a_particles, const std::string a_realm);

  /*!
    @brief Allocate a particle container
    @param[out] a_container Particle container to be allocated
    @param[in]  a_pvrBuffer The particle valid region buffer. 
    @param[out] a_realm     Realm on which the particles will be allocated. 
    @details Particles that live within a_pvrBuffer on the fine side of the refinement boundary will be mapped to the coarse leve. 
  */
  template <typename T>
  void allocate(ParticleContainer<T>& a_container,
		const int              a_pvrBuffer,
		const std::string      a_realm);

  /*!
    @brief Allocate pointer but not any memory.
    @param[out] a_data Data structure for which pointers will be constructed. 
    @details This sets all the elements in the vector to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void allocatePointer(Vector<RefCountedPtr<T> >& a_data);

  /*!
    @brief Allocate pointer but not any memory.
    @param[out] a_data       Data structure for which pointers will be constructed. 
    @param[in] a_finestLevel Finest grid level on which a_data will be defined. 
    @details This sets all the elements in the vector to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void allocatePointer(Vector<RefCountedPtr<T> >& a_data, const int a_finestLevel);

  /*!
    @brief Allocate pointer but not any memory.
    @param[out] a_data Data structure for which pointers will be constructed. 
    @details This sets all the elements in the vector to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void allocatePointer(EBAMRData<T>& a_data);

  /*!
    @brief Allocate pointer but not any memory.
    @param[out] a_data       Data structure for which pointers will be constructed. 
    @param[in] a_finestLevel Finest grid level on which a_data will be defined. 
    @details This sets all the elements in the vector to be RefCountedPtr<T>(nullptr)
  */
  template <typename T>
  void allocatePointer(EBAMRData<T>& a_data, const int a_finestLevel);

  /*!
    @brief Register an operator over a Realm and a phase. 
    @param[in] a_operator Operator to be registered. 
    @param[in] a_realm    Realm name
    @param[in] a_phase    Phase (gas or solid). 
  */
  void registerOperator(const std::string a_operator, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Register a boolean mask over a Realm. 
    @param[in] a_mask     Mask name to be registered. 
    @param[in] a_phase    Phase (gas or solid). 
    @param[in] a_buffer   Relevant buffer for the mask.
    @note The meaning of a_buffer differs for different masks. 
  */
  void registerMask(const std::string a_mask, const int a_buffer, const std::string a_realm);

  /*!
    @brief Query if a Realm exists
    @param[in] a_realm Name of the Realm. 
    @return True if the Realm exists and false otherwise. 
  */
  bool queryRealm(const std::string a_realm) const;

  /*!
    @brief Register a new Realm
    @param[in] a_realm Realm name
  */
  void registerRealm(const std::string a_realm);

  /*!
    @brief Alias multifluid data up to a specified level. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @param[in]  a_finestLevel Maximum grid level for aliasing operation. 
  */
  void alias(EBAMRCellData&           a_data,
	     const phase::which_phase a_phase,
	     const MFAMRCellData&     a_mfdata,
	     const int                a_finestLevel);

  /*!
    @brief Alias multifluid data up to a specified level. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @param[in]  a_finestLevel Maximum grid level for aliasing operation. 
  */
  void alias(EBAMRFluxData&           a_data,
	     const phase::which_phase a_phase,
	     const MFAMRFluxData&     a_mfdata,
	     const int                a_finestLevel);

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
  */
  void alias(EBAMRCellData& a_data, const phase::which_phase a_phase, const MFAMRCellData& a_mfdata);

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
  */
  void alias(EBAMRFluxData& a_data, const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata);

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[out] a_data        Single-phased data on specific phase
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
  */
  void alias(EBAMRIVData& a_data, const phase::which_phase a_phase, const MFAMRIVData& a_mfdata);

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @return A single-phase alias of the multifluid data, having the same size as the input data. 
  */
  EBAMRCellData alias(const phase::which_phase a_phase, const MFAMRCellData& a_mfdata);

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @return A single-phase alias of the multifluid data, having the same size as the input data. 
  */
  EBAMRFluxData alias(const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata);

  /*!
    @brief Alias an entire hiearchy of multifluid data. 
    @param[in]  a_phase       Phase (gas or solid).
    @param[in]  a_mfdata      Multifluid data to be aliased
    @return A single-phase alias of the multifluid data, having the same size as the input data. 
  */
  EBAMRIVData alias(const phase::which_phase a_phase, const MFAMRIVData& a_mfdata);

  /*!
    @brief Compute cell-centered gradient for a grid level.
    @param[out] a_gradient Cell centered gradient. 
    @param[in]  a_phi      The scalar for which the gradient is computed. 
    @param[in]  a_realm    Name of the Realm where the data lives. 
    @param[in]  a_phase    Phase on which the data lives. 
    @param[in]  a_lvl      Grid level where the data lives. 
    @note This routine will reach into ghost cells and across refinement boundaries. The user must make sure that 
    ghost cells updated before using this routine. 
  */
  void computeGradient(LevelData<EBCellFAB>&       a_gradient,
		       const LevelData<EBCellFAB>& a_phi,
		       const std::string           a_realm,
		       const phase::which_phase    a_phase,
		       const int                   a_lvl);

  /*!
    @brief Compute cell-centered gradient over an AMR hierarchy. 
    @param[out] a_gradient Cell centered gradient. 
    @param[in]  a_phi      The scalar for which the gradient is computed. 
    @param[in]  a_realm    Name of the Realm where the data lives. 
    @param[in]  a_phase    Phase on which the data lives. 
    @note This routine will reach into ghost cells and across refinement boundaries. The user must make sure that 
    ghost cells updated before using this routine. 
  */
  void computeGradient(EBAMRCellData&           a_gradient,
		       const EBAMRCellData&     a_phi,
		       const std::string        a_realm,  
		       const phase::which_phase a_phase);

  /*!
    @brief Compute cell-centered gradient over an AMR hierarchy. 
    @param[out] a_gradient Cell centered gradient. 
    @param[in]  a_phi      The scalar for which the gradient is computed. 
    @param[in]  a_realm    Name of the Realm where the data lives. 
    @note This routine will reach into ghost cells and across refinement boundaries. The user must make sure that 
    ghost cells updated before using this routine. 
  */
  void computeGradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi, const std::string a_realm);

  /*!
    @brief Allocate a data holder over a specific Realm
    @param[out] a_data
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_phase  Phase (gas or solid)
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void allocate(EBAMRCellData&           a_data,
		const std::string        a_realm,
		const phase::which_phase a_phase,
		const int                a_nComp,
		const int                a_nGhost = -1);

  /*!
    @brief Allocate a data holder over a specific Realm
    @param[out] a_data
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_phase  Phase (gas or solid)
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void allocate(EBAMRFluxData&           a_data,
		const std::string        a_realm,
		const phase::which_phase a_phase,
		const int                a_nComp,
		const int                a_nGhost = -1);

  /*!
    @brief Allocate a data holder over a specific Realm
    @param[out] a_data
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_phase  Phase (gas or solid)
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void allocate(EBAMRIVData&             a_data,
		const std::string        a_realm,
		const phase::which_phase a_phase,
		const int                a_nComp,
		const int                a_nGhost = -1);

  /*!
    @brief Allocate a data holder over a specific Realm
    @param[out] a_data
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void allocate(EBAMRIFData&             a_data,
		const std::string        a_realm,
		const phase::which_phase a_phase,
		const int                a_nComp,
		const int                a_nGhost = -1);

  /*!
    @brief Allocate Boolean data over a specific Realm
    @param[out] a_data
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void allocate(EBAMRBool& a_data, const std::string a_realm, const int a_nComp, const int a_nGhost = -1);
  
  /*!
    @brief Allocate a data holder over a specific Realm.
    @param[out] a_data
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void allocate(MFAMRCellData& a_data, const std::string a_realm, const int a_nComp, const int a_nGhost = -1);

  /*!
    @brief Allocate a data holder over a specific Realm.
    @param[out] a_data
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void allocate(MFAMRFluxData& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Allocate a data holder over a specific Realm
    @param[out] a_data
    @param[in]  a_realm  Realm of the name where the data will be allocated. 
    @param[in]  a_nComp  Number of components in a_data
    @param[in]  a_nGhost Number of ghost cells for a_data
    @note If a_nGhost < 0, this routine will use the default number of ghost cells. 
  */
  void allocate(MFAMRIVData& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void reallocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void reallocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void reallocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_phase Phase over the reallocation.
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void reallocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void reallocate(EBAMRBool& a_data, const int a_lmin);

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void reallocate(MFAMRCellData& a_data, const int a_lmin);

  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void reallocate(MFAMRFluxData& a_data, const int a_lmin);
  
  /*!
    @brief Reallocate data. 
    @param[out] a_data  Data to be reallocated
    @param[in]  a_lmin  The lowest grid level that will be reallocated. 
    @details This will reallocate data for all grid levels lvl >= a_lmin. 
  */
  void reallocate(MFAMRIVData& a_data, const int a_lmin);

  /*!
    @brief Average down multifluid data over a Realm. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
  */
  void averageDown(MFAMRCellData& a_data, const std::string a_realm);

  /*!
    @brief Average down multifluid data over a Realm. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
  */
  void averageDown(MFAMRFluxData& a_data, const std::string a_realm);

  /*!
    @brief Average down on specific Realm and phase. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void averageDown(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Average down on specific Realm and phase. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
    @param[in]    a_level The coarse level.
    @details This coarsens data from level a_lvl+1 to level a_lvl
  */
  void averageDown(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase, const int a_lvl);

  /*!
    @brief Average down on specific Realm and phase. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void averageDown(EBAMRFluxData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Average down on specific Realm and phase. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void averageDown(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Do conservative coarsening for data defined on cut-cells. 
    @param[inout] a_data  Data to be coarsened.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void conservativeAverage(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Build new internal AMR grids.
    @param[inout] a_tags Sets of cell tags used for the grid generation.
    @param[in]    a_lmin The finest grid level which does not change. 
    @param[in]    a_lmax This is a dead parameter
    @param[in]    a_hardcap Hardcap for the maximum grid level which can be generated. If a_hardcap < 0 there is no restriction beyond the AmrMesh restrictions.
    @details This will call the specified grid generation method and fill m_grids with the new grids. 
  */
  void buildGrids(Vector<IntVectSet>& a_tags, const int a_lmin, const int a_lmax, const int a_hardcap = -1);

  /*!
    @brief Interpolate ghost vectors over a Realm, using the default ghost cell interpolation metho. 
    @param[inout] a_data  Data to be interpolated.
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void interpGhost(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Interpolate ghost cells over a Realm, using the default ghost cell interpolation method on a specific level. 
    @param[inout] a_fineData Fine grid data
    @param[inout] a_coarData Coarse grid data
    @param[in]    a_level    The grid level corresponding to a_fineData
    @param[in]    a_realm    Realm name
    @param[in]    a_phase    Phase (gas or solid)
  */
  void interpGhost(LevelData<EBCellFAB>&       a_fineData,
		   const LevelData<EBCellFAB>& a_coarData,
		   const int                   a_level,
		   const std::string           a_realm,
		   const phase::which_phase    a_phase);

  /*!
    @brief Interpolate ghost cells over a Realm. Calls the default ghost cell interpolation method. 
    @param[inout] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
  */
  void interpGhost(MFAMRCellData& a_data, const std::string a_realm);

  /*!
    @brief Interpolate ghost cells over a Realm and phase. This uses quadratic interpolation.
    @param[inout] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void interpGhostQuad(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Interpolate ghost cells over a Realm and phase. This uses piecewise linear interpolation.
    @param[inout] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
  */
  void interpGhostPwl(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Interpolate data to centroids on Realm and phase
    @param[inout] a_data  Data to be interpolated
    @param[in]    a_realm Realm name
    @param[in]    a_phase Phase (gas or solid)
    @details This leaves regular data intact, while the values in irregular are replaced by the corresponding interpolation to centroids. 
  */
  void interpToCentroids(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);
  
  /*!
    @brief Sets multifluid index space. 
    @param[in] a_multiFluidIndexSpace Multifluid index space wrapper. 
  */
  void setMultifluidIndexSpace(const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace);

  /*!
    @brief Set implicit function for a specific phase. Need e.g. for level-sets. 
    @param[in] a_phase  Phase
    @param[in] a_baseif The level-set function describing the phase. 
  */
  void setBaseImplicitFunction(const phase::which_phase a_phase, const RefCountedPtr<BaseIF>& a_baseif);

  /*!
    @brief Parse the low/high corners of the computational domain. 
  */
  void parseProbLoHiCorners();

  /*!
    @brief Set ghost cell interpolation type. Either piecewise linear or quadratic. 
  */
  void parseGhostInterpolation();

  /*!
    @brief Build domains
  */
  void buildDomains();

  /*!
    @brief Parse the grid generation algorithm. 
  */
  void parsegridGeneration();

  /*!
    @brief Parse the verbosity for AmrMesh. 
  */
  void parseVerbosity();

  /*!
    @brief Parse the coarsest domain grid.
  */
  void parseCoarsestLevelNumCells();

  /*!
    @brief Parse the maximum permitted AMR depth
  */
  void parseMaxAmrDepth();

  /*!
    @brief Set maximum simulation depth. 
  */
  void parseMaxSimulationDepth();

  /*!
    @brief Parse whether or not AmrMesh to forget about ebcf (user must guarantee that this is the case)
  */
  void parseEbCf();

  /*!
    @brief Parse refinement ratios
  */
  void parseRefinementRatios();

  /*!
    @brief Set the refinement ratios
  */
  void setRefinementRatios(const Vector<int> a_ref_ratios);

  /*!
    @brief Parse buffer size for Berger-Rigoutsous grid algorithm. 
  */
  void parseBrBufferSize();

  /*!
    @brief Parse how much to grow cut-cell geometric tags. 
  */
  void parseIrregTagGrowth();

  /*!
    @brief Parse the Berger-Rigoutsos fill ratio. 
  */
  void parseBrFillRatio();

  /*!
    @brief Set the finest level
  */
  void setFinestLevel(const int a_finestLevel);

  /*!
    @brief Parse the maximum permitted box size. 
  */
  void parseMaxBoxSize();

  /*!
    @brief Parse the maximum permitted box size. 
    @details Only relevant for the grid generation step. 
  */
  void parseMaxEbisBoxSize();

  /*!
    @brief Parse the B-R blocking factor. For tiled mesh refinement this parses the tile size. 
  */
  void parseBlockingFactor();

  /*!
    @brief Parse number of ghost cells for eb stuff
    @details Only matters for the grid generation step. 
  */ 
  void parseEbGhostCells();

  /*!
    @brief Parse the number of ghost cells 
  */ 
  void parseNumGhostCells();

  /*!
    @brief Parse the default redistribution radius
  */
  void parseRedistributionRadius();

  /*!
    @brief Parse centroid interpolation stencils
  */
  void parseCentroidStencils();

  /*!
    @brief Parse EB interpolation (or extrapolation) stencils
  */
  void parseEbCentroidStencils();

  /*!
    @brief Regrid function. Generate new grids and operators.
    @param[in] a_tags    Cell tags which will generate the grid. 
    @param[in] a_lmin    Coarsest grid level allowed to change. 
    @param[in] a_lmax    Dead parameter. 
    @param[in] a_regsize Redistribution register size. 
    @param[in] a_hardcap Grid generation hardcap. 
  */
  void regrid(const Vector<IntVectSet>& a_tags,
	      const int                 a_lmin,
	      const int                 a_lmax,
	      const int                 a_regsize=1, 
	      const int                 a_hardcap=-1);

  /*!
    @brief Regrid AMR. This versions generates the grids and Realms, but not the operator. 
    @param[in] a_tags    Cell tags which will generate the grid. 
    @param[in] a_lmin    Coarsest grid level allowed to change. 
    @param[in] a_lmax    Dead parameter. 
    @param[in] a_hardcap Grid generation hardcap. 
  */
  void regridAmr(const Vector<IntVectSet>& a_tags,
		 const int                 a_lmin,
		 const int                 a_lmax,
		 const int                 a_hardcap=-1);

  /*!
    @brief Define a Realm. This generates the grids for the ralm, but not the operators. 
    @param[in] a_realm Realm name
    @param[in] a_procs Processor IDs
    @param[in] a_boxes Grid boxes
    @param[in] a_lmin  Coarsest level allowed to change. 
  */
  void regridRealm(const std::string           a_realm,
		   const Vector<Vector<int> >& a_procs,
		   const Vector<Vector<Box> >& a_boxes,
		   const int                   a_lmin);

  /*!
    @brief Set grids from boxes and computational loads.
    @param[in] a_boxes Grid boxes
    @param[in] a_loads Computational loads
    @details This will run a load balancing algorithm and set m_grids accordingly. 
  */
  void setGrids(const Vector<Vector<Box> >& a_boxes, const std::map<std::string, Vector<Vector<long int> > >& a_loads);

  /*!
    @brief Regrid AMR operators. This is done for all Realms. 
    @param[in] a_lmin    Coarsest grid level that has changed. 
    @param[in] a_lmax    Dead parameter. 
    @param[in] a_regsize Dead parameter
  */
  void regridOperators(const int a_lmin,
		       const int a_lmax,
		       const int a_regsize=1);

  /*!
    @brief Do a sanity check to make sure everything is set up correctly
  */
  void sanityCheck() const;

  /*!
    @brief True if there is ebcf crossing
  */
  bool getEbCf() const;

  /*!
    @brief Get finest grid level. 
  */
  int getFinestLevel() const;

  /*!
    @brief Get irregular cell growth.
  */
  int getIrregTagGrowth() const;

  /*!
    @brief Get maximum permitted amr depth.
  */
  int getMaxAmrDepth() const;

  /*!
    @brief Get maximum permitted simulation depth.
  */
  int getMaxSimulationDepth() const;

  /*!
    @brief Get blocking factor (i.e. the smallest possible box).
  */
  int getBlockingFactor() const;

  /*!
    @brief Get maximum permitted box size
  */
  int getMaxBoxSize() const;

  /*!
    @brief Return buffer for B-R mesh refinement algorithm. 
  */
  int getBrBuffer() const;

  /*!
    @brief Get maximum box size for EBIS generation. 
  */
  int getMaxEbisBoxSize() const;

  /*!
    @brief Get the default number of ghost cells
  */
  int getNumberOfGhostCells() const;

  /*!
    @brief Get number of ghost cells used for EB grid generation. 
  */
  int getNumberOfEbGhostCells() const;

  /*!
    @brief Get default redistribution radius
  */
  int getRedistributionRadius() const;
  
  /*!
    @brief Get resolution on the finest grid level. 
  */
  Real getFinestDx() const;

  /*!
    @brief Get lower-left corner of computational domain. 
  */
  RealVect getProbLo() const;

  /*!
    @brief Get upper-left corner of computational domain. 
  */
  RealVect getProbHi() const;

  /*!
    @brief Get finest domain
  */
  ProblemDomain getFinestDomain() const;

  /*!
    @brief Get spatial resolutions
  */
  const Vector<Real>& getDx() const;

  /*!
    @brief Get refinement ratios
  */
  const Vector<int>& getRefinementRatios() const;

  /*!
    @brief Get refinement factor between two levels
    @param[in] a_level1 First level
    @param[in] a_level2 Second level
    @details This will work correctly in reverse as well.
  */
  int getRefinementRatio(const int a_level1, const int a_level2) const ;

  /*!
    @brief Get implicit function for a specific phase
    @param[in] a_phase Phase 8gas or solid)
  */
  const RefCountedPtr<BaseIF>& getBaseImplicitFunction(const phase::which_phase a_phase) const;

  /*!
    @brief Get all irregular tags.
    @details This collects the tags form EBIndexSpace, and the tags are local to each rank.
  */
  Vector<IntVectSet> getIrregularTags() const;

  /*!
    @brief Get the "proxy" grids in AmrMesh.
    @details Returns m_grids
  */
  const Vector<DisjointBoxLayout>& getProxyGrids() const;

  /*!
    @brief Get the grids
  */
  const Vector<DisjointBoxLayout>& getGrids(const std::string a_realm) const;

  /*!
    @brief Get domains
  */
  const Vector<ProblemDomain>& getDomains() const;

  /*!
    @brief Get a registered mask.
    @param[in] a_mask   Mask name
    @param[in] a_buffer Buffer (used when registering mask)
    @param[in] a_realm  Realm name
  */
  const AMRMask& getMask(const std::string a_mask, const int a_buffer, const std::string a_realm) const;

  /*!
    @brief Get the EBLevelGrid for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const Vector<RefCountedPtr<EBLevelGrid> >& getEBLevelGrid(const std::string a_realm, const phase::which_phase a_phase) const;
  
  /*!
    @brief Get EBISLayouts for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const Vector<EBISLayout>& getEBISLayout(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get EBISLayouts for a Realm
    @param[in] a_realm Realm name
  */
  const Vector<RefCountedPtr<MFLevelGrid> >& getMFLevelGrid(const std::string a_realm) const;

  /*!
    @brief Get vof iterators for a Realm and phase. This has the capability of iterating through cut-cells. 
    @details Not const because we need to reset the vofiterator
  */
  Vector<RefCountedPtr<LayoutData<VoFIterator> > >& getVofIterator(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get neighboring boxes
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& getNeighbors(const std::string        a_realm,
										const phase::which_phase a_phase) const;

  /*!
    @brief Get levelset function, allocated over a grid for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const EBAMRFAB& getLevelset(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get the coarsening utility
  */
  Vector<RefCountedPtr<ebcoarseaverage> >& getCoarseAverage(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the ghost cloud particle magic
  */
  Vector<RefCountedPtr<EbGhostCloud> >& getGhostCloud(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the ghost cell interpolation utility
  */
  Vector<RefCountedPtr<NwoEbQuadCfInterp> >& getNWOEBQuadCFInterp(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the old interpolation utility (needed for EBConductivityOps)
  */
  Vector<RefCountedPtr<EBQuadCFInterp> >& getEBQuadCFInterp(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get piecewise linear ghost cell interpolator
  */
  Vector<RefCountedPtr<AggEBPWLFillPatch> >& getFillPatch(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get interpolator 
  */
  Vector<RefCountedPtr<EBPWLFineInterp> >& getPwlInterpolator(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get interpolator 
  */
  Vector<RefCountedPtr<EBMGInterp> >& getEBMGInterp(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get flux register
  */
  Vector<RefCountedPtr<EBFluxRegister> >&  getFluxRegister(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get redistribution
  */
  Vector<RefCountedPtr<EBLevelRedist> >&  getLevelRedist(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get coar-to-fine redistribution for a Realm and phase
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBCoarToFineRedist> >&  getCoarToFineRedist(const std::string        a_realm,
								   const phase::which_phase a_phase);
  
  /*!
    @brief Get coar-to-coar redistribution objects
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBCoarToCoarRedist> >&  getCoarToCoarRedist(const std::string        a_realm,
								   const phase::which_phase a_phase);

  /*!
    @brief Get coar-to-coar redistribution redistribution objects
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBFineToCoarRedist> >&  getFineToCoarRedist(const std::string        a_realm,
								   const phase::which_phase a_phase);
  
  /*!
    @brief Get interpolation stencils for going from cell center to cell centroid
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const IrregAmrStencil<CentroidInterpolationStencil>& getCentroidInterpolationStencils(const std::string a_realm,
									     const phase::which_phase a_phase) const;

  /*!
    @brief Get interpolation stencils for going from cell center to cell boundary centroid
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const IrregAmrStencil<EbCentroidInterpolationStencil>& getEbCentroidInterpolationStencilStencils(const std::string        a_realm,
										  const phase::which_phase a_phase) const;

  /*!
    @brief Get the nonconservative divergence stencils
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  const IrregAmrStencil<NonConservativeDivergenceStencil>& getNonConservativeDivergenceStencils(const std::string a_realm, const phase::which_phase a_phase) const;

  /*!
    @brief Get copier
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<Copier> >& getCopier(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get reverse copier
    @param[in] a_realm Realm name
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<Copier> >& getReverseCopier(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the name of all Realms
    @return Names of all Realms
  */
  std::vector<std::string> getRealms() const;

  /*!
    @brief Get box sorting method
    @return Box sorting method.
  */
  BoxSorting getBoxSorting() const;
  
protected:

  /*!
    @brief Enum for having understandable notation for grid generation
  */
  enum class GridGenerationMethod{
    BergerRigoutsous,
    Tiled,
  };

  /*!
    @brief These are all the Realms
  */
  mutable std::map<std::string, RefCountedPtr<Realm> > m_realms;

  /*!
    @brief Implicit functions
  */
  std::map<phase::which_phase, RefCountedPtr<BaseIF> > m_baseif;

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
    @brief Interpolation type
  */
  GhostInterpolation m_ghostCellInterpolationMethod;

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
  RealVect m_prob_hi;

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
    @brief How much to grow irregular cells
  */
  int m_irregTagGrowth;

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
    @brief Grids have EB x CF crossing
  */
  bool m_hasEbCf;

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
  void defineRealms();
};

#include <CD_NamespaceFooter.H>

#include <CD_AmrMeshImplem.H>

#endif
