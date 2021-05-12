/*!
  @file AmrMesh.H
  @brief Declaration of a class that holds grids and useful EBAMR operators (for both phases)
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _AMR_MESH_
#define _AMR_MESH_

#include "EBAMRData.H"
#include "computational_geometry.H"
#include "mfis.H"
#include "irreg_amr_stencil.H"
#include "load_balance.H"
#include "nwoebquadcfinterp.H"
#include "MFLevelGrid.H"
#include "ebcoarseaverage.H"
#include "DomainFluxIFFAB.H"
#include "DomainFluxIFFABFactory.H"
#include "particle_container.H"
#include "EBGhostCloud.H"
#include "noncons_div.H"
#include "realm.H"
#include "load_balance.H"

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


#include "CD_NamespaceHeader.H"
  
/*!
  @brief Class for handling spatial operations for streamer equations.
  @details This class handles most (if not all) spatial operations. Examples are grid generation, coarsening of data, interpolation of data, memory allocations etc. 
*/
class AmrMesh {
protected:
  
  /*!
    @brief Silly enum for default ghost interpolation
  */
  enum class ghost_interpolation {
    pwl,
    quad,
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
    @brief Parse options
  */
  void parse_options();

  /*!
    @brief Parse options
  */
  void parse_runtime_options();

  /*!
    @brief Deallocate data. You should never need this, smart pointers are you friends :)
  */
  template<typename T>
  void deallocate(Vector<T*>& a_data);

  /*!
    @brief Deallocate data
  */
  template<typename T>
  void deallocate(Vector<RefCountedPtr<T> >& a_data);

  /*!
    @brief Deallocate data
  */
  template<typename T>
  void deallocate(EBAMRData<T>& a_data);

  /*!
    @brief Turn smart-pointer data structure into regular-pointer data structure
  */
  template<typename T>
  void alias(Vector<T*>& a_alias, const Vector<RefCountedPtr<T> >& a_data);

  /*!
    @brief Turn smart-pointer data structure into regular-pointer data structure
  */
  template<typename T, typename S>
  void alias(Vector<T*>& a_alias, const EBAMRData<S>& a_data);

  /*!
    @brief Template class for generic allocation of particles
  */
  template <typename T >
  void allocate(Vector<RefCountedPtr<ParticleData<T> > >& a_particles, const std::string a_realm);

  /*!
    @brief Allocate a particle container
  */
  template <typename T>
  void allocate(particle_container<T>& a_container,
		const int              a_pvr_buffer,
		const std::string      a_realm);

  /*!
    @brief Allocate pointer but not memory
  */
  template <typename T>
  void allocatePointer(Vector<RefCountedPtr<T> >& a_data);

    /*!
    @brief Allocate pointer but not memory
  */
  template <typename T>
  void allocatePointer(Vector<RefCountedPtr<T> >& a_data, const int a_finestLevel);

  /*!
    @brief Allocate pointer but not memory
  */
  template <typename T>
  void allocatePointer(EBAMRData<T>& a_data);

    /*!
    @brief Allocate pointer but not memory
  */
  template <typename T>
  void allocatePointer(EBAMRData<T>& a_data, const int a_finestLevel);

  /*!
    @brief Register an operator over a realm. This if for operators that are defined for a specific phase only. 
  */
  void registerOperator(const std::string a_operator, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Register a particle halo on a realm. 
  */
  void registerMask(const std::string a_mask, const int a_buffer, const std::string a_realm);

  /*!
    @brief Query if a realm exists
  */
  bool queryRealm(const std::string a_realm);

  /*!
    @brief Register a new realm
  */
  void registerRealm(const std::string a_realm);

  /*!
    @brief Alias an entire hierarchy of data. 
  */
  void alias(EBAMRCellData&           a_data,
	     const phase::which_phase a_phase,
	     const MFAMRCellData&     a_mfdata,
	     const int                a_finestLevel);

    /*!
    @brief Alias an entire hierarchy of data. 
  */
  void alias(EBAMRFluxData&           a_data,
	     const phase::which_phase a_phase,
	     const MFAMRFluxData&     a_mfdata,
	     const int                a_finestLevel);

  /*!
    @brief Alias an entire hierarchy of data
  */
  void alias(EBAMRCellData& a_data, const phase::which_phase a_phase, const MFAMRCellData& a_mfdata);

  /*!
    @brief Alias an entire hierarchy of data
  */
  void alias(EBAMRFluxData& a_data, const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata);

  /*!
    @brief Alias an entire hierarchy of data
  */
  void alias(EBAMRIVData& a_data, const phase::which_phase a_phase, const MFAMRIVData& a_mfdata);

  /*!
    @brief Alias data on phase. 
  */
  EBAMRCellData alias(const phase::which_phase a_phase, const MFAMRCellData& a_mfdata);

  /*!
    @brief Alias data on phase. 
  */
  EBAMRFluxData alias(const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata);

  /*!
    @brief Alias data on phase. 
  */
  EBAMRIVData alias(const phase::which_phase a_phase, const MFAMRIVData& a_mfdata);

  /*!
    @brief Compute gradient on level
  */
  void computeGradient(LevelData<EBCellFAB>&       a_gradient,
			const LevelData<EBCellFAB>& a_phi,
			const std::string           a_realm,
			const phase::which_phase    a_phase,
			const int                   a_lvl);
  /*!
    @brief Compute cell-centered gradient. 
    @details User must make sure ghost cells are updated when entering this routine. On output, the result is not averaged down, 
    nor are ghost cells updated. 
  */
  void computeGradient(EBAMRCellData&           a_gradient,
			const EBAMRCellData&     a_phi,
			const std::string        a_realm,  
			const phase::which_phase a_phase);

  /*!
    @brief Compute cell-centered gradient. 
    @details User must make sure ghost cells are updated when entering this routine. On output, the result is not averaged down, 
    nor are ghost cells updated. 
  */
  void computeGradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi, const std::string a_realm);

  /*!
    @brief Allocate a data holder over a specific realm
  */
  void allocate(EBAMRCellData&           a_data,
		const std::string        a_realm,
		const phase::which_phase a_phase,
		const int                a_ncomp,
		const int                a_ghost = -1);

  /*!
    @brief Allocate a data holder over a specific realm
  */
  void allocate(EBAMRFluxData&           a_data,
		const std::string        a_realm,
		const phase::which_phase a_phase,
		const int                a_ncomp,
		const int                a_ghost = -1);

  /*!
    @brief Allocate a data holder over a specific realm
  */
  void allocate(EBAMRIVData&             a_data,
		const std::string        a_realm,
		const phase::which_phase a_phase,
		const int                a_ncomp,
		const int                a_ghost = -1);

  /*!
    @brief Allocate a data holder over a specific realm
  */
  void allocate(EBAMRIFData&             a_data,
		const std::string        a_realm,
		const phase::which_phase a_phase,
		const int                a_ncomp,
		const int                a_ghost = -1);

  /*!
    @brief Allocate Boolean data over a specific realm
  */
  void allocate(EBAMRBool& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost);
  
  /*!
    @brief Allocate a data holder over a specific realm
  */
  void allocate(MFAMRCellData& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Allocate a data holder over a specific realm.
  */
  void allocate(MFAMRFluxData& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Allocate a data holder over a specific realm
  */
  void allocate(MFAMRIVData& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Reallocate data. 
  */
  void reallocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data. 
  */
  void reallocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data above a_base
  */
  void reallocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data above a_base
  */
  void reallocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data above a_base
  */
  void reallocate(EBAMRBool& a_data, const int a_lmin);

  /*!
    @brief Reallocate data above a_base
  */
  void reallocate(MFAMRCellData& a_data, const int a_lmin);

  /*!
    @brief Reallocate data above a_base
  */
  void reallocate(MFAMRFluxData& a_data, const int a_lmin);
  
  /*!
    @brief Reallocate data above a_base
  */
  void reallocate(MFAMRIVData& a_data, const int a_lmin);

  /*!
    @brief Average down both phases on specific realm
  */
  void averageDown(MFAMRCellData& a_data, const std::string a_realm);

  /*!
    @brief Average down on specific realm
  */
  void averageDown(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Average down level
  */
  void averageDown(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase, const int a_lvl);

  /*!
    @brief Average down both phases on specific realm
  */
  void averageDown(MFAMRFluxData& a_data, const std::string a_realm);

  /*!
    @brief Average down on specific realm
  */
  void averageDown(EBAMRFluxData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Average down on realm
  */
  void averageDown(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Average down on realm
  */
  void conservativeAverage(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Build new grids. This also load balances. 
  */
  void buildGrids(Vector<IntVectSet>& a_tags, const int a_lmin, const int a_lmax, const int a_hardcap = -1);

  /*!
    @brief Interpolate ghost vectors over a realm
  */
  void interpGhost(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Interpolate ghost cells. Realm, phase, and level version.
  */
  void interpGhost(LevelData<EBCellFAB>&       a_fineData,
		    const LevelData<EBCellFAB>& a_coarData,
		    const int                   a_level,
		    const std::string           a_realm,
		    const phase::which_phase    a_phase);

  /*!
    @brief Interpolate ghost cells over a realm
  */
  void interpGhost(MFAMRCellData& a_data, const std::string a_realm);

  /*!
    @brief Fill ghost cells using quadratic interpolation
  */
  void interpGhost_quad(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Fill ghost cells using piecewise linear interpolation. Specific realm. 
  */
  void interpGhost_pwl(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Interpolate data to centroids on realm
  */
  void InterpToCentroids(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);
  
  /*!
    @brief Set mfis
  */
  void setMultifluidIndexSpace(const RefCountedPtr<mfis>& a_mfis);

  /*!
    @brief Set implicit function
  */
  void setBaseImplicitFunction(const phase::which_phase a_phase, const RefCountedPtr<BaseIF>& a_baseif);

  /*!
    @brief Parse domain
  */
  void parseDomain();

  /*!
    @brief Set ghost cell interpolation type. Either piecewise linear or quadratic. 
  */
  void parseGhostInterpolation();

  /*!
    @brief Build domains
  */
  void buildDomains();

  /*!
    @brief Parse the grid generation algorithm
  */
  void parsegridGeneration();

  /*!
    @brief Set verbosity
  */
  void parseVerbosity();

  /*!
    @brief Set coarsest level box; this is the coarsest domain (well, except for lower levels GMG)
  */
  void parseCoarsestLevelNumCells();

  /*!
    @brief Maximum AMR depth
  */
  void parseMaxAmrDepth();

  /*!
    @brief Set maximum simulation depth. 
  */
  void parse_max_simulation_depth();

  /*!
    @brief Tell AmrMesh to forget about ebcf (user must guarantee that this is the case)
  */
  void parseEbCf();

  /*!
    @brief Set refinement ratio
  */
  void parseRefinementRatios();

  /*!
    @brief Set refinement ratios
  */
  void setRefinementRatios(const Vector<int> a_ref_ratios);

  /*!
    @brief Set buffer size
  */
  void parseBrBufferSize();

  /*!
    @brief Set irregular growth
  */
  void parseIrregTagGrowth();

  /*!
    @brief Set fill ratio
  */
  void parseBrFillRatio();

  /*!
    @brief Set the finest level
  */
  void setFinestLevel(const int a_finest_level);

  /*!
    @brief Set the maximum box size
  */
  void parseMaxBoxSize();

  /*!
    @brief Set the maximum EBIS generation box size
  */
  void parseMaxEbisBoxSize();

  /*!
    @brief Set the blocking factor
  */
  void parseBlockingFactor();

  /*!
    @brief Set number of ghost cells for eb stuff
  */ 
  void parseEbGhostCells();

  /*!
    @brief Set number of ghost cells 
  */ 
  void parseNumGhostCells();

  /*!
    @brief Set redistribution radius
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
    @brief Set preferred stencil type
  */
  void setIrregularInterpolationStencilType(const stencil_type a_type);

  /*!
    @brief Set irregular stencil order
  */
  void setIrregularInterpolationStencilOrder(const int a_irreg_sten_order);

  /*!
    @brief Set maximum stencil radius for irregular stencils
  */
  void setIrregularInterpolationStencilRadius(const int a_irreg_sten_radius);

  /*!
    @brief Regrid function. Generate new grids. 
    @details lmin is the finest level that changed. 
  */
  void regrid(const Vector<IntVectSet>& a_tags,
	      const int a_lmin,
	      const int a_lmax,
	      const int a_regsize=1, 
	      const int a_hardcap=123456);

  /*!
    @brief Regrid AMR. This builds the grids and sets them. 
  */
  void regridAmr(const Vector<IntVectSet>& a_tags,
		  const int a_lmin,
		  const int a_lmax,
		  const int a_hardcap=123456);

  /*!
    @brief Regrid AMR. This takes a set of grids. 
    @todo Remove. 
  */
  void regridAmr(const Vector<Vector<int> >& a_procs, const Vector<Vector<Box> >& a_boxes, const int a_lmin);

  /*!
    @brief Define a realm
  */
  void regridRealm(const std::string a_realm,
		    const Vector<Vector<int> >& a_procs,
		    const Vector<Vector<Box> >& a_boxes,
		    const int a_lmin);

  /*!
    @brief Set grids from boxes and computational loads
  */
  void setGrids(const Vector<Vector<Box> >& a_boxes, const std::map<std::string, Vector<Vector<long int> > >& a_loads);

  /*!
    @brief Regrid AMR operators
  */
  void regridOperators(const int a_lmin,
			const int a_lmax,
			const int a_regsize=1);



  /*!
    @brief Do a sanity check to make sure everything is set up correctly
  */
  void sanityCheck();

  /*!
    @brief True if there is ebcf crossing
  */
  bool getEbCf();

  /*!
    @brief Get refinement ratios
  */
  int getFinestLevel();

  /*!
    @brief Get irregular cell growth
  */
  int getIrregTagGrowth();

  /*!
    @brief Get maximum amr depth
  */
  int getMaxAmrDepth();

  /*!
    @brief Get maximum simulation depth
  */
  int getMaxSimulationDepth();

  /*!
    @brief Get maximum  depth
  */
  int getRefineAllLevelsDepth();

  /*!
    @brief Get blocking factor
  */
  int getBlockingFactor();

  /*!
    @brief Get maximum box size
  */
  int getMaxBoxSize();

  /*!
    @brief Return bffer
  */
  int get_buffer();

  /*!
    @brief Get maximum box size for EBIS generation
  */
  int get_max_ebis_box_size();

  /*!
    @brief Get number of ghost cells
  */
  int getNumberOfGhostCells();

  /*!
    @brief Get number of ghost cells used for EB
  */
  int getNumberOfEbGhostCells();

  /*!
    @brief Get redistribution radius
  */
  int getRedistributionRadius();
  
  /*!
    @brief Get finest level
  */
  Real getFinestDx();

  /*!
    @brief Get origin
  */
  RealVect getProbLo();

  /*!
    @brief Get high end
  */
  RealVect getProbHi();

  /*!
    @brief Get finest domain
  */
  ProblemDomain getFinestDomain();

  /*!
    @brief Get resolutions
  */
  Vector<Real>& getDx();

  /*!
    @brief Get refinement ratios
  */
  Vector<int>& getRefinementRatios();

  /*!
    @brief Get refinement factor between two levels
  */
  int getRefinementRatio(const int a_level1, const int a_level2);

  /*!
    @brief Get implicit function
  */
  const RefCountedPtr<BaseIF>& getBaseImplicitFunction(const phase::which_phase a_phase);

  /*!
    @brief Get ALL irregular tags
  */
  Vector<IntVectSet> getIrregularTags() const;

  /*!
    @brief Get the grids
  */
  Vector<DisjointBoxLayout>& getProxyGrids();

  /*!
    @brief Get the grids
  */
  Vector<DisjointBoxLayout>& getGrids(const std::string a_realm);

  /*!
    @brief Get domains
  */
  Vector<ProblemDomain>& getDomains();

  /*!
    @brief Get mask
  */
  AMRMask& getMask(const std::string a_mask, const int a_buffer, const std::string a_realm);

  /*!
    @brief Get the EBLevelGrid 
  */
  Vector<RefCountedPtr<EBLevelGrid> >& getEBLevelGrid(const std::string a_realm, const phase::which_phase a_phase);
  
  /*!
    @brief Get the EBISLayouts
  */
  Vector<EBISLayout>& getEBISLayout(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the EBLevelGrid 
  */
  Vector<RefCountedPtr<MFLevelGrid> >& getMFLevelGrid(const std::string a_realm);

  /*!
    @brief Get the vof iterator
  */
  Vector<RefCountedPtr<LayoutData<VoFIterator> > > getVofIterator(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get neighbor boxes
  */
  Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& getNeighbors(const std::string a_realm,
									   const phase::which_phase a_phase);

  /*!
    @brief Get levelset function
  */
  EBAMRFAB& getLevelset(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the coarsening utility
  */
  Vector<RefCountedPtr<ebcoarseaverage> >& getCoarseAverage(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the ghost cloud particle magic
  */
  Vector<RefCountedPtr<EBGhostCloud> >& getGhostCloud(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the ghost cell interpolation utility
  */
  Vector<RefCountedPtr<nwoebquadcfinterp> >& getNWOEBQuadCFInterp(const std::string a_realm, const phase::which_phase a_phase);

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
    @brief Get coar-to-fine redistribution
  */
  Vector<RefCountedPtr<EBCoarToFineRedist> >&  getCoarToFineRedist(const std::string        a_realm,
								       const phase::which_phase a_phase);
  
  /*!
    @brief Get coar-to-coar redistribution
  */
  Vector<RefCountedPtr<EBCoarToCoarRedist> >&  getCoarToCoarRedist(const std::string        a_realm,
								       const phase::which_phase a_phase);

  /*!
    @brief Get coar-to-coar redistribution
  */
  Vector<RefCountedPtr<EBFineToCoarRedist> >&  getFineToCoarRedist(const std::string        a_realm,
								       const phase::which_phase a_phase);
  
  /*!
    @brief Get interpolation stencils
  */
  irreg_amr_stencil<centroid_interp>& getCentroidInterpolationStencils(const std::string a_realm,
								   const phase::which_phase a_phase);

  /*!
    @brief Get interpolation stencils
  */
  irreg_amr_stencil<eb_centroid_interp>& getEbCentroidInterpolationStencils(const std::string        a_realm,
									 const phase::which_phase a_phase);

  /*!
    @brief Get Nonconservative divergence stencils
  */
  irreg_amr_stencil<noncons_div>& getNonConservativeDivergenceStencils(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get copier
  */
  Vector<RefCountedPtr<Copier> >& getCopier(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get reverse copier
  */
  Vector<RefCountedPtr<Copier> >& getReverseCopier(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Build tiles from a box
  */
  Vector<Box> makeTiles(const Box a_box, const IntVect a_tilesize);

  /*!
    @brief Get the name of all realms
  */
  std::vector<std::string> getRealms() const;

  /*!
    @brief Get box sorting method
  */
  box_sorting getBoxSorting() const;
  
protected:

  /*!
    @brief Enum for having understandable notation for grid generation
  */
  enum grid_generation{
    berger_rigoustous = 0,
    tiled             = 1,
  };

  /*!
    @brief These are all the realms
  */
  std::map<std::string, RefCountedPtr<realm> > m_realms;

  /*!
    @brief Implicit functions
  */
  std::map<phase::which_phase, RefCountedPtr<BaseIF> > m_baseif;

  /*!
    @brief Grid generation method
  */
  grid_generation m_gridGenerationMethod;

  /*!
    @brief Box sorting
  */
  box_sorting m_boxSort;

  /*!
    @brief mfis
  */
  RefCountedPtr<mfis> m_multifluidIndexSpace;

  /*!
    @brief Stencil type
  */
  stencil_type m_stencilType;

  /*!
    @brief Stencil type for centroid interpolation
  */
  stencil_type m_centroidStencilType;

  /*!
    @brief Stencil type for EB interpolation
  */
  stencil_type m_ebCentroidStencilType;

  /*!
    @brief Interpolation type
  */
  ghost_interpolation m_ghostCellInterpolationMethod;

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
    @brief Uniform refinement depth
  */
  int m_refineAllLevelsDepth;

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
  int m_centroid_sten_rad;

  /*!
    @brief Order for interpolation stencilst to EB
  */
  int m_ebCentroidStencilOrder;

  /*!
    @brief Radius for EB interpolation stencil
  */
  int m_eb_sten_rad;

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
    @brief Define realms
  */
  void define_realms();
};
#include "CD_NamespaceFooter.H"

#include <CD_AmrMeshImplem.H>

#endif
