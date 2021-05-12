/*!
  @file amr_mesh.H
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
class amr_mesh {
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
  amr_mesh();

  /*!
    @brief Destructor
  */
  ~amr_mesh();

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
  void allocate_ptr(Vector<RefCountedPtr<T> >& a_data);

    /*!
    @brief Allocate pointer but not memory
  */
  template <typename T>
  void allocate_ptr(Vector<RefCountedPtr<T> >& a_data, const int a_finestLevel);

  /*!
    @brief Allocate pointer but not memory
  */
  template <typename T>
  void allocate_ptr(EBAMRData<T>& a_data);

    /*!
    @brief Allocate pointer but not memory
  */
  template <typename T>
  void allocate_ptr(EBAMRData<T>& a_data, const int a_finestLevel);

  /*!
    @brief Register an operator over a realm. This if for operators that are defined for a specific phase only. 
  */
  void register_operator(const std::string a_operator, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Register a particle halo on a realm. 
  */
  void register_mask(const std::string a_mask, const int a_buffer, const std::string a_realm);

  /*!
    @brief Query if a realm exists
  */
  bool query_realm(const std::string a_realm);

  /*!
    @brief Register a new realm
  */
  void register_realm(const std::string a_realm);

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
  void compute_gradient(LevelData<EBCellFAB>&       a_gradient,
			const LevelData<EBCellFAB>& a_phi,
			const std::string           a_realm,
			const phase::which_phase    a_phase,
			const int                   a_lvl);
  /*!
    @brief Compute cell-centered gradient. 
    @details User must make sure ghost cells are updated when entering this routine. On output, the result is not averaged down, 
    nor are ghost cells updated. 
  */
  void compute_gradient(EBAMRCellData&           a_gradient,
			const EBAMRCellData&     a_phi,
			const std::string        a_realm,  
			const phase::which_phase a_phase);

  /*!
    @brief Compute cell-centered gradient. 
    @details User must make sure ghost cells are updated when entering this routine. On output, the result is not averaged down, 
    nor are ghost cells updated. 
  */
  void compute_gradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi, const std::string a_realm);

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
  void average_down(MFAMRCellData& a_data, const std::string a_realm);

  /*!
    @brief Average down on specific realm
  */
  void average_down(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Average down level
  */
  void average_down(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase, const int a_lvl);

  /*!
    @brief Average down both phases on specific realm
  */
  void average_down(MFAMRFluxData& a_data, const std::string a_realm);

  /*!
    @brief Average down on specific realm
  */
  void average_down(EBAMRFluxData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Average down on realm
  */
  void average_down(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Average down on realm
  */
  void conservative_average(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Build new grids. This also load balances. 
  */
  void build_grids(Vector<IntVectSet>& a_tags, const int a_lmin, const int a_lmax, const int a_hardcap = -1);

  /*!
    @brief Interpolate ghost vectors over a realm
  */
  void interp_ghost(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Interpolate ghost cells. Realm, phase, and level version.
  */
  void interp_ghost(LevelData<EBCellFAB>&       a_fineData,
		    const LevelData<EBCellFAB>& a_coarData,
		    const int                   a_level,
		    const std::string           a_realm,
		    const phase::which_phase    a_phase);

  /*!
    @brief Interpolate ghost cells over a realm
  */
  void interp_ghost(MFAMRCellData& a_data, const std::string a_realm);

  /*!
    @brief Fill ghost cells using quadratic interpolation
  */
  void interp_ghost_quad(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Fill ghost cells using piecewise linear interpolation. Specific realm. 
  */
  void interp_ghost_pwl(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Interpolate data to centroids on realm
  */
  void interpolate_to_centroids(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase);
  
  /*!
    @brief Set mfis
  */
  void set_mfis(const RefCountedPtr<mfis>& a_mfis);

  /*!
    @brief Set implicit function
  */
  void set_baseif(const phase::which_phase a_phase, const RefCountedPtr<BaseIF>& a_baseif);

  /*!
    @brief Parse domain
  */
  void parse_domain();

  /*!
    @brief Set ghost cell interpolation type. Either piecewise linear or quadratic. 
  */
  void parse_ghost_interpolation();

  /*!
    @brief Build domains
  */
  void build_domains();

  /*!
    @brief Parse the grid generation algorithm
  */
  void parse_grid_generation();

  /*!
    @brief Set verbosity
  */
  void parse_verbosity();

  /*!
    @brief Set coarsest level box; this is the coarsest domain (well, except for lower levels GMG)
  */
  void parse_coarsest_num_cells();

  /*!
    @brief Maximum AMR depth
  */
  void parse_max_amr_depth();

  /*!
    @brief Set maximum simulation depth. 
  */
  void parse_max_simulation_depth();

  /*!
    @brief Tell amr_mesh to forget about ebcf (user must guarantee that this is the case)
  */
  void parse_ebcf();

  /*!
    @brief Set refinement ratio
  */
  void parse_refinement_ratio();

  /*!
    @brief Set refinement ratios
  */
  void set_refinement_ratios(const Vector<int> a_ref_ratios);

  /*!
    @brief Set buffer size
  */
  void parse_buffer_size();

  /*!
    @brief Set irregular growth
  */
  void parse_irreg_growth();

  /*!
    @brief Set fill ratio
  */
  void parse_fill_ratio();

  /*!
    @brief Set the finest level
  */
  void set_finest_level(const int a_finest_level);

  /*!
    @brief Set the maximum box size
  */
  void parse_max_box_size();

  /*!
    @brief Set the maximum EBIS generation box size
  */
  void parse_max_ebis_box_size();

  /*!
    @brief Set the blocking factor
  */
  void parse_blocking_factor();

  /*!
    @brief Set number of ghost cells for eb stuff
  */ 
  void parse_eb_ghost();

  /*!
    @brief Set number of ghost cells 
  */ 
  void parse_num_ghost();

  /*!
    @brief Set redistribution radius
  */
  void parse_redist_rad();

  /*!
    @brief Parse centroid interpolation stencils
  */
  void parse_centroid_stencils();

  /*!
    @brief Parse EB interpolation (or extrapolation) stencils
  */
  void parse_eb_stencils();

  /*!
    @brief Set preferred stencil type
  */
  void set_irreg_sten_type(const stencil_type a_type);

  /*!
    @brief Set irregular stencil order
  */
  void set_irreg_sten_order(const int a_irreg_sten_order);

  /*!
    @brief Set maximum stencil radius for irregular stencils
  */
  void set_irreg_sten_radius(const int a_irreg_sten_radius);

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
  void regrid_amr(const Vector<IntVectSet>& a_tags,
		  const int a_lmin,
		  const int a_lmax,
		  const int a_hardcap=123456);

  /*!
    @brief Regrid AMR. This takes a set of grids. 
    @todo Remove. 
  */
  void regrid_amr(const Vector<Vector<int> >& a_procs, const Vector<Vector<Box> >& a_boxes, const int a_lmin);

  /*!
    @brief Define a realm
  */
  void regrid_realm(const std::string a_realm,
		    const Vector<Vector<int> >& a_procs,
		    const Vector<Vector<Box> >& a_boxes,
		    const int a_lmin);

  /*!
    @brief Set grids from boxes and computational loads
  */
  void set_grids(const Vector<Vector<Box> >& a_boxes, const std::map<std::string, Vector<Vector<long int> > >& a_loads);

  /*!
    @brief Regrid AMR operators
  */
  void regrid_operators(const int a_lmin,
			const int a_lmax,
			const int a_regsize=1);



  /*!
    @brief Do a sanity check to make sure everything is set up correctly
  */
  void sanity_check();

  /*!
    @brief True if there is ebcf crossing
  */
  bool get_ebcf();

  /*!
    @brief Get refinement ratios
  */
  int get_finest_level();

  /*!
    @brief Get irregular cell growth
  */
  int get_irreg_growth();

  /*!
    @brief Get maximum amr depth
  */
  int get_max_amr_depth();

  /*!
    @brief Get maximum simulation depth
  */
  int get_max_sim_depth();

  /*!
    @brief Get maximum  depth
  */
  int get_refine_all_depth();

  /*!
    @brief Get blocking factor
  */
  int get_blocking_factor();

  /*!
    @brief Get maximum box size
  */
  int get_max_box_size();

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
  int get_num_ghost();

  /*!
    @brief Get number of ghost cells used for EB
  */
  int get_eb_ghost();

  /*!
    @brief Get redistribution radius
  */
  int get_redist_rad();
  
  /*!
    @brief Get finest level
  */
  Real get_finest_dx();

  /*!
    @brief Get origin
  */
  RealVect get_prob_lo();

  /*!
    @brief Get high end
  */
  RealVect get_prob_hi();

  /*!
    @brief Get finest domain
  */
  ProblemDomain get_finest_domain();

  /*!
    @brief Get resolutions
  */
  Vector<Real>& get_dx();

  /*!
    @brief Get refinement ratios
  */
  Vector<int>& get_ref_rat();

  /*!
    @brief Get refinement factor between two levels
  */
  int get_refinement_ratio(const int a_level1, const int a_level2);

  /*!
    @brief Get implicit function
  */
  const RefCountedPtr<BaseIF>& get_baseif(const phase::which_phase a_phase);

  /*!
    @brief Get ALL irregular tags
  */
  Vector<IntVectSet> get_irreg_tags() const;

  /*!
    @brief Get the grids
  */
  Vector<DisjointBoxLayout>& get_proxy_grids();

  /*!
    @brief Get the grids
  */
  Vector<DisjointBoxLayout>& get_grids(const std::string a_realm);

  /*!
    @brief Get domains
  */
  Vector<ProblemDomain>& get_domains();

  /*!
    @brief Get mask
  */
  AMRMask& get_mask(const std::string a_mask, const int a_buffer, const std::string a_realm);

  /*!
    @brief Get the EBLevelGrid 
  */
  Vector<RefCountedPtr<EBLevelGrid> >& get_eblg(const std::string a_realm, const phase::which_phase a_phase);
  
  /*!
    @brief Get the EBISLayouts
  */
  Vector<EBISLayout>& get_ebisl(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the EBLevelGrid 
  */
  Vector<RefCountedPtr<MFLevelGrid> >& get_mflg(const std::string a_realm);

  /*!
    @brief Get the vof iterator
  */
  Vector<RefCountedPtr<LayoutData<VoFIterator> > > get_vofit(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get neighbor boxes
  */
  Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& get_neighbors(const std::string a_realm,
									   const phase::which_phase a_phase);

  /*!
    @brief Get levelset function
  */
  EBAMRFAB& get_levelset(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the coarsening utility
  */
  Vector<RefCountedPtr<ebcoarseaverage> >& get_coarave(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the ghost cloud particle magic
  */
  Vector<RefCountedPtr<EBGhostCloud> >& get_ghostcloud(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the ghost cell interpolation utility
  */
  Vector<RefCountedPtr<nwoebquadcfinterp> >& get_quadcfi(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get the old interpolation utility (needed for EBConductivityOps)
  */
  Vector<RefCountedPtr<EBQuadCFInterp> >& get_old_quadcfi(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get piecewise linear ghost cell interpolator
  */
  Vector<RefCountedPtr<AggEBPWLFillPatch> >& get_fillpatch(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get interpolator 
  */
  Vector<RefCountedPtr<EBPWLFineInterp> >& get_eb_pwl_interp(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get interpolator 
  */
  Vector<RefCountedPtr<EBMGInterp> >& get_eb_mg_interp(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get flux register
  */
  Vector<RefCountedPtr<EBFluxRegister> >&  get_flux_reg(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get redistribution
  */
  Vector<RefCountedPtr<EBLevelRedist> >&  get_level_redist(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get coar-to-fine redistribution
  */
  Vector<RefCountedPtr<EBCoarToFineRedist> >&  get_coar_to_fine_redist(const std::string        a_realm,
								       const phase::which_phase a_phase);
  
  /*!
    @brief Get coar-to-coar redistribution
  */
  Vector<RefCountedPtr<EBCoarToCoarRedist> >&  get_coar_to_coar_redist(const std::string        a_realm,
								       const phase::which_phase a_phase);

  /*!
    @brief Get coar-to-coar redistribution
  */
  Vector<RefCountedPtr<EBFineToCoarRedist> >&  get_fine_to_coar_redist(const std::string        a_realm,
								       const phase::which_phase a_phase);
  
  /*!
    @brief Get interpolation stencils
  */
  irreg_amr_stencil<centroid_interp>& get_centroid_interp_stencils(const std::string a_realm,
								   const phase::which_phase a_phase);

  /*!
    @brief Get interpolation stencils
  */
  irreg_amr_stencil<eb_centroid_interp>& get_eb_centroid_interp_stencils(const std::string        a_realm,
									 const phase::which_phase a_phase);

  /*!
    @brief Get Nonconservative divergence stencils
  */
  irreg_amr_stencil<noncons_div>& get_noncons_div_stencils(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get copier
  */
  Vector<RefCountedPtr<Copier> >& get_copier(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Get reverse copier
  */
  Vector<RefCountedPtr<Copier> >& get_reverse_copier(const std::string a_realm, const phase::which_phase a_phase);

  /*!
    @brief Build tiles from a box
  */
  Vector<Box> make_tiles(const Box a_box, const IntVect a_tilesize);

  /*!
    @brief Get the name of all realms
  */
  std::vector<std::string> get_realms() const;

  /*!
    @brief Get box sorting method
  */
  box_sorting get_box_sorting() const;
  
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
  grid_generation m_gridgen;

  /*!
    @brief Box sorting
  */
  box_sorting m_boxsort;

  /*!
    @brief mfis
  */
  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief Stencil type
  */
  stencil_type m_stencil_type;

  /*!
    @brief Stencil type for centroid interpolation
  */
  stencil_type m_centroid_stencil;

  /*!
    @brief Stencil type for EB interpolation
  */
  stencil_type m_eb_stencil;

  /*!
    @brief Interpolation type
  */
  ghost_interpolation m_interp_type;

  /*!
    @brief Coarsest box where we compute
  */
  IntVect m_num_cells;

  /*!
    @brief Fill ratio
  */
  Real m_fill_ratio;

  /*!
    @brief Domain simulation corner
  */
  RealVect m_prob_lo;

  /*!
    @brief Domain simulation corner
  */
  RealVect m_prob_hi;

  /*!
    @brief Refinement ratio
  */
  int m_ref_ratio;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Finest level
  */
  int m_finest_level;

  /*!
    @brief Maximum amr depth
  */
  int m_max_amr_depth;

  /*!
    @brief Maximum allowed depth for simulation
  */
  int m_max_sim_depth;

  /*!
    @brief Uniform refinement depth
  */
  int m_refine_all_depth;

  /*!
    @brief Max box size
  */
  int m_max_box_size;

  /*!
    @brief Maximum box size for EBIS generation
  */
  int m_max_ebis_box_size;

  /*!
    @brief Set buffer size
  */
  int m_buffer_size;

  /*!
    @brief How much to grow irregular cells
  */
  int m_irreg_growth;

  /*!
    @brief Blocking factor
  */
  int m_blocking_factor;

  /*!
    @brief Number of ghost cells to use for eb stuff
  */
  int m_ebghost;

  /*!
    @brief Number of ghost cells
  */
  int m_num_ghost;

  /*!
    @brief Number of ghost cells to use when writing level-set to grid
  */
  int m_lsf_ghost;

  /*!
    @brief Redistribution radius
  */
  int m_redist_rad;

  /*!
    @brief Order for interpolation stencils to centroids
  */
  int m_centroid_sten_order;

  /*!
    @brief Radius for centroid interpolation stencil
  */
  int m_centroid_sten_rad;

  /*!
    @brief Order for interpolation stencilst to EB
  */
  int m_eb_sten_order;

  /*!
    @brief Radius for EB interpolation stencil
  */
  int m_eb_sten_rad;

  /*!
    @brief Stencil order for irregular stencils
  */
  int m_irreg_sten_order;

  /*!
    @brief radius for irregular stencils
  */
  int m_irreg_sten_radius;

  /*!
    @brief Grids have EB x CF crossing
  */
  bool m_ebcf;

  /*!
    @brief Has grids or not
  */
  bool m_has_grids;

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
  Vector<int> m_ref_ratios;

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

#include "amr_meshImplem.H"

#endif
