/*!
  @file amr_mesh.H
  @brief Declaration of a class that holds grids and useful EBAMR operators (for both phases)
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _AMR_MESH_
#define _AMR_MESH_

#include "type_definitions.H"
#include "computational_geometry.H"
#include "mfis.H"
#include "irreg_amr_stencil.H"
#include "load_balance.H"
#include "nwoebquadcfinterp.H"
#include "MFLevelGrid.H"
#include "ebcoarseaverage.H"
#include "DomainFluxIFFAB.H"
#include "DomainFluxIFFABFactory.H"
#include "EBLevelConcentrationRedist.H"
#include "particle_container.H"
#include "EBGhostCloud.H"
#include "noncons_div.H"
#include "mf_realm.H"

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

/*!
  @brief Class for handling spatial operations for streamer equations.
  @details This class handles most (if not all) spatial operations. Examples are grid generation, coarsening of data, interpolation of data, memory allocations etc. 
*/
class amr_mesh {
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
    @brief Deallocate data. You should never need this, smart pointers are you friends :)
  */
  template<typename T> void deallocate(Vector<T*>& a_data);

  /*!
    @brief Deallocate data
  */
  template<typename T> void deallocate(Vector<RefCountedPtr<T> >& a_data);

  /*!
    @brief Turn smart-pointer data structure into regular-pointer data structure
  */
  template<typename T> void alias(Vector<T*>& a_alias, const Vector<RefCountedPtr<T> >& a_data);

  /*!
    @brief Template class for generic allocation of particles
  */
  template <typename T > void allocate(Vector<RefCountedPtr<ParticleData<T> > >& a_particles);

  /*!
    @brief Allocate a particle container
  */
  template <typename T> void allocate(particle_container<T>& a_container, const int a_buffer);

  /*!
    @brief Allocate pointer but not memory
  */
  template <typename T> void allocate_ptr(Vector<RefCountedPtr<T> >& a_data);

  /*!
    @brief Query operator
  */
  virtual bool query_operator(const std::string a_operator, const phase::which_phase a_phase);

  /*!
    @brief Register an operator
  */
  virtual void register_operator(const std::string a_operator, const phase::which_phase a_phase);

  /*!
    @brief Query if a realm exists
  */
  virtual bool query_realm(const std::string a_realm);

  /*!
    @brief Register a new realm
  */
  virtual void register_realm(const std::string a_realm);

  /*!
    @brief Alias an entire hierarchy of data. 
  */
  virtual void alias(EBAMRCellData&           a_data,
		     const phase::which_phase a_phase,
		     const MFAMRCellData&     a_mfdata,
		     const int                a_finest_level);

  /*!
    @brief Alias an entire hierarchy of data
  */
  virtual void alias(EBAMRCellData& a_data, const phase::which_phase a_phase, const MFAMRCellData& a_mfdata);

  /*!
    @brief Alias an entire hierarchy of data
  */
  virtual void alias(EBAMRFluxData& a_data, const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata);

    /*!
    @brief Alias an entire hierarchy of data
  */
  virtual void alias(EBAMRIVData& a_data, const phase::which_phase a_phase, const MFAMRIVData& a_mfdata);

  /*!
    @brief Compute gradient on level
  */
  virtual void compute_gradient(LevelData<EBCellFAB>& a_gradient,
				const LevelData<EBCellFAB>& a_phi,
				const phase::which_phase a_phase,
				const int a_lvl);

  /*!
    @brief Compute cell-centered gradient. 
    @details User must make sure ghost cells are updated when entering this routine. On output, the result is not averaged down, 
    nor are ghost cells updated. 
  */
  virtual void compute_gradient(EBAMRCellData& a_gradient, const EBAMRCellData& a_phi, const phase::which_phase a_phase);

    /*!
    @brief Compute cell-centered gradient. 
    @details User must make sure ghost cells are updated when entering this routine. On output, the result is not averaged down, 
    nor are ghost cells updated. 
  */
  virtual void compute_gradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi);

  /*!
    @brief Allocate particle valid regions
  */
  virtual void allocate(AMRPVR& a_pvr, const int a_buffer);

  /*!
    @brief Allocate a data holder. 
  */
  virtual void allocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Allocate a data holder. 
  */
  virtual void allocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Allocate a data holder. 
  */
  virtual void allocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Allocate a data holder. 
  */
  virtual void allocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost = -1);

  /*!
    @Allocate Boolean data
  */
  virtual void allocate(EBAMRBool& a_data, const int a_ncomp, const int a_ghost);

  /*!
    @brief Allocate a data holder. 
    @note This routines _allocates_, so if a_data is not empty you're leaking memory
  */
  virtual void allocate(MFAMRCellData& a_data, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Allocate a data holder. 
    @note This routines _allocates_, so if a_data is not empty you're leaking memory
  */
  virtual void allocate(MFAMRFluxData& a_data, const int a_ncomp, const int a_ghost = -1);
  
  /*!
    @brief Allocate a data holder. 
    @note This routines _allocates_, so if a_data is not empty you're leaking memory
  */
  virtual void allocate(MFAMRIVData& a_data, const int a_ncomp, const int a_ghost = -1);
  
  /*!
    @brief Reallocate data. 
    @details This does not touch levels below or equal to a_base. If a_base is < 0, everything is reallocated. 
    Size of ghost vectors and components is maintained. 
  */
  virtual void reallocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data. 
  */
  virtual void reallocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data above a_base
  */
  virtual void reallocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data above a_base
  */
  virtual void reallocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_lmin);

  /*!
    @brief Reallocate data above a_base
  */
  virtual void reallocate(EBAMRBool& a_data, const int a_ncomp, const int a_lmin);

  /*!
    @brief Reallocate data above a_base
  */
  virtual void reallocate(MFAMRCellData& a_data, const int a_lmin);

  /*!
    @brief Reallocate data above a_base
  */
  virtual void reallocate(MFAMRFluxData& a_data, const int a_lmin);
  
  /*!
    @brief Reallocate data above a_base
  */
  virtual void reallocate(MFAMRIVData& a_data, const int a_lmin);

  /*!
    @brief Allocate a data holder over the intersecting regions of the mfis
    @note This routines _allocates_, so if a_data is not empty you're leaking memory
  */
  virtual void allocate_interface(EBAMRIVData& a_data, phase::which_phase a_phase, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Average down both phases
  */
  virtual void average_down(MFAMRCellData& a_data);

  /*!
    @brief Average down
  */
  virtual void average_down(EBAMRCellData& a_data, phase::which_phase a_phase);

  /*!
    @brief Average down level
  */
  virtual void average_down(EBAMRCellData& a_data, phase::which_phase a_phase, const int a_lvl);

  /*!
    @brief Average down both phases
  */
  virtual void average_down(MFAMRFluxData& a_data);

  /*!
    @brief Average down
  */
  virtual void average_down(EBAMRFluxData& a_data, phase::which_phase a_phase);

  /*!
    @brief Average down
  */
  virtual void average_down(EBAMRIVData& a_data, phase::which_phase a_phase);

  /*!
    @brief Average down
  */
  virtual void conservative_average(EBAMRIVData& a_data, phase::which_phase a_phase);

  /*!
    @brief Build new grids. This also load balances. 
  */
  virtual void build_grids(Vector<IntVectSet>& a_tags, const int a_lmin, const int a_lmax, const int a_hardcap = -1);

  /*!
    @brief Interpolate ghost vectors
  */
  virtual void interp_ghost(EBAMRCellData& a_data, phase::which_phase a_phase);

  /*!
    @brief Interpolate ghost cells
  */
  virtual void interp_ghost(MFAMRCellData& a_data);

  /*!
    @brief Fill ghost cells using quadratic interpolation
  */
  virtual void interp_ghost_quad(EBAMRCellData& a_data, phase::which_phase a_phase);

  /*!
    @brief Fill ghost cells using piecewise linear interpolation
  */
  virtual void interp_ghost_pwl(EBAMRCellData& a_data, phase::which_phase a_phase);

  /*!
    @brief Interpolate data to centroids
  */
  virtual void interpolate_to_centroids(EBAMRCellData& a_data, phase::which_phase a_phase);
  
  /*!
    @brief Set mfis
  */
  virtual void set_mfis(const RefCountedPtr<mfis>& a_mfis);

  /*!
    @brief Set load balance
  */
  virtual void parse_load_balance();

  /*!
    @brief Parse domain
  */
  virtual void parse_domain();

  /*!
    @brief Set ghost cell interpolation type. Either piecewise linear or quadratic. 
  */
  virtual void parse_ghost_interpolation();

  /*!
    @brief Build domains
  */
  virtual void build_domains();

  /*!
    @brief Parse the grid generation algorithm
  */
  virtual void parse_grid_generation();

  /*!
    @brief Set verbosity
  */
  virtual void parse_verbosity();

  /*!
    @brief Set coarsest level box; this is the coarsest domain (well, except for lower levels GMG)
  */
  virtual void parse_coarsest_num_cells();

  /*!
    @brief Maximum AMR depth
  */
  virtual void parse_max_amr_depth();

  /*!
    @brief Set maximum simulation depth. 
  */
  virtual void parse_max_simulation_depth();

  /*!
    @brief Refine all cells down to this level
  */
  virtual void parse_refine_all_depth();

  /*!
    @brief Tell amr_mesh to forget about ebcf (user must guarantee that this is the case)
  */
  virtual void parse_ebcf();

  /*!
    @brief Set refinement ratio
  */
  virtual void parse_refinement_ratio();

  /*!
    @brief Set refinement ratios
  */
  virtual void set_refinement_ratios(const Vector<int> a_ref_ratios);

  /*!
    @brief Set buffer size
  */
  virtual void parse_buffer_size();

  /*!
    @brief Set irregular growth
  */
  virtual void parse_irreg_growth();

  /*!
    @brief Set fill ratio
  */
  virtual void parse_fill_ratio();

  /*!
    @brief Set the finest level
  */
  virtual void set_finest_level(const int a_finest_level);



  /*!
    @brief Set the maximum box size
  */
  virtual void parse_max_box_size();

  /*!
    @brief Set the maximum EBIS generation box size
  */
  virtual void parse_max_ebis_box_size();

  /*!
    @brief Set the blocking factor
  */
  virtual void parse_blocking_factor();

  /*!
    @brief Set number of ghost cells for eb stuff
  */ 
  virtual void parse_eb_ghost();

  /*!
    @brief Set number of ghost cells 
  */ 
  virtual void parse_num_ghost();

  /*!
    @brief Set redistribution radius
  */
  virtual void parse_redist_rad();

  /*!
    @brief Parse centroid interpolation stencils
*/
  virtual void parse_centroid_stencils();

  /*!
    @brief Parse EB interpolation (or extrapolation) stencils
  */
  virtual void parse_eb_stencils();

  /*!
    @brief Set preferred stencil type
  */
  virtual void set_irreg_sten_type(const stencil_type::which_type a_type);

  /*!
    @brief Set irregular stencil order
  */
  virtual void set_irreg_sten_order(const int a_irreg_sten_order);

  /*!
    @brief Set maximum stencil radius for irregular stencils
  */
  virtual void set_irreg_sten_radius(const int a_irreg_sten_radius);

  /*!
    @brief Regrid function. Generate new grids. 
    @details lmin is the finest level that changed. 
  */
  virtual void regrid(const Vector<IntVectSet>& a_tags,
		      const int a_lmin,
		      const int a_lmax,
		      const int a_regsize=1, 
		      const int a_hardcap=123456);

  /*!
    @brief Regrid AMR. This builds the grids and sets them. 
  */
  virtual void regrid_amr(const Vector<IntVectSet>& a_tags,
			  const int a_lmin,
			  const int a_lmax,
			  const int a_hardcap=123456);

  /*!
    @brief Regrid AMR. This takes a set of grids. 
  */
  virtual void regrid_amr(const Vector<Vector<int> >& a_procs, const Vector<Vector<Box> >& a_boxes, const int a_lmin);

  /*!
    @brief Set grids on level
  */
  virtual void set_grids(Vector<Vector<Box> >& boxes, const int a_regsize=1);


  /*!
    @brief Regrid AMR operators
  */
  virtual void regrid_operators(const int a_lmin,
				const int a_lmax,
				const int a_regsize=1);



  /*!
    @brief Do a sanity check to make sure everything is set up correctly
  */
  virtual void sanity_check();

  /*!
    @brief True if there is ebcf crossing
  */
  virtual bool get_ebcf();

    /*!
    @brief Get refinement ratios
  */
  virtual int get_finest_level();

  /*!
    @brief Get irregular cell growth
  */
  virtual int get_irreg_growth();

  /*!
    @brief Get maximum amr depth
  */
  virtual int get_max_amr_depth();

  /*!
    @brief Get maximum simulation depth
  */
  virtual int get_max_sim_depth();

  /*!
    @brief Get maximum  depth
  */
  virtual int get_refine_all_depth();

  /*!
    @brief Get blocking factor
  */
  virtual int get_blocking_factor();

  /*!
    @brief Get maximum box size
  */
  virtual int get_max_box_size();

  /*!
    @brief Return bffer
  */
  virtual int get_buffer();

  /*!
    @brief Get maximum box size for EBIS generation
  */
  virtual int get_max_ebis_box_size();

  /*!
    @brief Get number of ghost cells
  */
  virtual int get_num_ghost();

  /*!
    @brief Get number of ghost cells used for EB
  */
  virtual int get_eb_ghost();

  /*!
    @brief Get redistribution radius
  */
  virtual int get_redist_rad();
  
  /*!
    @brief Get finest level
  */
  virtual Real get_finest_dx();

  /*!
    @brief Get origin
  */
  virtual RealVect get_prob_lo();

  /*!
    @brief Get high end
  */
  virtual RealVect get_prob_hi();

  /*!
    @brief Get finest domain
  */
  virtual ProblemDomain get_finest_domain();

  /*!
    @brief Get resolutions
  */
  virtual Vector<Real>& get_dx();

  /*!
    @brief Get refinement ratios
  */
  virtual Vector<int>& get_ref_rat();

  /*!
    @brief Get refinement factor between two levels
  */
  virtual int get_refinement_ratio(const int a_level1, const int a_level2);

  /*!
    @brief Get ALL irregular tags
  */
  virtual Vector<IntVectSet> get_irreg_tags() const;

  /*!
    @brief Get the grids
  */
  virtual Vector<DisjointBoxLayout>& get_grids();

  /*!
    @brief Get the grids
  */
  virtual Vector<DisjointBoxLayout>& get_grids(const std::string a_realm);

  /*!
    @brief Get domains
  */
  virtual Vector<ProblemDomain>& get_domains();

  /*!
    @brief Get the EBLevelGrid 
  */
  virtual Vector<RefCountedPtr<EBLevelGrid> >& get_eblg(phase::which_phase a_phase);

  /*!
    @brief Get the EBLevelGrid 
  */
  virtual Vector<RefCountedPtr<MFLevelGrid> >& get_mflg();

  /*!
    @brief Get the EBISLayouts
  */
  virtual Vector<EBISLayout>& get_ebisl(phase::which_phase a_phase);

  /*!
    @brief Get the vof iterator
  */
  virtual Vector<RefCountedPtr<LayoutData<VoFIterator> > > get_vofit(phase::which_phase a_phase);

  /*!
    @brief Get neighbor boxes
  */
  virtual Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& get_neighbors(phase::which_phase a_phase);

  /*!
    @brief Get the coarsening utility
  */
  virtual Vector<RefCountedPtr<ebcoarseaverage> >& get_coarave(phase::which_phase a_phase);

  /*!
    @brief Get the ghost cloud particle magic
  */
  virtual Vector<RefCountedPtr<EBGhostCloud> >& get_ghostcloud(phase::which_phase a_phase);

  /*!
    @brief Get the ghost cell interpolation utility
  */
  virtual Vector<RefCountedPtr<nwoebquadcfinterp> >& get_quadcfi(phase::which_phase a_phase);

  /*!
    @brief Get the old interpolation utility (needed for EBConductivityOps)
  */
  virtual Vector<RefCountedPtr<EBQuadCFInterp> >& get_old_quadcfi(phase::which_phase a_phase);

  /*!
    @brief Get piecewise linear ghost cell interpolator
  */
  virtual Vector<RefCountedPtr<AggEBPWLFillPatch> >& get_fillpatch(phase::which_phase a_phase);

  /*!
    @brief Get interpolator 
  */
  virtual Vector<RefCountedPtr<EBPWLFineInterp> >& get_eb_pwl_interp(phase::which_phase a_phase);

  /*!
    @brief Get interpolator 
  */
  virtual Vector<RefCountedPtr<EBMGInterp> >& get_eb_mg_interp(phase::which_phase a_phase);

  /*!
    @brief Get flux register
  */
  virtual Vector<RefCountedPtr<EBFluxRegister> >&  get_flux_reg(phase::which_phase a_phase);

  /*!
    @brief Get redistribution
  */
  virtual Vector<RefCountedPtr<EBLevelRedist> >&  get_level_redist(phase::which_phase a_phase);

  /*!
    @brief Get coar-to-fine redistribution
  */
  virtual Vector<RefCountedPtr<EBCoarToFineRedist> >&  get_coar_to_fine_redist(phase::which_phase a_phase);

  /*!
    @brief Get coar-to-coar redistribution
  */
  virtual Vector<RefCountedPtr<EBCoarToCoarRedist> >&  get_coar_to_coar_redist(phase::which_phase a_phase);

  /*!
    @brief Get coar-to-coar redistribution
  */
  virtual Vector<RefCountedPtr<EBFineToCoarRedist> >&  get_fine_to_coar_redist(phase::which_phase a_phase);

  /*!
    @brief Get interpolation stencils
  */
  virtual irreg_amr_stencil<centroid_interp>& get_centroid_interp_stencils(phase::which_phase a_phase);

  /*!
    @brief Get interpolation stencils
  */
  virtual irreg_amr_stencil<eb_centroid_interp>& get_eb_centroid_interp_stencils(phase::which_phase a_phase);

  /*!
    @brief Get Nonconservative divergence stencils
  */
  virtual irreg_amr_stencil<noncons_div>& get_noncons_div_stencils(phase::which_phase a_phase);

  /*!
    @brief Get copier
  */
  virtual Vector<RefCountedPtr<Copier> >& get_copier(phase::which_phase a_phase);

  /*!
    @brief Get reverse copier
  */
  virtual Vector<RefCountedPtr<Copier> >& get_reverse_copier(phase::which_phase a_phase);

  /*!
    @brief Build tiles from a box
  */
  virtual Vector<Box> make_tiles(const Box a_box, const IntVect a_tilesize);

  
protected:

  /*!
    @brief Enum for having understandable notation for grid generation
  */
  enum grid_generation{
    berger_rigoustous = 0,
    tiled             = 1,
  };

  RefCountedPtr<mf_realm> m_realm;

  /*!
    @brief These are all the realms
  */
  std::map<std::string, RefCountedPtr<mf_realm> > m_realms;

  /*!
    @brief Grid generation method
  */
  grid_generation m_gridgen;

  /*!
    @brief mfis
  */
  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief Load balance type
  */
  load_balance::which_balance m_which_balance;

  /*!
    @brief Stencil type
  */
  stencil_type::which_type m_stencil_type;

  /*!
    @brief Stencil type for centroid interpolation
  */
  stencil_type::which_type m_centroid_stencil;

  /*!
    @brief Stencil type for EB interpolation
  */
  stencil_type::which_type m_eb_stencil;

  /*!
    @brief Interpolation type
  */
  ghost_interpolation::which_type m_interp_type;

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
    @brief Number of ghost cells to interpolate
  */
  int m_num_ghost;

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
    @brief Load balance grids
  */
  virtual void loadbalance(Vector<Vector<int> >& a_procs, Vector<Vector<Box> >& a_boxes);

  /*!
    @brief Define realms
  */
  virtual void define_realms();
};

#include "amr_meshI.H"

#endif
