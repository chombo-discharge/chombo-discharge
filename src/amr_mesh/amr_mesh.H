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

#include <DisjointBoxLayout.H>
#include <AggEBPWLFillPatch.H>
#include <ebcoarseaverage.H>
#include <EBFluxRegister.H>
#include <EBFluxRegister.H>
#include <EBLevelRedist.H>
#include <EBCoarToFineRedist.H>
#include <EBCoarToCoarRedist.H>
#include <EBFineToCoarRedist.H>
#include <EBPWLFineInterp.H>
#include <ProblemDomain.H>

/*!
  @brief Class for handling spatial operations for streamer equations.
  @details This class handles most (if not all) spatial operations. Examples are grid generation, coarsening of data, interpolation of data, memory allocations etc. 

  Class options
  -------------
  The following options control the behavior of this class. The arguments provided below are the default arguments. 

      amr.verbosity       = 10          # Controls verbosity. 
      amr.coarsest_domain = 128 128 128 # Number of cells on coarsest domain
      amr.max_amr_depth   = 0           # Maximum amr depth
      amr.max_sim_depth   = -1          # Maximum simulation depth
      amr.fill_ratio      = 1.0         # Fill ratio for grid generation
      amr.irreg_growth    = 2           # How much to grow irregular tagged cells
      amr.buffer_size     = 2           # Number of cells between grid levels
      amr.blocking_factor = 8           # Default blocking factor (16 in 3D)
      amr.max_box_size    = 16          # Maximum allowed box size
      amr.max_ebis_box    = 32          # Maximum allowed box size
      amr.ref_rat         = 2 2 2 2 2 2 # Refinement ratios
      amr.num_ghost       = 3           # Number of ghost cells. Default is 3
      amr.eb_ghost        = 4           # Set number of of ghost cells for EB stuff
      amr.redist_radius   = 1           # Redistribution radius for hyperbolic conservation laws
      amr.stencil_order   = 1           # Order for interpolation/extrapolation stencils. Must be 1 or 2
      amr.stencil_radius  = 1           # Radius for interpolation/extrapolation stencils. Must be 1 or 2
      amr.stencil_type    = linear      # Default stencil type. Valid options are 'linear', 'taylor', 'lsq'
      amr.ghost_interp    = pwl         # Ghost cell interpolation type. Valid options are 'pwl' or 'quad'
      amr.load_balance    = volume      # Load balancing algorithm. Valid options are 'volume' or 'elliptic'
      amr.ebcf            = false       # Tell amr to ignore EBCF-related code. 
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
    @brief Compute cell-centered gradient. 
    @details User must make sure ghost cells are updated when entering this routine. On output, the result is not averaged down, 
    nor are ghost cells updated. 
  */
  virtual void compute_gradient(EBAMRCellData& a_gradient, const EBAMRCellData& a_phi);

    /*!
    @brief Compute cell-centered gradient. 
    @details User must make sure ghost cells are updated when entering this routine. On output, the result is not averaged down, 
    nor are ghost cells updated. 
  */
  virtual void compute_gradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi);

  /*!
    @brief Create pointers
  */
  virtual void allocate_ptr(EBAMRCellData& a_data);

  /*!
    @brief Allocate a data holder. 
    @note This routines _allocates_, so if a_data is not empty you're leaking memory
  */
  virtual void allocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Allocate a data holder. 
    @note This routines _allocates_, so if a_data is not empty you're leaking memory
  */
  virtual void allocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Allocate a data holder. 
    @note This routines _allocates_, so if a_data is not empty you're leaking memory
  */
  virtual void allocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost = -1);

  /*!
    @brief Allocate a data holder. 
    @note This routines _allocates_, so if a_data is not empty you're leaking memory
  */
  virtual void allocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_dir, const int a_ncomp, const int a_ghost = -1);

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
  virtual void build_grids(Vector<IntVectSet>& a_tags, const int a_hardcap = -1);
  
  /*!
    @brief Define EBLevelGrid
  */
  virtual void define_eblevelgrid();

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
  virtual void set_balance(const load_balance::which_balance a_load);

  /*!
    @brief Set ghost cell interpolation type. Either piecewise linear or quadratic. 
  */
  virtual void set_ghost_interpolation(const ghost_interpolation::which_type a_interp);

  /*!
    @brief Build domains
  */
  virtual void build_domains();

  /*!
    @brief Set the physical domain
  */
  virtual void set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom);

  /*!
    @brief Set verbosity
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Set coarsest level box; this is the coarsest domain (well, except for lower levels GMG)
  */
  virtual void set_coarsest_num_cells(const IntVect a_num_cells);

  /*!
    @brief Maximum AMR depth
  */
  virtual void set_max_amr_depth(const int a_max_amr_depth);

  /*!
    @brief Set maximum simulation depth. 
  */
  virtual void set_max_simulation_depth(const int a_max_sim_depth);

  /*!
    @brief Tell amr_mesh to forget about ebcf (user must guarantee that this is the case)
  */
  virtual void set_ebcf(const bool a_ebcf);

  /*!
    @brief Set refinement ratio
  */
  virtual void set_refinement_ratio(const int a_refinement_ratio);

  /*!
    @brief Set refinement ratios
  */
  virtual void set_refinement_ratios(const Vector<int> a_ref_ratios);

  /*!
    @brief Set buffer size
  */
  virtual void set_buffer_size(const int a_buffer_size);

  /*!
    @brief Set irregular growth
  */
  virtual void set_irreg_growth(const int a_irreg_growth);

  /*!
    @brief Set fill ratio
  */
  virtual void set_fill_ratio(const Real a_fill_ratio);

  /*!
    @brief Set the finest level
  */
  virtual void set_finest_level(const int a_finest_level);

  /*!
    @brief Set grids on level
  */
  virtual void set_grids(Vector<Vector<Box> >& boxes);

  /*!
    @brief Set the maximum box size
  */
  virtual void set_max_box_size(const int a_max_box_size);

  /*!
    @brief Set the maximum EBIS generation box size
  */
  virtual void set_max_ebis_box_size(const int a_max_ebis_box_size);

  /*!
    @brief Set the blocking factor
  */
  virtual void set_blocking_factor(const int a_blocking_factor);

  /*!
    @brief Set number of ghost cells for eb stuff
  */ 
  virtual void set_eb_ghost(const int a_ebghost);

  /*!
    @brief Set number of ghost cells 
  */ 
  virtual void set_num_ghost(const int a_num_ghost);

  /*!
    @brief Set redistribution radius
  */
  virtual void set_redist_rad(const int a_redist_rad);

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
  */
  virtual void regrid(const Vector<IntVectSet>& a_tags, const int a_hardcap);

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
    @brief Get blocking factor
  */
  virtual int get_blocking_factor();

  /*!
    @brief Get maximum box size
  */
  virtual int get_max_box_size();

  /*!
    @brief Get maximum box size for EBIS generation
  */
  virtual int get_max_ebis_box_size();

  /*!
    @brief Get number of ghost cells
  */
  virtual int get_num_ghost();

  /*!
    @brief Get redistribution radius
  */
  virtual int get_redist_rad();
  
  /*!
    @brief Get finest level
  */
  virtual Real get_finest_dx();

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
    @brief Get the grids
  */
  virtual Vector<DisjointBoxLayout>& get_grids();

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
    @brief Get the coarsening utility
  */
  virtual Vector<RefCountedPtr<ebcoarseaverage> >& get_coarave(phase::which_phase a_phase);

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

  
protected:

  /*!
    @brief mfis
  */
  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief Physical domain (needed for computing dx)
  */
  RefCountedPtr<physical_domain> m_physdom;

  /*!
    @brief Load balance type
  */
  load_balance::which_balance m_which_balance;

  /*!
    @brief Stencil type
  */
  stencil_type::which_type m_stencil_type;

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
    @brief EBLevel grids. One for each phase.
  */
  Vector<Vector<RefCountedPtr<EBLevelGrid> > > m_eblg;

  /*!
    @brief MFLevelGrid objects. Cheap objects, and useful for getting interface region between EBIS'
  */
  Vector<RefCountedPtr<MFLevelGrid> > m_mflg;

    /*!
    @brief EBIS layout. One for each phase.
  */
  Vector<Vector<EBISLayout> > m_ebisl;

  /*!
    @brief ebcoarseaverage. One for each phase. 
  */
  Vector<Vector<RefCountedPtr<ebcoarseaverage> > > m_coarave;

  /*!
    @brief QuadCFInterp. One for each phase
  */
  Vector<Vector<RefCountedPtr<nwoebquadcfinterp> > > m_quadcfi;

  /*!
    @brief Old type quadcfi. Needed for EBConductivtyOps
  */
  Vector<Vector<RefCountedPtr<EBQuadCFInterp> > > m_old_quadcfi;

  /*!
    @brief Patcher for piecewise linear interpolation (one ghost cell layer only)
  */
  Vector<Vector<RefCountedPtr<AggEBPWLFillPatch> > > m_pwl_fillpatch;

  /*!
    @brief Piecewise linear interpolation 
  */
  Vector<Vector<RefCountedPtr<EBPWLFineInterp> > > m_pwl_interp;
  
  /*!
    @brief Flux register
  */
  Vector<Vector<RefCountedPtr<EBFluxRegister> > > m_flux_reg;

  /*!
    @brief Level redistribution
  */
  Vector<Vector<RefCountedPtr<EBLevelRedist> > > m_level_redist;

  /*!
    @brief Coarse to fine redistribution
  */
  Vector<Vector<RefCountedPtr<EBCoarToFineRedist> > > m_coar_to_fine_redist;

  /*!
    @brief Coarse to coarse redistribution
  */
  Vector<Vector<RefCountedPtr<EBCoarToCoarRedist> > > m_coar_to_coar_redist;
  
  /*!
    @brief Fine to coarse redistribution
  */
  Vector<Vector<RefCountedPtr<EBFineToCoarRedist> > > m_fine_to_coar_redist;

  /*!
    @brief Interpolation stencil for interpolating to cell centroid
  */
  Vector<RefCountedPtr<irreg_amr_stencil<centroid_interp> > > m_centroid_interp;

  /*!
    @brief Interpolation stencil for interpolating to boundary centroid
  */
  Vector<RefCountedPtr<irreg_amr_stencil<eb_centroid_interp> > > m_eb_centroid_interp;

  /*!
    @brief Load balance grids
  */
  virtual void loadbalance(Vector<Vector<int> >& a_procs, Vector<Vector<Box> >& a_boxes);

  /*!
    Define MFLevelGrid
  */
  virtual void define_mflevelgrid();

  /*!
    @brief Define ebcoarseaverage
  */
  virtual void define_eb_coar_ave();

  /*!
    @brief Define nwoebquadcfinterp
  */
  virtual void define_eb_quad_cfi();

  /*!
    @brief Define operators for piecewise linear interpolation of ghost cells
  */
  virtual void define_fillpatch();

  /*!
    @brief Define interpolator for regridding
  */
  virtual void define_ebpwl_interp();

  /*!
    @brief Define redistribution operators
  */
  virtual void define_redist_oper();


  /*!
    @brief Define flux registers
  */
  virtual void define_flux_reg();

  /*!
    @brief Define irregular stencils
  */
  virtual void define_irreg_sten();
};

#include "amr_meshI.H"

#endif
