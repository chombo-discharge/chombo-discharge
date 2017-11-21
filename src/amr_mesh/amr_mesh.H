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

#include <DisjointBoxLayout.H>
#include <ProblemDomain.H>
#include <EBCoarseAverage.H>
#include <nwoebquadcfinterp.H>
#include <EBFastFR.H>
#include <EBLevelRedist.H>
#include <EBCoarToFineRedist.H>
#include <EBCoarToCoarRedist.H>
#include <EBFineToCoarRedist.H>

/*!
  @brief Class for handling spatial operations for streamer equations.
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
    @brief Average down
  */
  virtual void average_down(EBAMRCellData& a_data, Phase::WhichPhase a_phase);

  /*!
    @brief Average down
  */
  virtual void average_down(EBAMRFluxData& a_data, Phase::WhichPhase a_phase);

  /*!
    @brief Average down
  */
  virtual void average_down(EBAMRIVData& a_data, Phase::WhichPhase a_phase);

  /*!
    @brief Interpolate ghost vectors
  */
  virtual void interp_ghost(EBAMRCellData& a_data, Phase::WhichPhase a_phase);
  
  /*!
    @brief Set mfis
  */
  virtual void set_mfis(const RefCountedPtr<mfis>& a_mfis);

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
    @brief Tell amr_mesh to forget about ebcf (user must guarantee that this is the case)
  */
  virtual void set_ebcf(const bool a_ebcf);

  /*!
    @brief Set refinement ratio
  */
  virtual void set_refinement_ratio(const int a_refinement_ratio);

  /*!
    @brief Set buffer size
  */
  virtual void set_buffer_size(const int a_buffer_size);

  /*!
    @brief Set fill ratio
  */
  virtual void set_fill_ratio(const Real a_fill_ratio);

  /*!
    @brief Set the maximum box size
  */
  virtual void set_max_box_size(const int a_max_box_size);

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
    @brief Regrid function. Generate new grids. 
  */
  virtual void regrid(const Vector<IntVectSet>& a_tags);

  /*!
    @brief Do a sanity check to make sure everything is set up correctly
  */
  virtual void sanity_check();

    /*!
    @brief Get refinement ratios
  */
  virtual int get_finest_level();

  /*!
    @brief Get maximum amr depth
  */
  virtual int get_max_amr_depth();

  /*!
    @brief Get blocking factor
  */
  virtual int get_blocking_factor();

  /*!
    @brief Get maximum box size
  */
  virtual int get_max_box_size();

  /*!
    @brief Get refinement ratio
  */
  virtual int get_refinement_ratio();

  /*!
    @brief Get finest domain
  */
  virtual ProblemDomain get_finest_domain();

  /*!
    @brief Get finest level
  */
  virtual Real get_finest_dx();

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
  virtual Vector<RefCountedPtr<EBLevelGrid> >& get_eblg(Phase::WhichPhase a_phase);

    /*!
    @brief Get the EBISLayouts
  */
  virtual Vector<EBISLayout>& get_ebisl(Phase::WhichPhase a_phase);

  /*!
    @brief Get the coarsening utility
  */
  virtual Vector<RefCountedPtr<EBCoarseAverage> >& get_coarave(Phase::WhichPhase a_phase);

  /*!
    @brief Get the ghost cell interpolation utility
  */
  virtual Vector<RefCountedPtr<nwoebquadcfinterp> >& get_quadcfi(Phase::WhichPhase);

  /*!
    @brief Get flux register
  */
  virtual Vector<RefCountedPtr<EBFastFR> >&  get_flux_reg(Phase::WhichPhase a_phase);

  /*!
    @brief Get redistribution
  */
  virtual Vector<RefCountedPtr<EBLevelRedist> >&  get_level_redist(Phase::WhichPhase a_phase);

  /*!
    @brief Get coar-to-fine redistribution
  */
  virtual Vector<RefCountedPtr<EBCoarToFineRedist> >&  get_coar_to_fine_redist(Phase::WhichPhase a_phase);

  /*!
    @brief Get coar-to-coar redistribution
  */
  virtual Vector<RefCountedPtr<EBCoarToCoarRedist> >&  get_coar_to_coar_redist(Phase::WhichPhase a_phase);

    /*!
    @brief Get coar-to-coar redistribution
  */
  virtual Vector<RefCountedPtr<EBFineToCoarRedist> >&  get_fine_to_coar_redist(Phase::WhichPhase a_phase);

  
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
    @brief Max box size
  */
  int m_max_box_size;

  /*!
    @brief Set buffer size
  */
  int m_buffer_size;

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
    @brief Grids have EB x CF crossing
  */
  bool m_ebcf;

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
    @brief EBIS layout. One for each phase.
  */
  Vector<Vector<EBISLayout> > m_ebisl;

  /*!
    @brief EBCoarseAverage. One for each phase. 
  */
  Vector<Vector<RefCountedPtr<EBCoarseAverage> > > m_coarave;

  /*!
    @brief QuadCFInterp. One for each phase
  */
  Vector<Vector<RefCountedPtr<nwoebquadcfinterp> > > m_quadcfi;
  
  /*!
    @brief Flux register
  */
  Vector<Vector<RefCountedPtr<EBFastFR> > > m_flux_reg;

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
    @brief Build new grids. This also load balances. 
  */
  virtual void build_grids(Vector<IntVectSet>& a_tags);

  /*!
    @brief Load balance grids
  */
  virtual void load_balance(Vector<int>& a_proc_assign, Vector<Box>& a_boxes, const int a_lvl);

  /*!
    @brief Define EBLevelGrid
  */
  virtual void define_eblevelgrid();

  /*!
    @brief Define EBCoarseAverage
  */
  virtual void define_eb_coar_ave();

  /*!
    @brief Define nwoebquadcfinterp
  */
  virtual void define_eb_quad_cfi();

  /*!
    @brief Define redistribution operators
  */
  virtual void define_redist_oper();

  /*!
    @brief Define flux registers
  */
  virtual void define_flux_reg();
};

#endif
