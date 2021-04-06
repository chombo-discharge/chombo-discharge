/*!
  @file   realm.H
  @brief  Multi-fluid realms
  @author Robert Marskar
  @date   July 2020
*/

#ifndef _REALM_H_
#define _REALM_H_

#include "phase_realm.H"
#include "mfis.H"
#include "MFLevelGrid.H"

static const std::string s_particle_halo         = "particle_halo";
static const std::string s_particle_valid_region = "pvr";

using AMRMask = Vector<RefCountedPtr<LevelData<BaseFab<bool> > > >;

/*!
  @brief Class for holding multifluid realms, which are the realms that we actually use. 
*/
class realm {
public:

  // Identifier for default realm
  const static std::string primal;

  // ctor, dtor
  realm();
  ~realm();

  // Define function
  void define(const Vector<DisjointBoxLayout>& a_grids,
	      const Vector<ProblemDomain>& a_domains,
	      const Vector<int>& a_ref_rat,
	      const Vector<Real>& a_dx,
	      const RealVect a_prob_lo,
	      const int a_finest_level,
	      const int a_ebghost,
	      const int a_num_ghost,
	      const int a_lsf_ghost,
	      const int a_redist_rad,
	      const bool a_ebcf,
	      const stencil_type a_centroid_stencil,
	      const stencil_type a_eb_stencil,
	      const std::map<phase::which_phase, RefCountedPtr<BaseIF> > a_baseif,
	      const RefCountedPtr<mfis>& a_mfis);

  // realm copy function
  template <class T>
  static void copy(Vector<RefCountedPtr<LevelData<T> > >& a_dst, const Vector<RefCountedPtr<LevelData<T> > >& a_src);

  // Base and operator regrid
  void set_grids(const Vector<DisjointBoxLayout>& a_grids, const int a_finest_level);
  void regrid_base(const int a_lmin);
  void regrid_operators(const int a_lmin, const int a_lmax, const int a_regsize);

  // Register and query operator
  void register_operator(const std::string a_operator, const phase::which_phase a_phase);
  bool query_operator(const std::string a_operator, const phase::which_phase a_phase);

  // Register and query masks
  void register_mask(const std::string a_mask, const int a_buffer);
  bool query_mask(const std::string a_mask, const int a_buffer) const;

  // Get functions for neutral stuff
  Vector<int>& get_ref_rat();
  Vector<Real>& get_dx();
  Vector<DisjointBoxLayout>& get_grids();
  Vector<ProblemDomain>& get_domains();
  Vector<RefCountedPtr<MFLevelGrid> >& get_mflg();
  
  const RefCountedPtr<EBIndexSpace>& get_ebis(const phase::which_phase a_phase);
  Vector<EBISLayout>& get_ebisl(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBLevelGrid> >& get_eblg(const phase::which_phase a_phase);
  Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& get_neighbors(const phase::which_phase a_phase);
  Vector<RefCountedPtr<LayoutData<VoFIterator> > >& get_vofit(const phase::which_phase a_phase);
  irreg_amr_stencil<centroid_interp>& get_centroid_interp_stencils(const phase::which_phase a_phase);
  irreg_amr_stencil<eb_centroid_interp>& get_eb_centroid_interp_stencils(const phase::which_phase a_phase);
  irreg_amr_stencil<noncons_div>& get_noncons_div_stencils(const phase::which_phase a_phase);
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > >& get_gradsten(const phase::which_phase a_phase);

  // Get functions for stuff that are specific to phase, but is not necessarily const
  Vector<RefCountedPtr<ebcoarseaverage> >& get_coarave(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBGhostCloud> >& get_ghostcloud(const phase::which_phase a_phase);
  Vector<RefCountedPtr<nwoebquadcfinterp> >& get_quadcfi(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBQuadCFInterp> >& get_old_quadcfi(const phase::which_phase a_phase);
  Vector<RefCountedPtr<AggEBPWLFillPatch> >& get_fillpatch(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBPWLFineInterp> >& get_eb_pwl_interp(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBMGInterp> >& get_eb_mg_interp(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBFluxRegister> >&  get_flux_reg(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBLevelRedist> >&  get_level_redist(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBCoarToFineRedist> >&  get_coar_to_fine_redist(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBCoarToCoarRedist> >&  get_coar_to_coar_redist(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBFineToCoarRedist> >&  get_fine_to_coar_redist(const phase::which_phase a_phase);
  Vector<RefCountedPtr<Copier> >& get_copier(const phase::which_phase a_phase);
  Vector<RefCountedPtr<Copier> >& get_reverse_copier(const phase::which_phase a_phase);

  EBAMRFAB& get_levelset(const phase::which_phase a_phase);

  AMRMask& get_mask(const std::string a_mask, const int a_buffer);
  
protected:

  bool m_defined;
  int m_verbosity;
  int m_finest_level;

  RealVect m_prob_lo;

  Vector<Real> m_dx;
  Vector<int> m_ref_ratios;
  Vector<DisjointBoxLayout> m_grids;
  Vector<ProblemDomain> m_domains;
  Vector<RefCountedPtr<MFLevelGrid> > m_mflg;
  RefCountedPtr<mfis> m_mfis;

  // Phase realms and implicit functions
  std::map<phase::which_phase, RefCountedPtr<phase_realm> > m_realms;
  std::map<phase::which_phase, RefCountedPtr<BaseIF> > m_baseif;

  // Map for checking which masks are present. 
  std::map<std::pair<std::string, int>, bool>    m_masks_map;
  std::map<std::pair<std::string, int>, AMRMask> m_masks;

  void define_mflevelgrid(const int a_lmin);

  // Mask definitions. 
  void define_masks(const int a_lmin);
  void define_halo_masks(const int a_lmin);
  void define_halo_mask(LevelData<BaseFab<bool> >& a_coarMask,
			const ProblemDomain&       a_domainCoar,
			const ProblemDomain&       a_domainFine,
			const DisjointBoxLayout&   a_gridsCoar,
			const DisjointBoxLayout&   a_gridsFine,
			const int                  a_buffer,
			const int                  a_ref_rat);

  // Get realm
  phase_realm& get_realm(const phase::which_phase a_phase);


  
};
#endif
