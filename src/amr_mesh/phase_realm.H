/*!
  @file   phase_realm.H
  @brief  Declaration of the phase_realm class
  @author Robert Marskar
  @date   July 2020
*/

#ifndef _PHASE_REALM_H_
#define _PHASE_REALM_H_

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
#include "EBFasterFR.H"

#include <DisjointBoxLayout.H>
#include <AggEBPWLFillPatch.H>
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

// These are operator that can be defined.
static const std::string s_eb_gradient     = "eb_gradient";
static const std::string s_eb_irreg_interp = "eb_irreg_interp";
static const std::string s_eb_coar_ave     = "eb_coar_ave";
static const std::string s_eb_quad_cfi     = "eb_quad_cfi";
static const std::string s_eb_fill_patch   = "eb_fill_patch";
static const std::string s_eb_pwl_interp   = "eb_pwl_interp";
static const std::string s_eb_mg_interp    = "eb_mg_interp";
static const std::string s_eb_flux_reg     = "eb_flux_reg";
static const std::string s_eb_redist       = "eb_redist";
static const std::string s_eb_copier       = "eb_copier";
static const std::string s_eb_ghostcloud   = "eb_ghostcloud";
static const std::string s_eb_noncons_div  = "eb_non_cons_div";
static const std::string s_eb_fast_fr      = "eb_fast_fr";

/*!
  @brief Class that holds important things for an AMR phase_realm. I.e. grids, operators, etc. 
*/
class phase_realm {
public:

  // ctor, dtor
  phase_realm();
  ~phase_realm();

  // Define function
  void define(const Vector<DisjointBoxLayout>& a_grids,
	      const Vector<ProblemDomain>& a_domains,
	      const Vector<int>& a_ref_rat,
	      const Vector<Real>& a_dx,
	      const int a_finest_level,
	      const int a_ebghost,
	      const int a_num_ghost,
	      const int a_redist_rad,
	      const stencil_type a_centroid_stencil,
	      const stencil_type a_eb_stencil,
	      const bool a_ebcf,
	      const RefCountedPtr<EBIndexSpace>& a_ebis);

  // Base and operator regrid
  void set_grids(const Vector<DisjointBoxLayout>& a_grids, const int a_finest_level);
  void regrid_base(const int a_lmin);
  void regrid_operators(const int a_lmin, const int a_lmax, const int a_regsize);

  // Register and query operator
  void register_operator(const std::string a_operator);
  bool query_operator(const std::string a_operator);

  // Get functions
  const RefCountedPtr<EBIndexSpace>& get_ebis();
  Vector<int>& get_ref_rat();
  Vector<Real>& get_dx();
  Vector<DisjointBoxLayout>& get_grids();
  Vector<ProblemDomain>& get_domains();
  Vector<EBISLayout>& get_ebisl();
  Vector<RefCountedPtr<EBLevelGrid> >& get_eblg();
  Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& get_neighbors();
  Vector<RefCountedPtr<LayoutData<VoFIterator> > >& get_vofit();
  irreg_amr_stencil<centroid_interp>& get_centroid_interp_stencils();
  irreg_amr_stencil<eb_centroid_interp>& get_eb_centroid_interp_stencils();
  irreg_amr_stencil<noncons_div>& get_noncons_div_stencils();
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > >& get_gradsten();

  Vector<RefCountedPtr<ebcoarseaverage> >& get_coarave();
  Vector<RefCountedPtr<EBGhostCloud> >& get_ghostcloud();
  Vector<RefCountedPtr<nwoebquadcfinterp> >& get_quadcfi();
  Vector<RefCountedPtr<EBQuadCFInterp> >& get_old_quadcfi();
  Vector<RefCountedPtr<AggEBPWLFillPatch> >& get_fillpatch();
  Vector<RefCountedPtr<EBPWLFineInterp> >& get_eb_pwl_interp();
  Vector<RefCountedPtr<EBMGInterp> >& get_eb_mg_interp();
  Vector<RefCountedPtr<EBFluxRegister> >&  get_flux_reg();
  Vector<RefCountedPtr<EBFasterFR> >&  get_eb_fast_fr();
  Vector<RefCountedPtr<EBLevelRedist> >&  get_level_redist();
  Vector<RefCountedPtr<EBCoarToFineRedist> >&  get_coar_to_fine_redist();
  Vector<RefCountedPtr<EBCoarToCoarRedist> >&  get_coar_to_coar_redist();
  Vector<RefCountedPtr<EBFineToCoarRedist> >&  get_fine_to_coar_redist();
  Vector<RefCountedPtr<Copier> >& get_copier();
  Vector<RefCountedPtr<Copier> >& get_reverse_copier();
  
protected:

  bool m_empty;
  bool m_defined;
  bool m_ebcf;
  int m_finest_level;
  int m_mg_coarsen;
  int m_verbosity;
  int m_ebghost;
  int m_num_ghost;
  int m_redist_rad;
  stencil_type m_centroid_stencil;
  stencil_type m_eb_stencil;
  Vector<Real> m_dx;
  Vector<int> m_ref_ratios;

  RefCountedPtr<EBIndexSpace> m_ebis;

  std::map<std::string, bool> m_operator_map; // Map that holds operators

  // AMR grids
  Vector<DisjointBoxLayout> m_grids;

  // Domains
  Vector<ProblemDomain> m_domains;

  // EB stuff
  Vector<EBISLayout> m_ebisl;
  Vector<RefCountedPtr<EBLevelGrid> > m_eblg;

  // Neighbors
  Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > > m_neighbors;

  // Iterator for vofs
  Vector<RefCountedPtr<LayoutData<VoFIterator> > > m_vofiter;

  // Copiers
  Vector<RefCountedPtr<Copier> > m_copier;
  Vector<RefCountedPtr<Copier> > m_reverse_copier;

  // Various operators
  Vector<RefCountedPtr<ebcoarseaverage> > m_coarave;
  Vector<RefCountedPtr<nwoebquadcfinterp> > m_quadcfi;
  Vector<RefCountedPtr<EBGhostCloud> > m_ghostclouds;
  Vector<RefCountedPtr<EBQuadCFInterp> > m_old_quadcfi;
  Vector<RefCountedPtr<AggEBPWLFillPatch> > m_pwl_fillpatch;
  Vector<RefCountedPtr<EBPWLFineInterp> > m_pwl_interp;
  Vector<RefCountedPtr<EBMGInterp> > m_ebmg_interp;
  Vector<RefCountedPtr<EBFluxRegister> > m_flux_reg;
  Vector<RefCountedPtr<EBFasterFR> > m_eb_fast_fr;
  Vector<RefCountedPtr<EBLevelRedist> > m_level_redist;
  Vector<RefCountedPtr<EBCoarToFineRedist> > m_coar_to_fine_redist;
  Vector<RefCountedPtr<EBCoarToCoarRedist> > m_coar_to_coar_redist;
  Vector<RefCountedPtr<EBFineToCoarRedist> > m_fine_to_coar_redist;

  // Various stencils for irregular cells
  RefCountedPtr<irreg_amr_stencil<centroid_interp> > m_centroid_interp;
  RefCountedPtr<irreg_amr_stencil<eb_centroid_interp> > m_eb_centroid_interp;
  RefCountedPtr<irreg_amr_stencil<noncons_div> > m_noncons_div;
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > > m_gradsten;

  // Base defines
  void define_eblevelgrid(const int a_lmin);
  void define_vofiter(const int a_lmin);
  void define_neighbors(const int a_lmin);
  
  // Operator defines
  void define_eb_coar_ave(const int a_lmin);                       // Conservative coarsening
  void define_eb_quad_cfi(const int a_lmin);                       // Quadratic ghost cell interpolation
  void define_fillpatch(const int a_lmin);                         // Define operator for ghost cell interpolation
  void define_ebpwl_interp(const int a_lmin);                      // Define interpolator for piecewise interpolation
  void define_ebmg_interp(const int a_lmin);                       // Define interpolator used for e.g. multigrid
  void define_flux_reg(const int a_lmin, const int a_regsize);     // Define flux register
  void define_eb_fast_fr(const int a_lmin, const int a_regsize);   // Define EBFastFR
  void define_redist_oper(const int a_lmin, const int a_regsize);  // Define redistribution (phase::gas only)
  void define_gradsten(const int a_lmin);                          // Make stencils for computing gradients
  void define_copier(const int a_lmin);                            // Make stencils for copier
  void define_ghostcloud(const int a_lmin);                        // Make stencils for ghost clouds with particle depositions
  void define_irreg_sten();                                        // Make stencils for doing interpolation to centroids
  void define_noncons_sten();                                      // Make stencils for nonconservative averaging
};

#endif
