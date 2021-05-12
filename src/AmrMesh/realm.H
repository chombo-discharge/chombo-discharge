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

#include "CD_NamespaceHeader.H"

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
	      const RealVect a_probLo,
	      const int a_finest_level,
	      const int a_ebghost,
	      const int a_num_ghost,
	      const int a_lsf_ghost,
	      const int a_redist_rad,
	      const bool a_ebcf,
	      const IrregStencil::StencilType a_centroid_stencil,
	      const IrregStencil::StencilType a_eb_stencil,
	      const std::map<phase::which_phase, RefCountedPtr<BaseIF> > a_baseif,
	      const RefCountedPtr<mfis>& a_mfis);

  // realm copy function
  template <class T>
  static void copy(Vector<RefCountedPtr<LevelData<T> > >& a_dst, const Vector<RefCountedPtr<LevelData<T> > >& a_src);

  // Base and operator regrid
  void setGrids(const Vector<DisjointBoxLayout>& a_grids, const int a_finest_level);
  void regrid_base(const int a_lmin);
  void regridOperators(const int a_lmin, const int a_lmax, const int a_regsize);

  // Register and query operator
  void registerOperator(const std::string a_operator, const phase::which_phase a_phase);
  bool query_operator(const std::string a_operator, const phase::which_phase a_phase);

  // Register and query masks
  void registerMask(const std::string a_mask, const int a_buffer);
  bool query_mask(const std::string a_mask, const int a_buffer) const;

  // Get functions for neutral stuff
  Vector<int>& getRefinementRatios();
  Vector<Real>& getDx();
  const Vector<Real>& getDx() const;
  Vector<DisjointBoxLayout>& getGrids();
  const Vector<DisjointBoxLayout>& getGrids() const;
  Vector<ProblemDomain>& getDomains();
  const Vector<ProblemDomain>& getDomains() const;
  Vector<RefCountedPtr<MFLevelGrid> >& getMFLevelGrid();
  
  const RefCountedPtr<EBIndexSpace>& get_ebis(const phase::which_phase a_phase);
  Vector<EBISLayout>& getEBISLayout(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBLevelGrid> >& getEBLevelGrid(const phase::which_phase a_phase);
  Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& getNeighbors(const phase::which_phase a_phase);
  Vector<RefCountedPtr<LayoutData<VoFIterator> > >& getVofIterator(const phase::which_phase a_phase);
  irreg_amr_stencil<CentroidInterpolationStencil>& getCentroidInterpolationStencils(const phase::which_phase a_phase);
  irreg_amr_stencil<EbCentroidInterpolationStencil>& getEbCentroidInterpolationStencilStencils(const phase::which_phase a_phase);
  irreg_amr_stencil<NonConservativeDivergenceStencil>& getNonConservativeDivergenceStencils(const phase::which_phase a_phase);
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > >& get_gradsten(const phase::which_phase a_phase);

  // Get functions for stuff that are specific to phase, but is not necessarily const
  Vector<RefCountedPtr<ebcoarseaverage> >& getCoarseAverage(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBGhostCloud> >& getGhostCloud(const phase::which_phase a_phase);
  Vector<RefCountedPtr<nwoebquadcfinterp> >& getNWOEBQuadCFInterp(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBQuadCFInterp> >& getEBQuadCFInterp(const phase::which_phase a_phase);
  Vector<RefCountedPtr<AggEBPWLFillPatch> >& getFillPatch(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBPWLFineInterp> >& getPwlInterpolator(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBMGInterp> >& getEBMGInterp(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBFluxRegister> >&  getFluxRegister(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBLevelRedist> >&  getLevelRedist(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBCoarToFineRedist> >&  getCoarToFineRedist(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBCoarToCoarRedist> >&  getCoarToCoarRedist(const phase::which_phase a_phase);
  Vector<RefCountedPtr<EBFineToCoarRedist> >&  getFineToCoarRedist(const phase::which_phase a_phase);
  Vector<RefCountedPtr<Copier> >& getCopier(const phase::which_phase a_phase);
  Vector<RefCountedPtr<Copier> >& getReverseCopier(const phase::which_phase a_phase);

  EBAMRFAB& getLevelset(const phase::which_phase a_phase);

  AMRMask& getMask(const std::string a_mask, const int a_buffer);
  
protected:

  bool m_defined;
  int m_verbosity;
  int m_finestLevel;

  RealVect m_probLo;

  Vector<Real> m_dx;
  Vector<int> m_refinementRatios;
  Vector<DisjointBoxLayout> m_grids;
  Vector<ProblemDomain> m_domains;
  Vector<RefCountedPtr<MFLevelGrid> > m_mflg;
  RefCountedPtr<mfis> m_multifluidIndexSpace;

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
#include "CD_NamespaceFooter.H"
#endif
