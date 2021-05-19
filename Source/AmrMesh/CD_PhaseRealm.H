/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_PhaseRealm.H
  @brief  Declaration of the PhaseRealm class
  @author Robert Marskar
*/

#ifndef CD_PhaseRealm_H
#define CD_PhaseRealm_H

// Chombo includes
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

// Our includes
#include <CD_ComputationalGeometry.H>
#include <mfis.H>
#include <CD_IrregAmrStencil.H>
#include <CD_LoadBalancing.H>
#include <CD_NwoEbQuadCfInterp.H>
#include <MFLevelGrid.H>
#include <ebcoarseaverage.H>
#include <DomainFluxIFFAB.H>
#include <DomainFluxIFFABFactory.H>
#include <CD_ParticleContainer.H>
#include <CD_EbGhostCloud.H>
#include <CD_NonConservativeDivergenceStencil.H>
#include <CD_EbCentroidInterpolationStencil.H>
#include <CD_CentroidInterpolationStencil.H>
#include <CD_EbFastFluxRegister.H>
#include <CD_NamespaceHeader.H>

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
static const std::string s_noncons_div  = "eb_non_cons_div";
static const std::string s_levelset        = "levelset";

/*!
  @brief Class that holds important things for an AMR PhaseRealm. I.e. grids, operators, etc. 
*/
class PhaseRealm {
public:

  // ctor, dtor
  PhaseRealm();
  ~PhaseRealm();

  // Define function
  void define(const Vector<DisjointBoxLayout>& a_grids,
	      const Vector<ProblemDomain>& a_domains,
	      const Vector<int>& a_ref_rat,
	      const Vector<Real>& a_dx,
	      const RealVect a_probLo,
	      const int a_finestLevel,
	      const int a_ebghost,
	      const int a_num_ghost,
	      const int a_lsf_ghost,
	      const int a_redist_rad,
	      const IrregStencil::StencilType a_centroid_stencil,
	      const IrregStencil::StencilType a_eb_stencil,
	      const bool a_ebcf,
	      const RefCountedPtr<BaseIF>&       a_baseif,
	      const RefCountedPtr<EBIndexSpace>& a_ebis);

  // Base and operator regrid
  void setGrids(const Vector<DisjointBoxLayout>& a_grids, const int a_finestLevel);
  void regridBase(const int a_lmin);
  void regridOperators(const int a_lmin, const int a_lmax, const int a_regsize);

  // Register and query operator
  void registerOperator(const std::string a_operator);
  bool queryOperator(const std::string a_operator);

  // Get functions
  const RefCountedPtr<EBIndexSpace>& getEBIndexSpace();
  Vector<int>& getRefinementRatios();
  Vector<Real>& getDx();
  Vector<DisjointBoxLayout>& getGrids();
  const Vector<DisjointBoxLayout>& getGrids() const;
  Vector<ProblemDomain>& getDomains();
  Vector<EBISLayout>& getEBISLayout();
  Vector<RefCountedPtr<EBLevelGrid> >& getEBLevelGrid();
  Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& getNeighbors();
  Vector<RefCountedPtr<LayoutData<VoFIterator> > >& getVofIterator();
  IrregAmrStencil<CentroidInterpolationStencil>& getCentroidInterpolationStencils();
  IrregAmrStencil<EbCentroidInterpolationStencil>& getEbCentroidInterpolationStencilStencils();
  IrregAmrStencil<NonConservativeDivergenceStencil>& getNonConservativeDivergenceStencils();
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > >& getGradientStencils();

  Vector<RefCountedPtr<ebcoarseaverage> >& getCoarseAverage();
  Vector<RefCountedPtr<EbGhostCloud> >& getGhostCloud();
  Vector<RefCountedPtr<NwoEbQuadCfInterp> >& getNWOEBQuadCFInterp();
  Vector<RefCountedPtr<EBQuadCFInterp> >& getEBQuadCFInterp();
  Vector<RefCountedPtr<AggEBPWLFillPatch> >& getFillPatch();
  Vector<RefCountedPtr<EBPWLFineInterp> >& getPwlInterpolator();
  Vector<RefCountedPtr<EBMGInterp> >& getEBMGInterp();
  Vector<RefCountedPtr<EBFluxRegister> >&  getFluxRegister();
  Vector<RefCountedPtr<EBLevelRedist> >&  getLevelRedist();
  Vector<RefCountedPtr<EBCoarToFineRedist> >&  getCoarToFineRedist();
  Vector<RefCountedPtr<EBCoarToCoarRedist> >&  getCoarToCoarRedist();
  Vector<RefCountedPtr<EBFineToCoarRedist> >&  getFineToCoarRedist();
  Vector<RefCountedPtr<Copier> >& getCopier();
  Vector<RefCountedPtr<Copier> >& getReverseCopier();
  EBAMRFAB& getLevelset();
  
protected:

  bool m_empty;
  bool m_defined;
  bool m_hasEbCf;
  int m_finestLevel;
  int m_mg_coarsen;
  int m_verbosity;
  int m_numGhostCells;
  int m_numEbGhostsCells;
  int m_numLsfGhostCells;
  int m_redistributionRadius;
  RealVect m_probLo;
  IrregStencil::StencilType m_centroidStencilType;
  IrregStencil::StencilType m_ebCentroidStencilType;
  Vector<Real> m_dx;
  Vector<int> m_refinementRatios;

  RefCountedPtr<EBIndexSpace> m_ebis;
  RefCountedPtr<BaseIF> m_baseif;

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

  // Iterator for vofs. 
  Vector<RefCountedPtr<LayoutData<VoFIterator> > > m_vofiter;

  // Copiers
  Vector<RefCountedPtr<Copier> > m_copier;
  Vector<RefCountedPtr<Copier> > m_reverse_copier;

  // Levelset
  EBAMRFAB m_levelset;

  // Various operators
  Vector<RefCountedPtr<ebcoarseaverage> > m_coarave;
  Vector<RefCountedPtr<NwoEbQuadCfInterp> > m_quadcfi;
  Vector<RefCountedPtr<EbGhostCloud> > m_ghostclouds;
  Vector<RefCountedPtr<EBQuadCFInterp> > m_old_quadcfi;
  Vector<RefCountedPtr<AggEBPWLFillPatch> > m_pwl_fillpatch;
  Vector<RefCountedPtr<EBPWLFineInterp> > m_pwl_interp;
  Vector<RefCountedPtr<EBMGInterp> > m_ebmg_interp;
  Vector<RefCountedPtr<EBFluxRegister> > m_flux_reg;
  Vector<RefCountedPtr<EBLevelRedist> > m_level_redist;
  Vector<RefCountedPtr<EBCoarToFineRedist> > m_coar_to_fine_redist;
  Vector<RefCountedPtr<EBCoarToCoarRedist> > m_coar_to_coar_redist;
  Vector<RefCountedPtr<EBFineToCoarRedist> > m_fine_to_coar_redist;

  // Various stencils for irregular cells
  RefCountedPtr<IrregAmrStencil<CentroidInterpolationStencil> > m_CentroidInterpolationStencil;
  RefCountedPtr<IrregAmrStencil<EbCentroidInterpolationStencil> > m_EbCentroidInterpolationStencil;
  RefCountedPtr<IrregAmrStencil<NonConservativeDivergenceStencil> > m_NonConservativeDivergenceStencil;
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > > m_gradsten;

  // Base defines
  void define_eblevelgrid(const int a_lmin);
  void defineVofIterator(const int a_lmin);
  void define_neighbors(const int a_lmin);
  
  // Operator defines
  void define_eb_coar_ave(const int a_lmin);                       // Conservative coarsening
  void define_eb_quad_cfi(const int a_lmin);                       // Quadratic ghost cell interpolation
  void define_fillpatch(const int a_lmin);                         // Define operator for ghost cell interpolation
  void define_ebpwl_interp(const int a_lmin);                      // Define interpolator for piecewise interpolation
  void define_ebmg_interp(const int a_lmin);                       // Define interpolator used for e.g. multigrid
  void define_flux_reg(const int a_lmin, const int a_regsize);     // Define flux register
  void define_redist_oper(const int a_lmin, const int a_regsize);  // Define redistribution (phase::gas only)
  void define_gradsten(const int a_lmin);                          // Make stencils for computing gradients
  void define_copier(const int a_lmin);                            // Make stencils for copier
  void define_ghostcloud(const int a_lmin);                        // Make stencils for ghost clouds with particle depositions
  void define_irreg_sten();                                        // Make stencils for doing interpolation to centroids
  void define_noncons_sten();                                      // Make stencils for nonconservative averaging
  void define_levelset(const int a_lmin, const int a_numGhost);    // Define levelset
};

#include <CD_NamespaceFooter.H>

#endif
