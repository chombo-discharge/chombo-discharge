/*!
  @file poisson_multifluid_gmg.cpp
  @brief Implementation of poisson_multifluid_gmg.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "poisson_multifluid_gmg.H"
#include "data_ops.H"
#include "MFQuadCFInterp.H"
#include "MFInterfaceFAB.H"
#include "jump_bc.H"

#include <Stencils.H>
#include <MFCellFAB.H>
#include <LayoutData.H>
#include <MFLevelDataOps.H>

poisson_multifluid_gmg::poisson_multifluid_gmg(){
  m_needs_setup = true;
}

poisson_multifluid_gmg::~poisson_multifluid_gmg(){

}

int poisson_multifluid_gmg::query_ghost() const {
  return 2; // Need this many cells
}

void poisson_multifluid_gmg::solve(){
  if(m_needs_setup){
    this->setup_gmg();
  }

  
}

void poisson_multifluid_gmg::solve(MFAMRCellData& a_state, const MFAMRCellData& a_source){

}

void poisson_multifluid_gmg::set_coefficients(){
  CH_TIME("poisson_multifluid_gmg::set_coefficients");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::set_coefficients" << endl;
  }

  const int ncomps = 1;
  const int ghosts = 1;
  const Real eps0  = m_compgeom->get_eps0();
  
  m_amr->allocate(m_aco,       ncomps, ghosts);
  m_amr->allocate(m_bco,       ncomps, ghosts);
  m_amr->allocate(m_bco_irreg, ncomps, ghosts);

  data_ops::set_value(m_aco,       0.0);   // Always zero for poisson equation
  data_ops::set_value(m_bco,       eps0); // Set equal to something large to detect bugs
  data_ops::set_value(m_bco_irreg, eps0); // Set equal to something large to detect bugs

  this->set_permittivities(m_compgeom->get_dielectrics());
}

void poisson_multifluid_gmg::set_permittivities(const Vector<dielectric>& a_dielectrics){
  CH_TIME("poisson_multifluid_gmg::set_permittivities");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::set_permittivities" << endl;
  }

  if(a_dielectrics.size() > 0){
    const RealVect origin  = m_physdom->get_prob_lo();
    const Vector<Real> dx  = m_amr->get_dx();
    const int finest_level = m_amr->get_finest_level();

    for (int lvl = 0; lvl <= finest_level; lvl++){

      LevelData<EBFluxFAB> bco;
      LevelData<BaseIVFAB<Real> > bco_irreg;

      mfalias::aliasMF(bco,       phase::solid, *m_bco[lvl]);
      mfalias::aliasMF(bco_irreg, phase::solid, *m_bco_irreg[lvl]);

      for (DataIterator dit = bco.dataIterator(); dit.ok(); ++dit){
	EBFluxFAB& perm          = bco[dit()];
	BaseIVFAB<Real>& perm_eb = bco_irreg[dit()];

	this->set_face_perm(perm,  origin, dx[lvl], a_dielectrics);
	this->set_eb_perm(perm_eb, origin, dx[lvl], a_dielectrics);
      }
    }
  }
}

void poisson_multifluid_gmg::set_face_perm(EBFluxFAB&                a_perm,
					   const RealVect&           a_origin,
					   const Real&               a_dx,
					   const Vector<dielectric>& a_dielectrics){
  CH_TIME("poisson_multifluid_gmg::set_face_perm");
  if(m_verbosity > 10){
    pout() << "poisson_multifluid_gmg::set_face_perm" << endl;
  }

  const int comp         = 0;
  const IntVectSet ivs   = IntVectSet(a_perm.getRegion());
  const EBISBox& ebisbox = a_perm.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  FaceStop::WhichFaces stop_crit = FaceStop::SurroundingWithBoundary;
  
  for (int dir = 0; dir < SpaceDim; dir++){
    for (FaceIterator faceit(ivs, ebgraph, dir, stop_crit); faceit.ok(); ++faceit){
      const FaceIndex& face  = faceit();
      const IntVect iv       = face.gridIndex(Side::Lo);
      const RealVect pos     = a_origin + a_dx*iv + 0.5*a_dx*BASISV(dir);
      
      Real dist   = 1.E99;
      int closest = 0;
      for (int i = 0; i < a_dielectrics.size(); i++){
	const RefCountedPtr<BaseIF> func = a_dielectrics[i].get_function();

	const Real cur_dist = func->value(pos);
	
	if(cur_dist <= dist){
	  dist = cur_dist;
	  closest = i;
	}
      }

      a_perm[dir](face, comp) = a_dielectrics[closest].get_permittivity();
    }
  }
}


void poisson_multifluid_gmg::set_eb_perm(BaseIVFAB<Real>&          a_perm,
					 const RealVect&           a_origin,
					 const Real&               a_dx,
					 const Vector<dielectric>& a_dielectrics){
  CH_TIME("poisson_multifluid_gmg::set_eb_perm");
  if(m_verbosity > 10){
    pout() << "poisson_multifluid_gmg::set_eb_perm" << endl;
  }
  
  const int comp         = 0;
  const IntVectSet ivs   = a_perm.getIVS();
  const EBGraph& ebgraph = a_perm.getEBGraph();

  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect pos  = EBArith::getVofLocation(vof, a_dx, a_origin); // This is strictly speaking not on the boundary...
      
    Real dist   = 1.E99;
    int closest = 0;
    for (int i = 0; i < a_dielectrics.size(); i++){
      const RefCountedPtr<BaseIF> func = a_dielectrics[i].get_function();

      const Real cur_dist = func->value(pos);
	
      if(cur_dist <= dist){
	dist = cur_dist;
	closest = i;
      }
    }

    a_perm(vof, comp) = a_dielectrics[closest].get_permittivity();
  }
}

void poisson_multifluid_gmg::setup_gmg(){
  CH_TIME("poisson_multifluid_gmg::setup_gmg");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::setup_gmg" << endl;
  }
  
  this->set_coefficients(); // Set coefficients

}

void poisson_multifluid_gmg::base_tests(){

  this->base_tests();
  const int finest_level = m_amr->get_finest_level();
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids();

  m_state.resize(1 + finest_level);
  m_source.resize(1 + finest_level);

  // Create data holders
  for (int lvl = 0; lvl <= finest_level; lvl++){
    Vector<EBISLayout> ebisl_phases(2);
    Vector<int> components(2);
    ebisl_phases[0] = m_amr->get_ebisl(phase::gas)[lvl];
    ebisl_phases[1] = m_amr->get_ebisl(phase::solid)[lvl];
    components[0] = 1;
    components[1] = 1;
    
    MFCellFactory cellfact(ebisl_phases, components);

    m_state[lvl] = RefCountedPtr<LevelData<MFCellFAB> > (new LevelData<MFCellFAB>(grids[lvl], 1, 3*IntVect::Unit, cellfact));
    m_source[lvl] =RefCountedPtr<LevelData<MFCellFAB> > (new LevelData<MFCellFAB>(grids[lvl], 1, 3*IntVect::Unit, cellfact));

    MFLevelDataOps::setVal(*m_state[lvl], 0.0);
    MFLevelDataOps::setVal(*m_source[lvl], 0.0);
  }

  // Things to pass to the operator factory
  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);
  const Vector<int>& refinement_ratios  = m_amr->get_ref_rat();
  const Vector<ProblemDomain>& domains  = m_amr->get_domains();
  const Vector<Real>& dx                = m_amr->get_dx();
  const RealVect& origin                = m_physdom->get_prob_lo();

  
  RefCountedPtr<BaseDomainBCFactory> domfact = RefCountedPtr<BaseDomainBCFactory> (NULL);
  RefCountedPtr<BaseEBBCFactory>     ebcfact = RefCountedPtr<BaseEBBCFactory> (NULL);

  EBAMRIVData   sigma;
  Real          alpha =  0.0;
  Real          beta  = -1.0;

  const int ncomps = 1;
  const int ghosts = 1;
  m_amr->allocate(sigma, phase::gas, 1, 0);


  // Create two MFBaseIVFABs and copy between them
  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<EBISLayout> ebisl(2);
    Vector<int> comps(2);

    ebisl[phase::gas] = m_amr->get_ebisl(phase::gas)[lvl];
    ebisl[phase::solid] = m_amr->get_ebisl(phase::solid)[lvl];

    comps[phase::gas] = 1;
    comps[phase::solid] = 1;

    MFBaseIVFABFactory fact(ebisl, comps);

    const int ignored = 0;
    LevelData<MFBaseIVFAB>* mfiv1 = new LevelData<MFBaseIVFAB> (grids[lvl], ignored, 2*IntVect::Unit, fact);
    LevelData<MFBaseIVFAB>* mfiv2 = new LevelData<MFBaseIVFAB> (grids[lvl], ignored, 2*IntVect::Unit, fact);

    data_ops::set_value(*mfiv1, 1.0);
    data_ops::set_value(*mfiv2, 2.0);

    mfiv1->copyTo(*mfiv2);

    Vector<EBLevelGrid> eblg;
    eblg.push_back(*(m_amr->get_eblg(phase::gas)[lvl]));
    eblg.push_back(*(m_amr->get_eblg(phase::solid)[lvl]));
    LayoutData<MFInterfaceFAB<Real> >       bco(grids[lvl]);
    LayoutData<MFInterfaceFAB<Real> >       weights(grids[lvl]);
    LayoutData<MFInterfaceFAB<VoFStencil> > mfstencil(grids[lvl]);

    MFLevelGrid mflg(m_mfis, eblg);
    for (DataIterator dit = mfstencil.dataIterator(); dit.ok(); ++dit){
      MFInterfaceFAB<VoFStencil>& sten = mfstencil[dit()];
      sten.define(mflg, dit());
    }

    LayoutData<IntVectSet> cfivs;
    EBArith::defineCFIVS(cfivs, grids[lvl], domains[lvl]);

    jump_bc* jump = new jump_bc(mflg, *m_bco_irreg[lvl], dx[lvl], 2, &cfivs);
    
  }

  // Allocate coefficients
  EBAMRIVData   ajump;
  EBAMRIVData   bjump;
  EBAMRIVData   cjump;

  Vector<RefCountedPtr<LayoutData<IntVectSet> > > jump_cells(1 + finest_level);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    jump_cells[lvl] = RefCountedPtr<LayoutData<IntVectSet> > (new LayoutData<IntVectSet> (grids[lvl]));
    const IntVectSet interface_ivs = m_mfis->interface_region(domains[lvl]);
    for (DataIterator dit = jump_cells[lvl]->dataIterator(); dit.ok(); ++dit){
      (*jump_cells[lvl])[dit()] = interface_ivs & grids[lvl].get(dit());
    }
  }



    
  Vector<MFLevelGrid> mfeblg(1 + finest_level);
  Vector<MFQuadCFInterp> mfquadcfi(1 + finest_level);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    Vector<EBLevelGrid> eblg_phases(phase::num_phases);
    Vector<RefCountedPtr<EBQuadCFInterp> > quadcfi_phases(phase::num_phases);

    
    eblg_phases[phase::gas] = *(m_amr->get_eblg(phase::gas)[lvl]);
    eblg_phases[phase::solid] = *(m_amr->get_eblg(phase::solid)[lvl]);

    quadcfi_phases[phase::gas] =   (m_amr->get_old_quadcfi(phase::gas)[lvl]);
    quadcfi_phases[phase::solid] = (m_amr->get_old_quadcfi(phase::solid)[lvl]);
    
    mfeblg[lvl].define(m_mfis, eblg_phases);
    mfquadcfi[lvl].define(quadcfi_phases);
  }

  m_opfact = RefCountedPtr<mf_helmholtz_opfactory> (new mf_helmholtz_opfactory(m_mfis,
									       mfeblg,
									       mfquadcfi,
									       refinement_ratios,
									       grids,
									       m_aco,
									       m_bco,
									       m_bco_irreg,
									       alpha,
									       beta,
									       dx[0],
									       domains[0],
  									       domfact,
									       ebcfact,
  									       origin,
  									       3*IntVect::Unit,
  									       3*IntVect::Unit));

  m_opfact->set_jump(1.0, 1.0);


  // Test the jump_bc stuff
  // 1. Allocate surface potential
  const int lvl = finest_level;
  EBAMRIVData phi_bc;
  m_amr->allocate(phi_bc, phase::gas, 1, 0);

  LayoutData<IntVectSet> cfivs;
  EBArith::defineCFIVS(cfivs, grids[lvl], domains[lvl]);
  jump_bc* jump = new jump_bc(mfeblg[lvl], *m_bco_irreg[lvl], dx[lvl], 2, &cfivs);

  EBLevelDataOps::setVal(*sigma[lvl], 1.0);
  MFLevelDataOps::setVal(*m_state[lvl], 1.0);

  jump->match_bc(*phi_bc[lvl], *sigma[lvl], *m_state[lvl]);
}
