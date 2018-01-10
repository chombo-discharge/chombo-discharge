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
#include "amr_mesh.H"
#include "domain_bc.H"
#include "poissonconductivitydomainbcfactory.H"

#include <Stencils.H>
#include <MFCellFAB.H>
#include <LayoutData.H>
#include <MFLevelDataOps.H>
#include <DirichletConductivityDomainBC.H>
#include <DirichletConductivityEBBC.H>

poisson_multifluid_gmg::poisson_multifluid_gmg(){
  m_needs_setup = true;

  this->set_gmg_solver_parameters();
  this->set_bottom_solver(0);
  this->set_botsolver_smooth(32);
  this->set_bottom_drop(16);
}

poisson_multifluid_gmg::~poisson_multifluid_gmg(){

}

bool poisson_multifluid_gmg::solve(const bool a_zerophi){
  CH_TIME("poisson_multifluid_gmg::solve");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::solve";
  }
  
  const bool converged = this->solve(m_state, m_source, a_zerophi);

  return converged;
}

bool poisson_multifluid_gmg::solve(MFAMRCellData& a_state, const MFAMRCellData& a_source, const bool a_zerophi){
  CH_TIME("poisson_multifluid_gmg::solve(mfamrcell, mfamrcell");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::solve(mfamrcell, mfamrcell)" << endl;
  }

  if(m_needs_setup){
    this->setup_gmg(); // This does everything, allocates coefficients, gets bc stuff and so on
  }

  const int finest_level = m_amr->get_finest_level();

  Vector<LevelData<MFCellFAB>* > phi, rhs, res;
  m_amr->alias(phi, a_state);
  m_amr->alias(rhs, a_source);
  m_amr->alias(res, m_resid);

  m_gmg_solver.init(phi, rhs, finest_level, 0);
  m_gmg_solver.solveNoInitResid(phi, res, rhs, finest_level, 0, a_zerophi);
  m_gmg_solver.revert(phi, rhs, finest_level, 0);

  m_amr->average_down(a_state);
  m_amr->interp_ghost(a_state);

  bool converged = false;

  const int status = m_gmg_solver.m_exitStatus;   // 1 => Initial norm sufficiently reduced
  if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
    converged = true;
  }

  return converged;
}

int poisson_multifluid_gmg::query_ghost() const {
  return 3; // Need this many cells
}

void poisson_multifluid_gmg::set_bottom_solver(const int a_whichsolver){
  if(a_whichsolver == 0 || a_whichsolver == 1){
    m_bottomsolver = a_whichsolver;
  }
  else{
    MayDay::Abort("poisson_multifluid_gmg::set_bottom_solver - Unsupported solver type requested");
  }
}

void poisson_multifluid_gmg::set_botsolver_smooth(const int a_numsmooth){
  CH_TIME("poisson_multifluid_gmg::set_botsolver_smooth");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::set_botsolver_smooth" << endl;
  }
  CH_assert(a_numsmooth > 0);
  
  m_numsmooth = a_numsmooth;
}

void poisson_multifluid_gmg::set_bottom_drop(const int a_bottom_drop){
  m_bottom_drop = a_bottom_drop;
}

void poisson_multifluid_gmg::set_gmg_solver_parameters(relax::which_relax a_relax_type,
						       amrmg::which_mg a_gmg_type,      
						       const int a_verbosity,          
						       const int a_pre_smooth,         
						       const int a_post_smooth,       
						       const int a_bot_smooth,         
						       const int a_max_iter,
						       const int a_min_iter,
						       const Real a_eps,               
						       const Real a_hang,              
						       const Real a_norm_thresh){
  CH_TIME("poisson_multifluid_gmg::set_gmg_solver_parameters");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::set_gmg_solver_parameters" << endl;
  }

  m_gmg_relax_type  = a_relax_type;
  m_gmg_type        = a_gmg_type;
  m_gmg_verbosity   = a_verbosity;
  m_gmg_pre_smooth  = a_pre_smooth;
  m_gmg_post_smooth = a_post_smooth;
  m_gmg_bot_smooth  = a_bot_smooth;
  m_gmg_max_iter    = a_max_iter;
  m_gmg_min_iter    = a_min_iter;
  m_gmg_eps         = a_eps;
  m_gmg_hang        = a_hang;
  m_gmg_norm_thresh = a_norm_thresh;
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

  data_ops::set_value(m_aco,       0.0);  // Always zero for poisson equation, but that is done from alpha. 
  data_ops::set_value(m_bco,       eps0); // Will override this later
  data_ops::set_value(m_bco_irreg, eps0); // Will override this later

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

      a_perm[dir](face, comp) = a_dielectrics[closest].get_permittivity(pos);
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

    a_perm(vof, comp) = a_dielectrics[closest].get_permittivity(pos);
  }
}

void poisson_multifluid_gmg::setup_gmg(){
  CH_TIME("poisson_multifluid_gmg::setup_gmg");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::setup_gmg" << endl;
  }
  // const int finest_level = m_amr->get_finest_level();
  
  // for (int lvl = 0; lvl <= finest_level; lvl++){
  //   const ProblemDomain dom = m_amr->get_domains()[lvl];
  //   const RefCountedPtr<EBIndexSpace>& ebis1 = m_mfis->get_ebis(phase::gas);
  //   const RefCountedPtr<EBIndexSpace>& ebis2 = m_mfis->get_ebis(phase::solid);
  //   const int whichLev = ebis1->getLevel(dom);
  //   pout() << "level = " << lvl << "\t domain = " << dom << "\t ebis1 irreg = " << ebis1->irregCells(whichLev).numPts() << endl;
  //   pout() << "level = " << lvl << "\t domain = " << dom << "\t ebis2 irreg = " << ebis2->irregCells(whichLev).numPts() << endl;
  //   pout() << "level = " << lvl << "\t domain = " << dom << "\t isect_cells = " << m_mfis->interface_region(dom).numPts() << endl;
  //   pout() << endl;
  // }
  
  this->set_coefficients();       // Set coefficients
  this->setup_operator_factory(); // Set the operator factory
  this->setup_solver();           // Set up the AMR multigrid solver
  
#if 0 // Testing stuff
  this->do_ebcond_test();
#endif
#if 0 // Testing stuff
  this->base_tests();
  MayDay::Abort("stop");
#endif
  


}

void poisson_multifluid_gmg::setup_operator_factory(){
  CH_TIME("poisson_multifluid_gmg::setup_operator_factory");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::setup_operator_factory" << endl;
  }

  const int finest_level                 = m_amr->get_finest_level();
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids();
  const Vector<int>& refinement_ratios   = m_amr->get_ref_rat();
  const Vector<ProblemDomain>& domains   = m_amr->get_domains();
  const Vector<Real>& dx                 = m_amr->get_dx();
  const RealVect& origin                 = m_physdom->get_prob_lo();


  // This stuff is needed for the operator factory
  Vector<MFLevelGrid>    mflg(1 + finest_level);
  Vector<MFQuadCFInterp> mfquadcfi(1 + finest_level);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    Vector<EBLevelGrid>                    eblg_phases(phase::num_phases);
    Vector<RefCountedPtr<EBQuadCFInterp> > quadcfi_phases(phase::num_phases);

    eblg_phases[phase::gas]   = *(m_amr->get_eblg(phase::gas)[lvl]);
    eblg_phases[phase::solid] = *(m_amr->get_eblg(phase::solid)[lvl]);

    quadcfi_phases[phase::gas]   = (m_amr->get_old_quadcfi(phase::gas)[lvl]);
    quadcfi_phases[phase::solid] = (m_amr->get_old_quadcfi(phase::solid)[lvl]);
    
    mflg[lvl].define(m_mfis, eblg_phases);
    mfquadcfi[lvl].define(quadcfi_phases);
  }

  // Appropriate coefficients for poisson equation
  const Real alpha =  0.0;
  const Real beta  = -1.0;

  RefCountedPtr<BaseDomainBCFactory> domfact = RefCountedPtr<BaseDomainBCFactory> (NULL);

  const IntVect ghost_phi = this->query_ghost()*IntVect::Unit;
  const IntVect ghost_rhs = this->query_ghost()*IntVect::Unit;

  poissonconductivitydomainbcfactory* bcfact = new poissonconductivitydomainbcfactory();
  bcfact->set_value(0.0);
  domfact = RefCountedPtr<BaseDomainBCFactory> (bcfact);
  

  m_opfact = RefCountedPtr<mfconductivityopfactory> (new mfconductivityopfactory(m_mfis,
									       mflg,
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
									       origin,
									       ghost_phi,
									       ghost_rhs,
									       1 + finest_level));

  m_opfact->set_relax_type(2);
  m_opfact->set_bottom_drop(m_bottom_drop);
  m_opfact->set_jump(0.0, -1.0);
  m_opfact->set_electrodes(m_compgeom->get_electrodes());

  domain_bc dombc(m_wallbc);
}

void poisson_multifluid_gmg::setup_solver(){
  CH_TIME("poisson_multifluid_gmg::setup_solver");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::setup_solver" << endl;
  }

  const int finest_level       = m_amr->get_finest_level();
  const ProblemDomain coar_dom = m_amr->get_domains()[0];

  // Select bottom solver
  LinearSolver<LevelData<MFCellFAB> >* botsolver = NULL;
  if(m_bottomsolver == 0){
    m_mfsolver.setNumSmooths(m_numsmooth);
    botsolver = &m_mfsolver;
  }
  else{
    botsolver = &m_bicgstab;
  }

  m_gmg_solver.define(coar_dom, *m_opfact, botsolver, 1 + finest_level);
  m_gmg_solver.setSolverParameters(m_gmg_pre_smooth,
				   m_gmg_post_smooth,
				   m_gmg_bot_smooth,
				   m_gmg_type,
				   m_gmg_max_iter,
				   m_gmg_eps,
				   m_gmg_hang,
				   m_gmg_norm_thresh);
  m_gmg_solver.m_imin = m_gmg_min_iter;
  m_gmg_solver.m_verbosity = m_gmg_verbosity;
}
