/*!
  @file cdr_solver.cpp
  @brief Implementation of cdr_solver.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "cdr_solver.H"
#include "cdr_solverF_F.H"
#include "data_ops.H"

#include <EBAMRIO.H>
#include <EBArith.H>

cdr_solver::cdr_solver(){
  m_name = "cdr_solver";

  this->set_verbosity(-1);
  this->set_phase(phase::gas);
  this->set_time(0, 0., 0.);
  this->set_mass_redist(false);
}

cdr_solver::~cdr_solver(){

}

std::string cdr_solver::get_name(){
  return m_name;
}

int cdr_solver::query_ghost() const {
  CH_TIME("cdr_solver::query_ghost");
  if(m_verbosity > 5){
    pout() << m_name + "::query_ghost" << endl;
  }

  return 3;
}

void cdr_solver::advance(const Real& a_dt){
  CH_TIME("cdr_solver::advance(dt)");
  if(m_verbosity > 1){
    pout() << m_name + "::advance(dt)" << endl;
  }
  
  this->advance(m_state, a_dt);
}

void cdr_solver::advance(EBAMRCellData& a_state, const Real& a_dt){
  CH_TIME("cdr_solver::advance(state, dt)");
  if(m_verbosity > 1){
    pout() << m_name + "::advance(state, dt)" << endl;
  }

  const Real midpoint = 0.5;
  this->advance_rk2(a_state, a_dt, midpoint);

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::floor(*a_state[lvl], 0.0); // This adds mass!

    a_state[lvl]->exchange();
  }
  
  m_amr->average_down(a_state, m_phase);
  m_amr->interp_ghost(a_state, m_phase);

  m_time += a_dt;
  m_step++;
}

void cdr_solver::advance_rk2(EBAMRCellData& a_state, const Real a_dt, const Real a_alpha){
  CH_TIME("cdr_solver::advance_rk2");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_rk2" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  // Compute right-hand-side
  EBAMRCellData k1, k2, phi;
  m_amr->allocate(k1,  m_phase, ncomp);
  m_amr->allocate(k2,  m_phase, ncomp);
  m_amr->allocate(phi, m_phase, ncomp);

  // Compute f(yn, tn);
  this->compute_rhs(k1, a_state, a_dt);

  //  phi = yn + a_alpha*dt*f(tn, yn)
  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_state[lvl]->copyTo(*phi[lvl]);
    data_ops::incr(*phi[lvl], *k1[lvl], a_alpha*a_dt);
  }

  // Compute f(tn + alpha*dt, yn + alpha*dt*f(yn, tn)) = f(tn + alpha*dt, phi)
  this->compute_rhs(k2, phi, a_dt);

  // Make y(n+1) = y(n) 
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::incr(*a_state[lvl], *k1[lvl], (1. - 1./(2*a_alpha))*a_dt);
    data_ops::incr(*a_state[lvl], *k2[lvl], 1./(2*a_alpha)*a_dt);
  }
}

void cdr_solver::allocate_internals(){
  CH_TIME("cdr_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  CH_assert(m_phase == phase::gas);

  const int sca = 1;
  const int vec = SpaceDim;

  m_amr->allocate(m_state,      m_phase, sca);
  m_amr->allocate(m_source,     m_phase, sca);
  m_amr->allocate(m_velo_face,  m_phase, sca);
  m_amr->allocate(m_velo_cell,  m_phase, vec); 
  m_amr->allocate(m_ebflux,     m_phase, sca);
  m_amr->allocate(m_diffco,     m_phase, sca);
  m_amr->allocate(m_diffco_eb,  m_phase, sca);
  m_amr->allocate(m_scratch,       m_phase, sca);

  data_ops::set_value(m_state,      0.0);
  data_ops::set_value(m_source,     0.0);
  data_ops::set_value(m_velo_face,  0.0);
  data_ops::set_value(m_velo_cell,  0.0);
  data_ops::set_value(m_ebflux,     0.0);
  data_ops::set_value(m_diffco,     0.0);
  data_ops::set_value(m_diffco_eb,  0.0);
  data_ops::set_value(m_scratch,    0.0);

  this->define_interp_stencils();
  this->define_divFnc_stencils();
}

void cdr_solver::average_velo_to_faces(EBAMRFluxData& a_velo_face, const EBAMRCellData& a_velo_cell){
  CH_TIME("cdr_solver::average_velo_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::average_velo_to_faces" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::average_cell_to_face(*a_velo_face[lvl], *a_velo_cell[lvl], m_amr->get_domains()[lvl]);
    a_velo_face[lvl]->exchange();
  }
}

void cdr_solver::cache_state(){
  CH_TIME("cdr_solver::cache_state");
  if(m_verbosity > 5){
    pout() << m_name + "::cache_state" << endl;
  }

  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();
  
  m_amr->allocate(m_cache, m_phase, ncomp);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    m_state[lvl]->copyTo(*m_cache[lvl]);
  }

}

void cdr_solver::coarse_fine_increment(const EBAMRIVData& a_mass_diff){
  CH_TIME("cdr_solver::coarse_fine_increment");
  if(m_verbosity > 5){
    pout() << m_name + "::coarse_fine_increment" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(0,0);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];

    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->get_fine_to_coar_redist(m_phase)[lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->get_coar_to_fine_redist(m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->get_coar_to_coar_redist(m_phase)[lvl];

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < 0;

    if(has_coar){
      fine2coar_redist->setToZero();

    }
    if(has_fine){
      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      if(has_coar){
	fine2coar_redist->increment((*a_mass_diff[lvl])[dit()], dit(), interv);
      }

      if(has_fine){
	coar2fine_redist->increment((*a_mass_diff[lvl])[dit()], dit(), interv);
	coar2coar_redist->increment((*a_mass_diff[lvl])[dit()], dit(), interv);
      }
    }
  }
}

void cdr_solver::coarse_fine_redistribution(EBAMRCellData& a_state){
  CH_TIME("cdr_solver::coarse_fine_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::coarse_fine_redistribution" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->get_dx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->get_coar_to_fine_redist(m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->get_coar_to_coar_redist(m_phase)[lvl];
    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->get_fine_to_coar_redist(m_phase)[lvl];
    
    if(has_coar){
      fine2coar_redist->redistribute(*a_state[lvl-1], interv);
      fine2coar_redist->setToZero();
    }

    if(has_fine){
      coar2fine_redist->redistribute(*a_state[lvl+1], interv);
      coar2coar_redist->redistribute(*a_state[lvl],   interv);

      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }
  }
}

void cdr_solver::compute_divF_irreg(LevelData<EBCellFAB>&              a_divF,
				    const LevelData<BaseIFFAB<Real> >  a_flux[SpaceDim],
				    const LevelData<BaseIVFAB<Real> >& a_ebflux,
				    const int                          a_lvl){
  CH_TIME("cdr_solver::compute_divF_irreg");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divF_irreg" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];
  const Real dx                = m_amr->get_dx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);

    EBCellFAB& divF               = a_divF[dit()];
    const BaseIVFAB<Real>& ebflux = a_ebflux[dit()];

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real area     = ebisbox.bndryArea(vof);

      // EB flux
      divF(vof, comp) = ebflux(vof,comp)*area;

      // Face fluxes
      for (int dir = 0; dir < SpaceDim; dir++){
	const BaseIFFAB<Real>& flux = a_flux[dir][dit()];

	for (SideIterator sit; sit.ok(); ++sit){
	  const int isign = sign(sit());
	  const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit());

	  for (int iface = 0; iface < faces.size(); iface++){
	    const FaceIndex face = faces[iface];
	    const Real face_area = ebisbox.areaFrac(face);
	    divF(vof, comp) += isign*face_area*flux(face, comp);
	  }
	}
      }

      // Scale divF by dx but not by kappa.
      divF(vof, comp) *= 1./dx;
    }
  }
}

void cdr_solver::compute_flux(EBAMRFluxData& a_flux, const EBAMRFluxData& a_face_state, const EBAMRFluxData& a_face_vel){
  CH_TIME("cdr_solver::compute_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_flux" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& flx       = (*a_flux[lvl])[dit()][dir];
	const EBFaceFAB& phi = (*a_face_state[lvl])[dit()][dir];
	const EBFaceFAB& vel = (*a_face_vel[lvl])[dit()][dir];

	const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingWithBoundary;
	for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();
	  flx(face, comp) = vel(face, comp)*phi(face, comp);
	}
      }
    }
  }
						
}

void cdr_solver::compute_rhs(EBAMRCellData& a_rhs, const Real& a_dt){
  CH_TIME("cdr_solver::compute_rhs(rhs, dt)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_rhs(rhs, dt)" << endl;
  }
  
  this->compute_rhs(a_rhs, m_state, a_dt);
}

void cdr_solver::compute_rhs(EBAMRCellData& a_rhs, const EBAMRCellData& a_state, const Real& a_dt){
  CH_TIME("cdr_solver::compute_rhs(rhs, state, dt)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_rhs(rhs, state, dt)" << endl;
  }
  
  this->compute_divJ(a_rhs, a_state, a_dt); // a_rhs = div(J)
  data_ops::scale(a_rhs, -1.0);             // a_rhs = -div(J)
  data_ops::incr(a_rhs, m_source, 1.0);     // a_rhs = div(J) + S

  m_amr->average_down(a_rhs, m_phase);
  m_amr->interp_ghost(a_rhs, m_phase);
}

void cdr_solver::conservative_divergence(EBAMRCellData& a_cons_div, const EBAMRFluxData& a_flux){
  CH_TIME("cdr_solver::conservative_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::conservative_divergence" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

    // Compute div(F) on regular cells
    this->consdiv_regular(*a_cons_div[lvl], *a_flux[lvl], lvl);

    LevelData<BaseIFFAB<Real> > flux[SpaceDim];
    this->setup_flux_interpolant(flux, *a_flux[lvl], lvl);
    this->interpolate_flux_to_centroids(flux, lvl);
    this->compute_divF_irreg(*a_cons_div[lvl], flux, *m_ebflux[lvl], lvl);

    a_cons_div[lvl]->exchange();
  }
}

void cdr_solver::conservative_divergence(EBAMRCellData&       a_cons_div,
					 const EBAMRFluxData& a_face_vel,
					 const EBAMRFluxData& a_face_state){
  CH_TIME("cdr_solver::conservative_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::conservative_divergence" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  // Use some transient storage
  EBAMRFluxData flux;
  m_amr->allocate(flux, m_phase, ncomp);

  this->compute_flux(flux, a_face_vel, a_face_state);
  this->conservative_divergence(a_cons_div, flux);
}

void cdr_solver::consdiv_regular(LevelData<EBCellFAB>& a_divJ, const LevelData<EBFluxFAB>& a_flux, const int a_lvl){
  CH_TIME("cdr_solver::consdiv_regular");
  if(m_verbosity > 5){
    pout() << m_name + "::consdiv_regular" << endl;
  }

  CH_assert(a_divJ.nComp() == 1);
  CH_assert(a_flux.nComp() == 1);

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain domain  = m_amr->get_domains()[a_lvl];
  const Real dx               = m_amr->get_dx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& divJ         = a_divJ[dit()];
    BaseFab<Real>& divJ_fab = divJ.getSingleValuedFAB();
    const Box box = dbl.get(dit());

    for (int dir = 0; dir < SpaceDim; dir++){
      const EBFaceFAB& flx         = a_flux[dit()][dir];
      const BaseFab<Real>& flx_fab = flx.getSingleValuedFAB();

      FORT_CONSDIV_REG(CHF_FRA1(divJ_fab, comp),
		       CHF_CONST_FRA1(flx_fab, comp),
		       CHF_CONST_INT(dir),
		       CHF_CONST_INT(ncomp),
		       CHF_CONST_REAL(dx),
		       CHF_BOX(box));
    }


    // Reset irregular cells - these contain bogus values and will be set elsewhere. 
    const EBISBox& ebisbox = divJ.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      divJ(vof, comp) = 0.0;
    }
  }

  a_divJ.exchange();
}

void cdr_solver::define_divFnc_stencils(){
  CH_TIME("cdr_solver::define_divFnc_stencils");
  if(m_verbosity > 5){
    pout() << m_name + "::define_divFnc_stencils" << endl;
  }

  const int rad   = 1;
  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  m_stencils_nc.resize(1 + finest_level);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

    m_stencils_nc[lvl] = RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > (new LayoutData<BaseIVFAB<VoFStencil> >(dbl));
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = ebisbox.getIrregIVS(box);

      BaseIVFAB<VoFStencil>& stens = (*m_stencils_nc[lvl])[dit()];
      stens.define(ivs, ebgraph, ncomp);
      
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	VoFStencil& sten = stens(vof, comp);
	sten.clear();
	
	Real norm = 0.;
	Vector<VolIndex> vofs;
	EBArith::getAllVoFsInMonotonePath(vofs, vof, ebisbox, rad);
	for (int i = 0; i < vofs.size(); i++){
	  const VolIndex& ivof = vofs[i];
	  const Real iweight   = ebisbox.volFrac(ivof);

	  norm += iweight;
	  sten.add(ivof, 1.0);
	}

	sten *= 1./norm;
      }
    }
  }
}

void cdr_solver::define_interp_stencils(){
  CH_TIME("cdr_solver::define_interp_stencils");
  if(m_verbosity > 5){
    pout() << m_name + "::define_interp_stencils" << endl;
  }

  const int comp            = 0;
  const int ncomp           = 1;
  const int finest_level    = m_amr->get_finest_level();
  FaceStop::WhichFaces stop = FaceStop::SurroundingWithBoundary;
  
  for (int dir = 0; dir < SpaceDim; dir++){
    (m_interp_stencils[dir]).resize(1 + finest_level);

    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const ProblemDomain& domain  = m_amr->get_domains()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

      m_interp_stencils[dir][lvl] = RefCountedPtr<LayoutData<BaseIFFAB<FaceStencil> > >
	(new LayoutData<BaseIFFAB<FaceStencil> >(dbl));

      LayoutData<IntVectSet> cfivs(dbl);
      EBArith::defineCFIVS(cfivs, dbl, domain);

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	BaseIFFAB<FaceStencil>& sten = (*m_interp_stencils[dir][lvl])[dit()];
	const Box box          = dbl.get(dit());
	const EBISBox& ebisbox = ebisl[dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();
	const IntVectSet ivs   = ebisbox.getIrregIVS(box);

	sten.define(ivs, ebisbox.getEBGraph(), dir, ncomp);
	
	for (FaceIterator faceit(ivs, ebgraph, dir, stop); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();
	  FaceStencil& facesten = sten(face,comp);

	  facesten = EBArith::getInterpStencil(face, cfivs[dit()], ebisbox, domain.domainBox());
	}
      }
    }
  }
}

void cdr_solver::increment_redist_flux(){
  CH_TIME("cdr_solver::increment_redist_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_redist_flux" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->get_dx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;
    
    if(has_fine){
      RefCountedPtr<EBFluxRegister>& fluxreg = m_amr->get_flux_reg(m_phase)[lvl];
      RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->get_coar_to_fine_redist(m_phase)[lvl];
      RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->get_coar_to_coar_redist(m_phase)[lvl];
      
      const Real scale = -dx;

      fluxreg->incrementRedistRegister(*coar2fine_redist, interv, scale);
      fluxreg->incrementRedistRegister(*coar2coar_redist, interv, scale);

    }
  }
}

void cdr_solver::initial_data(){
  CH_TIME("cdr_solver::initial_data");
  if(m_verbosity > 5){
    pout() << m_name + "::initial_data" << endl;
  }

  const RealVect origin  = m_physdom->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (DataIterator dit = m_state[lvl]->dataIterator(); dit.ok(); ++dit){
      EBCellFAB& state       = (*m_state[lvl])[dit()];
      const Box box          = m_state[lvl]->disjointBoxLayout().get(dit());
      const EBISBox& ebisbox = state.getEBISBox();
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = EBArith::getVofLocation(vof, m_amr->get_dx()[lvl]*RealVect::Unit, origin);
	
	for (int comp = 0; comp < state.nComp(); comp++){
	  state(vof, comp) = m_species->initial_data(pos, m_time);
	}
      }
    }
  }

  m_amr->average_down(m_state, m_phase);
  m_amr->interp_ghost(m_state, m_phase);
}

void cdr_solver::hybrid_divergence(EBAMRCellData&     a_hybrid_div,
				   EBAMRIVData&       a_mass_diff,
				   const EBAMRIVData& a_noncons_div){
  CH_TIME("cdr_solver::hybrid_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::hybrid_divergence" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      EBCellFAB& divH               = (*a_hybrid_div[lvl])[dit()];  // On input, this contains kappa*div(F)
      BaseIVFAB<Real>& deltaM       = (*a_mass_diff[lvl])[dit()];
      const BaseIVFAB<Real>& divNC  = (*a_noncons_div[lvl])[dit()]; 

      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = ebisbox.getIrregIVS(box);

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real kappa    = ebisbox.volFrac(vof);
	const Real dc       = divH(vof, comp);
	const Real dnc      = divNC(vof, comp);

	divH(vof, comp)   = dc + (1-kappa)*dnc;          // On output, contains hybrid divergence
	deltaM(vof, comp) = (1-kappa)*(dc - kappa*dnc);
      }
    }
  }
}

void cdr_solver::hyperbolic_redistribution(EBAMRCellData&       a_divF,
					   const EBAMRIVData&   a_mass_diff,
					   const EBAMRCellData& a_redist_weights){
  CH_TIME("cdr_solver::hyberbolic_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::hyperbolic_redistribution" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);


  for (int lvl = 0; lvl <= finest_level; lvl++){
    EBLevelRedist& level_redist = *(m_amr->get_level_redist(m_phase)[lvl]);
    if(m_mass_redist){
      level_redist.resetWeights(*a_redist_weights[lvl], comp);
    }
    level_redist.redistribute(*a_divF[lvl], interv);
    level_redist.setToZero();
  }
}

void cdr_solver::interpolate_flux_to_centroids(LevelData<BaseIFFAB<Real> >       a_flux[SpaceDim],
					       const int                         a_lvl){
  CH_TIME("cdr_solver::interpolate_flux_to_centroids");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_flux_to_centroids" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);
  
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];

  const FaceStop::WhichFaces stop = FaceStop::SurroundingWithBoundary;

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);
    
    for (int dir = 0; dir < SpaceDim; dir++){
      BaseIFFAB<Real>& flux = a_flux[dir][dit()];
      BaseIFFAB<Real> centroid_flux;

      // Compute face centroid flux
      centroid_flux.define(ivs, ebgraph, dir, ncomp);
      for (FaceIterator faceit(ivs, ebgraph, dir, stop); faceit.ok(); ++faceit){
	const FaceIndex& face   = faceit();
	const FaceStencil& sten = (*m_interp_stencils[dir][a_lvl])[dit()](face, comp);

	centroid_flux(face, comp) = 0.;
	for (int i = 0; i < sten.size(); i++){
	  const FaceIndex& iface = sten.face(i);
	  const Real iweight     = sten.weight(i);
	  
	  centroid_flux(face, comp) += iweight*flux(iface, comp);
	}
      }

      // Copy centroid flux into a_flux
      flux.setVal(0.0);
      flux.copy(box, interv, box, centroid_flux, interv);
    }
  }
}

void cdr_solver::increment_flux_register(const EBAMRFluxData& a_face_state, const EBAMRFluxData& a_velo_face){
  CH_TIME("cdr_solver::increment_flux_register");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_flux_register" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  Vector<RefCountedPtr<EBFluxRegister> >& fluxreg = m_amr->get_flux_reg(m_phase);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    if(has_fine){
      fluxreg[lvl]->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      const Box box          = dbl.get(dit());
      
      for (int dir = 0; dir < SpaceDim; dir++){
	const Real scale     = 1.;
	const EBFaceFAB& phi = (*a_face_state[lvl])[dit()][dir];
	const EBFaceFAB& vel = (*a_velo_face[lvl])[dit()][dir];

	// Compute flux
	EBFaceFAB flux(ebisbox, box, dir, ncomp);
	flux.setVal(0.0);
	flux += phi;
	flux *= vel;

	// Increment flux register for irregular/regular. Add both from coarse to fine and from fine to coarse
	if(has_fine){
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxreg[lvl]->incrementCoarseBoth(flux, scale, dit(), interv, dir, sit());
	  }
	}
	if(has_coar){
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxreg[lvl-1]->incrementFineBoth(flux, scale, dit(), interv, dir, sit());
	  }
	}
      }
    }
  }
}

void cdr_solver::increment_flux_register(const EBAMRFluxData& a_flux){
  CH_TIME("cdr_solver::increment_flux_register(flux)");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_flux_register(flux)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  Vector<RefCountedPtr<EBFluxRegister> >& fluxreg = m_amr->get_flux_reg(m_phase);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    if(has_fine){
      fluxreg[lvl]->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      const Box box          = dbl.get(dit());
      
      for (int dir = 0; dir < SpaceDim; dir++){
	const Real scale      = 1.0;
	const EBFaceFAB& flux = (*a_flux[lvl])[dit()][dir];

	// Increment flux register for irregular/regular. Add both from coarse to fine and from fine to coarse
	if(has_fine){
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxreg[lvl]->incrementCoarseBoth(flux, scale, dit(), interv, dir, sit());
	  }
	}
	if(has_coar){
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxreg[lvl-1]->incrementFineBoth(flux, scale, dit(), interv, dir, sit());
	  }
	}
      }
    }
  }
}

void cdr_solver::increment_redist(const EBAMRIVData& a_mass_diff){
  CH_TIME("cdr_solver::increment_redist");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_redist" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    
    EBLevelRedist& level_redist = *(m_amr->get_level_redist(m_phase)[lvl]);
    level_redist.setToZero();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      level_redist.increment((*a_mass_diff[lvl])[dit()], dit(), interv);
    }
  }
}

void cdr_solver::nonconservative_divergence(EBAMRIVData&         a_div_nc,
					    const EBAMRCellData& a_divF,
					    const EBAMRFluxData& a_face_state){
  CH_TIME("cdr_solver::nonconservative_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::nonconservative_divergence" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    
    for (DataIterator dit = a_div_nc[lvl]->dataIterator(); dit.ok(); ++dit){
      BaseIVFAB<Real>& div_nc = (*a_div_nc[lvl])[dit()];
      EBCellFAB& divF         = (*a_divF[lvl])[dit()];

      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = ebisbox.getIrregIVS(box);

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof    = vofit();
	const VoFStencil& sten = (*m_stencils_nc[lvl])[dit()](vof, comp);

	div_nc(vof, comp) = 0.;
	for (int i = 0; i < sten.size(); i++){
	  const VolIndex& ivof = sten.vof(i);
	  const Real& iweight  = sten.weight(i);
	  div_nc(vof,comp) += iweight*divF(ivof, comp);
	}
      }
    }
  }
}



void cdr_solver::regrid(const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("cdr_solver::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);

  this->allocate_internals();

  Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->get_eb_pwl_interp(m_phase);

  m_cache[0]->copyTo(*m_state[0]); // Base level should never change. 
  for (int lvl = 1; lvl <= a_new_finest_level; lvl++){
    interpolator[lvl]->interpolate(*m_state[lvl], *m_state[lvl-1], interv);
    if(lvl <= Min(a_old_finest_level, a_new_finest_level)){
      m_cache[lvl]->copyTo(*m_state[lvl]);
    }
  }
}

void cdr_solver::reflux(EBAMRCellData& a_state){
  CH_TIME("cdr_solver::reflux");
  if(m_verbosity > 5){
    pout() << m_name + "::reflux" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  Vector<RefCountedPtr<EBFluxRegister > >& fluxreg = m_amr->get_flux_reg(m_phase);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->get_dx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    if(has_fine){
      const Real scale = 1.0/dx;
      
      fluxreg[lvl]->reflux(*a_state[lvl], interv, scale);
      fluxreg[lvl]->setToZero();
    }
  }
}

void cdr_solver::sanity_check(){
  CH_TIME("cdr_solver::sanity_check");
  if(m_verbosity > 5){
    pout() << m_name + "::sanity_check" << endl;
  }

  CH_assert(!m_compgeom.isNull());
  CH_assert(!m_physdom.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_species.isNull());
  CH_assert(!m_ebis.isNull());
}

void cdr_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("cdr_solver::set_amr");
  if(m_verbosity > 5){
    pout() << m_name + "::set_amr" << endl;
  }

  m_amr = a_amr;
}

void cdr_solver::set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom){
  CH_TIME("cdr_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << m_name + "::set_computational_geometry" << endl;
  }
  m_compgeom = a_compgeom;

  const RefCountedPtr<mfis> mfis = m_compgeom->get_mfis();
  
  this->set_ebis(mfis->get_ebis(m_phase));
}

void cdr_solver::set_diffco(const EBAMRFluxData& a_diffco, const EBAMRIVData& a_diffco_eb){
  CH_TIME("cdr_solver::set_diffco");
  if(m_verbosity > 5){
    pout() << m_name + "::set_diffco" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_diffco[lvl]->copyTo(*m_diffco[lvl]);
    a_diffco_eb[lvl]->copyTo(*m_diffco_eb[lvl]);
  }

  m_amr->average_down(m_diffco, m_phase);
  m_amr->average_down(m_diffco_eb, m_phase);
}

void cdr_solver::set_diffco(const Real a_diffco){
  CH_TIME("cdr_solver::set_diffco");
  if(m_verbosity > 5){
    pout() << m_name + "::set_diffco" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::set_value(*m_diffco[lvl],    a_diffco);
    data_ops::set_value(*m_diffco_eb[lvl], a_diffco);

    m_diffco[lvl]->exchange();
  }

  m_amr->average_down(m_diffco, m_phase);
  m_amr->average_down(m_diffco_eb, m_phase);
}

void cdr_solver::set_ebflux(const EBAMRIVData& a_ebflux){
  CH_TIME("cdr_solver::set_ebflux(variable)");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebflux(variable)" << endl;
  }

  const int comp         = 0;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_ebflux[lvl]->copyTo(interv, *m_ebflux[lvl], interv);
  }
}

void cdr_solver::set_ebflux(const Real a_ebflux){
  CH_TIME("cdr_solver::set_ebflux(constant)");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebflux(constant)" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::set_value(*m_ebflux[lvl], a_ebflux);
  }
}

void cdr_solver::set_ebis(const RefCountedPtr<EBIndexSpace>& a_ebis){
  CH_TIME("cdr_solver::set_ebis");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebis" << endl;
  }

  m_ebis = a_ebis;
}

void cdr_solver::set_mass_redist(const bool a_mass_redist){
  CH_TIME("cdr_solver::set_mass_redist");
  if(m_verbosity > 5){
    pout() << m_name + "::set_mass_redist" << endl;
  }

  m_mass_redist = a_mass_redist;
}

void cdr_solver::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("cdr_solver::set_physical_domain");
  if(m_verbosity > 5){
    pout() << m_name + "::set_physical_domain" << endl;
  }

  m_physdom = a_physdom;
}

void cdr_solver::set_species(const RefCountedPtr<species> a_species){
  CH_TIME("cdr_solver::set_species");
  if(m_verbosity > 5){
    pout() << m_name + "::set_species" << endl;
  }

  m_species   = a_species;
  m_name      = m_species->get_name();
  m_diffusive = m_species->is_diffusive();
}

void cdr_solver::set_source(const EBAMRCellData& a_source){
  CH_TIME("cdr_solver::set_source");
  if(m_verbosity > 5){
    pout() << m_name + "::set_source" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_source[lvl]->copyTo(*m_source[lvl]);
  }

  m_amr->average_down(m_source, m_phase);
  m_amr->interp_ghost(m_source, m_phase);
}

void cdr_solver::set_source(const Real a_source){
  CH_TIME("cdr_solver::set_source");
  if(m_verbosity > 5){
    pout() << m_name + "::set_source" << endl;
  }

  const int comp = 0;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::set_value(*m_source[lvl], a_source, comp);
  }

  m_amr->average_down(m_source, m_phase);
  m_amr->interp_ghost(m_source, m_phase);
}

void cdr_solver::set_time(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("cdr_solver::set_time");
  if(m_verbosity > 5){
    pout() << m_name + "::set_time" << endl;
  }

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void cdr_solver::set_velocity(const EBAMRCellData& a_velo){
  CH_TIME("cdr_solver::set_velocity");
  if(m_verbosity > 5){
    pout() << m_name + "::set_velocity" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_velo[lvl]->copyTo(*m_velo_cell[lvl]);
  }

  m_amr->average_down(m_velo_cell, m_phase);
  m_amr->interp_ghost(m_velo_cell, m_phase);
}

void cdr_solver::set_velocity(const RealVect a_velo){
  CH_TIME("cdr_solver::set_velocity");
  if(m_verbosity > 5){
    pout() << m_name + "::set_velocity" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (int dir = 0; dir < SpaceDim; dir++){
      data_ops::set_value(*m_velo_cell[lvl], a_velo[dir], dir);
    }

    m_velo_cell[lvl]->exchange();
  }

  m_amr->average_down(m_velo_cell, m_phase);
  m_amr->interp_ghost(m_velo_cell, m_phase);
}

void cdr_solver::set_phase(const phase::which_phase a_phase){
  CH_TIME("cdr_solver::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::set_phase" << endl;
  }

  m_phase = a_phase;
}

void cdr_solver::set_verbosity(const int a_verbosity){
  CH_TIME("cdr_solver::set_verbosity");
  m_verbosity = a_verbosity;
  
  if(m_verbosity > 5){
    pout() << m_name + "::set_verbosity" << endl;
  }
}

void cdr_solver::setup_flux_interpolant(LevelData<BaseIFFAB<Real> >   a_interpolant[SpaceDim],
					const LevelData<EBFluxFAB>&   a_flux,
					const int                     a_lvl){
  CH_TIME("cdr_solver::setup_flux_interpolant");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_flux_interpolant" << endl;
  }

  const int ncomp = 1;
  const int comp  = 0;
  
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];

  for (int dir = 0; dir < SpaceDim; dir++){
    // Define interpolant
    LayoutData<IntVectSet> grown_set;
    EBArith::defineFluxInterpolant(a_interpolant[dir],
				   grown_set,
				   dbl,
				   ebisl,
				   domain,
				   ncomp,
				   dir);

    // Compute interpolant
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      BaseIFFAB<Real>& interpol   = a_interpolant[dir][dit()];
      const Box& box              = dbl.get(dit());
      const EBISBox& ebisbox      = ebisl[dit()];
      const EBGraph& ebgraph      = ebisbox.getEBGraph();
      const EBFaceFAB& flux       = a_flux[dit()][dir];
      const IntVectSet ivs        = grown_set[dit()] & box;
      const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingWithBoundary;

      for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	const FaceIndex& face = faceit();
	interpol(face, comp) = flux(face, comp);
      }
    }

    a_interpolant[dir].exchange();
  }
}

#ifdef CH_USE_HDF5
void cdr_solver::write_plot_file(){
  CH_TIME("cdr_solver::write_plot_file");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_file" << endl;
  }

  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_step, SpaceDim);

  const int ncomps = 2 + SpaceDim;
  Vector<string> names(ncomps);
  names[0] = "density";
  names[1] = "source";
  names[2] = "x-velocity";
  names[3] = "y-velocity";
  if(SpaceDim == 3){
    names[4] = "z-velocity";
  }


  EBAMRCellData output;
  m_amr->allocate(output, m_phase, ncomps);

  for (int lvl = 0; lvl < output.size(); lvl++){
    LevelData<EBCellFAB>& state  = *m_state[lvl];
    LevelData<EBCellFAB>& source = *m_source[lvl];
    LevelData<EBCellFAB>& velo   = *m_velo_cell[lvl];

    state.copyTo(Interval(0,0),           *output[lvl],  Interval(0,0));
    source.copyTo(Interval(0,0),          *output[lvl],  Interval(1,1));
    velo.copyTo(Interval(0,SpaceDim - 1), *output[lvl],  Interval(2, 2 + (SpaceDim-1)));
  }

#if 0 // Possibly removed. 
  // Transform to centroids
  irreg_amr_stencil<centroid_interp>& sten = m_amr->get_centroid_interp_stencils(phase::gas);
  sten.apply(output, true);
#endif


  Vector<LevelData<EBCellFAB>* > output_ptr;
  m_amr->alias(output_ptr, output);

#if 0 // Should be removed
  data_ops::floor(output, 0.0);
#endif
  
  Vector<Real> covered_values(ncomps, 0.0);
  string fname(file_char);
  writeEBHDF5(fname,
	      m_amr->get_grids(),
	      output_ptr,
	      names,
	      m_amr->get_domains()[0].domainBox(),
	      m_amr->get_dx()[0],
	      m_dt,
	      m_time,
	      m_amr->get_ref_rat(),
	      m_amr->get_finest_level() + 1,
	      true,
	      covered_values);
}
#endif

Real cdr_solver::compute_cfl_dt(){
  CH_TIME("cdr_solver::compute_cfl_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_cfl_dt" << endl;
  }

  Real min_dt = 1.E99;

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBCellFAB& velo  = (*m_velo_cell[lvl])[dit()];
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex vof = vofit();
	const RealVect u  = RealVect(D_DECL(velo(vof, 0), velo(vof, 1), velo(vof, 2)));
	const Real thisdt = dx/u.vectorLength();

	min_dt = Min(min_dt, thisdt);
      }
    }
  }
  
#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&min_dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("cdr_solver::compute_cfl_dt() - communication error on norm");
  }
  min_dt = tmp;
#endif


  return min_dt;
}

Real cdr_solver::compute_diffusive_dt(){
  CH_TIME("cdr_solver::compute_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_dt" << endl;
  }

  Real min_dt = 1.E99;

  if(m_diffusive){
    const int comp  = 0;
    const int ncomp = 1;
    const int finest_level = m_amr->get_finest_level();
    const Interval interv(comp, comp);

    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
      const Real dx                = m_amr->get_dx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const EBCellFAB& velo     = (*m_velo_cell[lvl])[dit()];
	const Box box             = dbl.get(dit());
	const EBISBox& ebisbox    = ebisl[dit()];
	const EBGraph& ebgraph    = ebisbox.getEBGraph();
	const IntVectSet& ivs     = ebisbox.getIrregIVS(box);
	FaceStop::WhichFaces stop = FaceStop::SurroundingWithBoundary;
	const IntVectSet norm_ivs(box);

	// Regular faces
	for (int dir = 0; dir < SpaceDim; dir++){
	  const EBFaceFAB& diffco = (*m_diffco[lvl])[dit()][dir];

	  for (FaceIterator faceit(norm_ivs, ebgraph, dir, stop); faceit.ok(); ++faceit){
	    const FaceIndex& face = faceit();
	    const Real thisdt = dx*dx/(2*SpaceDim*diffco(face, comp));

	    min_dt = Min(thisdt, min_dt);
	  }
	}

	// Irregular faces
	const BaseIVFAB<Real>& diffco = (*m_diffco_eb[lvl])[dit()];
	for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex vof = vofit();
	  const Real thisdt = dx*dx/(2*SpaceDim*diffco(vof, comp));

	  min_dt = Min(min_dt, thisdt);
	}
      }
    }
  
#ifdef CH_MPI
    Real tmp = 1.;
    int result = MPI_Allreduce(&min_dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
    if(result != MPI_SUCCESS){
      MayDay::Error("cdr_solver::compute_diffusive_dt() - communication error on norm");
    }
    min_dt = tmp;
#endif
  }


  return min_dt;
}

Real cdr_solver::compute_source_dt(){
  CH_TIME("cdr_solver::compute_source_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_source_dt" << endl;
  }

  const Real tolerance = 1.E-6;
  const int comp = 0;
  Real max, min;
  data_ops::get_max_min(max, min, m_state, comp);


  Real min_dt = 1.E99;

  if(max > 0.0){
    const int finest_level = m_amr->get_finest_level();
    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
      const Real dx                = m_amr->get_dx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const EBCellFAB& state  = (*m_state[lvl])[dit()];
	const EBCellFAB& source = (*m_source[lvl])[dit()];
	const Box box           = dbl.get(dit());
	const EBISBox& ebisbox  = ebisl[dit()];
	const EBGraph& ebgraph  = ebisbox.getEBGraph();
	const IntVectSet ivs(box);

	for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex vof = vofit();

	  const Real phi = state(vof, comp);
	  const Real src = source(vof, comp);

	  Real thisdt = 1.E99;
	  if(Abs(phi) > tolerance*max){
	    thisdt = Abs(phi/src);
	  }
	  min_dt = Min(min_dt, thisdt);
	}
      }
    }
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&min_dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("cdr_solver::compute_source_dt() - communication error on norm");
  }
  min_dt = tmp;
#endif

  
  return min_dt;

}

Real cdr_solver::compute_mass(){
  CH_TIME("cdr_solver::compute_mass");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_mass" << endl;
  }

  Real mass = 0.;
  const int base = 0;
  const Real dx = m_amr->get_dx()[base];
  m_amr->average_down(m_state, m_phase);
  
  data_ops::kappa_sum(mass, *m_state[base]);
  mass *= pow(dx, SpaceDim);
  
  return mass;
}

bool cdr_solver::is_diffusive(){
  CH_TIME("cdr_solver::is_diffusive");
  if(m_verbosity > 5){
    pout() << m_name + "::is_diffusive" << endl;
  }
  
  return m_diffusive;
}

EBAMRCellData& cdr_solver::get_state(){
  return m_state;
}

EBAMRCellData& cdr_solver::get_source(){
  return m_source;
}

EBAMRCellData& cdr_solver::get_velo_cell(){
  return m_velo_cell;
}

EBAMRFluxData& cdr_solver::get_velo_face(){
  return m_velo_face;
}

EBAMRIVData& cdr_solver::get_velo_eb(){
  return m_velo_eb;
}

EBAMRFluxData& cdr_solver::get_diffco_face(){
  return m_diffco;
}

EBAMRIVData& cdr_solver::get_diffco_eb(){
  return m_diffco_eb;
}

EBAMRIVData& cdr_solver::get_ebflux(){
  return m_ebflux;
}
