/*!
  @file   sigma_solver.cpp
  @brief  Implementation of sigma_solver.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "sigma_solver.H"
#include "data_ops.H"

#include <EBArith.H>

sigma_solver::sigma_solver(){

  this->set_verbosity(-1);
  this->set_phase(phase::gas);
  this->set_plot_variables();
}

sigma_solver::~sigma_solver(){

}

void sigma_solver::advance(const Real a_dt){
  CH_TIME("sigma_solver::advance(dt)");
  if(m_verbosity > 5){
    pout() << "sigma_solver::advance(dt)" << endl;
  }

  this->advance(m_state, a_dt);
}

void sigma_solver::advance(EBAMRIVData& a_state, const Real a_dt){
  CH_TIME("sigma_solver::advance(state, dt)");
  if(m_verbosity > 5){
    pout() << "sigma_solver::advance(state, dt)" << endl;
  }

  const Real midpoint = 0.5;

  this->advance_rk2(a_state, a_dt, midpoint);
}

void sigma_solver::advance_rk2(EBAMRIVData& a_state, const Real a_dt, const Real a_alpha){
  CH_TIME("sigma_solver::advance_rk2");
  if(m_verbosity > 5){
    pout() << "sigma_solver::advance_rk2" << endl;
  }

  MayDay::Abort("sigma_solver::advance_rk2 - not implemented (yet)");
}

void sigma_solver::allocate_internals(){
  CH_TIME("sigma_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << "sigma_solver::allocate_internals" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  
  m_amr->allocate(m_state, m_phase, ncomp);
  m_amr->allocate(m_flux,  m_phase, ncomp);
}

void sigma_solver::cache_state(){
  CH_TIME("sigma_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << "sigma_solver::allocate_internals" << endl;
  }

  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  
  m_amr->allocate(m_cache, m_phase, ncomp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    m_state[lvl]->localCopyTo(*m_cache[lvl]);
  }
}

void sigma_solver::compute_rhs(EBAMRIVData& a_rhs){
  CH_TIME("sigma_solver::compute_rhs");
  if(m_verbosity > 5){
    pout() << "sigma_solver::compute_rhs" << endl;
  }
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl ++){
    m_flux[lvl]->localCopyTo(*a_rhs[lvl]);
  }
}

void sigma_solver::deallocate_internals(){
  CH_TIME("sigma_solver::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "sigma_solver::deallocate_internals" << endl;
  }
  
  m_amr->deallocate(m_state);
  m_amr->deallocate(m_flux);
}

void sigma_solver::initial_data(){
  CH_TIME("sigma_solver::initial_data");
  if(m_verbosity > 5){
    pout() << "sigma_solver::initial_data" << endl;
  }

  const RealVect origin  = m_amr->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    
    for (DataIterator dit = m_state[lvl]->dataIterator(); dit.ok(); ++dit){
      BaseIVFAB<Real>& state = (*m_state[lvl])[dit()];

      const EBISBox& ebisbox = ebisl[dit()];
      const IntVectSet& ivs  = state.getIVS();
      const EBGraph& ebgraph = state.getEBGraph();
      
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = origin + vof.gridIndex()*dx + ebisbox.bndryCentroid(vof)*dx;
	
	for (int comp = 0; comp < state.nComp(); comp++){
	  state(vof, comp) = m_plaskin->initial_sigma(m_time, pos);
	}
      }
    }
  }

  m_amr->average_down(m_state, m_phase);
  this->reset_cells(m_state);
}

void sigma_solver::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("sigma_solver::regrid");
  if(m_verbosity > 5){
    pout() << "sigma_solver::regrid" << endl;
  }

  const RefCountedPtr<EBIndexSpace> ebis = m_mfis->get_ebis(m_phase);

  const int ebghost = 4; // m_amr->get_eb_ghost();
  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);

  this->allocate_internals();

  // These levels have never changed
  for (int lvl = 0; lvl <= Max(0, a_lmin-1); lvl++){
    m_cache[lvl]->copyTo(*m_state[lvl]); 
  }

  // These levels have changed
  for (int lvl = a_lmin; lvl <= a_new_finest_level; lvl++){
    const DisjointBoxLayout& fine_grid = m_amr->get_grids()[lvl];
    const ProblemDomain& fine_domain   = m_amr->get_domains()[lvl];
    const ProblemDomain& coar_domain   = m_amr->get_domains()[lvl-1];
    const EBISLayout& fine_ebisl       = m_amr->get_ebisl(m_phase)[lvl];
    const int nref                     = m_amr->get_ref_rat()[lvl-1];
    
    // Fill a coarsened grid and a coarsened ebisl
    EBISLayout coar_ebisl;
    DisjointBoxLayout coar_grid = DisjointBoxLayout();
    coarsen(coar_grid, fine_grid, nref);
    ebis->fillEBISLayout(coar_ebisl, coar_grid, coar_domain, ebghost);

    // Need extra storage for coarse data
    LayoutData<IntVectSet> sets(coar_grid);
    for (DataIterator dit = coar_grid.dataIterator(); dit.ok(); ++dit){
      Box box          = coar_grid.get(dit());
      box.grow(3);
      box &= coar_domain;
      
      const EBISBox& ebisbox = coar_ebisl[dit()];
      const IntVectSet ivs   = ebisbox.getIrregIVS(box);

      sets[dit()] = ivs;
    }

    // Allocate storage for coarsened data
    BaseIVFactory<Real> ivfact(coar_ebisl, sets);
    LevelData<BaseIVFAB<Real> > coarsened_fine_data(coar_grid, ncomp, m_state[0]->ghostVect(), ivfact);

    //
    m_state[lvl-1]->copyTo(coarsened_fine_data);


    // Loop through coarse grid and interpolate to fine grid
    for (DataIterator dit = coar_grid.dataIterator(); dit.ok(); ++dit){
      BaseIVFAB<Real>& fine_state       = (*m_state[lvl])[dit()];
      const BaseIVFAB<Real>& coar_state = coarsened_fine_data[dit()];
      const EBISBox& coar_ebisbox       = coar_ebisl[dit()];
      const EBISBox& fine_ebisbox       = fine_ebisl[dit()];

      const IntVectSet ivs   = fine_state.getIVS();
      const EBGraph& ebgraph = fine_state.getEBGraph();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
    	const VolIndex& fine_vof = vofit();
    	const VolIndex  coar_vof = fine_ebisl.coarsen(fine_vof, nref, dit());

    	// Get all the area fractions of the refine coar vof in order to find weighting fractions
    	Vector<VolIndex> refined_vofs = coar_ebisl.refine(coar_vof, nref, dit());
    	Real sumfrac = 0.0;
    	for (int i = 0; i < refined_vofs.size(); i++){
	  sumfrac += fine_ebisbox.bndryArea(refined_vofs[i]);
    	}

    	// Initialize conserved charge
    	const Real coarFrac = coar_ebisbox.bndryArea(coar_vof);
	fine_state(fine_vof, comp) = coar_state(coar_vof, comp)*nref/sumfrac;
      }
    }

    // If data already exists, it takes precedence
    if (lvl <= a_old_finest_level){
      m_cache[lvl]->copyTo(*m_state[lvl]);
    }
  }

  this->reset_cells(m_state);
}

void sigma_solver::reset_cells(EBAMRIVData& a_data){
  CH_TIME("sigma_solver::reset_cells");
  if(m_verbosity > 5){
    pout() << "sigma_solver::reset_cells" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    const MFLevelGrid& mflg      = *m_amr->get_mflg()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      BaseIVFAB<Real>& data  = (*a_data[lvl])[dit()];
      IntVectSet ivs         = data.getIVS();
      const EBGraph& ebgraph = data.getEBGraph();
      
      ivs -= mflg.interface_region(box, dit());

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	for (int comp = 0; comp < data.nComp(); comp++){
	  data(vof, comp) = 0.0;
	}
      }
    }
  }
}

void sigma_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("sigma_solver::set_amr");
  if(m_verbosity > 5){
    pout() << "sigma_solver::set_amr" << endl;
  }

  m_amr = a_amr;
}

void sigma_solver::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  CH_TIME("sigma_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << "sigma_solver::set_computational_geometry" << endl;
  }

  m_compgeom = a_compgeom;
  m_mfis     = m_compgeom->get_mfis();
}

void sigma_solver::set_plasma_kinetics(const RefCountedPtr<plasma_kinetics>& a_plaskin){
  CH_TIME("sigma_solver::set_plasma_kinetics");
  if(m_verbosity > 5){
    pout() << "sigma_solver::set_plasma_kinetics" << endl;
  }

  m_plaskin = a_plaskin;
}

void sigma_solver::set_phase(phase::which_phase a_phase){
  CH_TIME("sigma_solver::set_phase");
  if(m_verbosity > 5){
    pout() << "sigma_solver::set_phase" << endl;
  }

  m_phase = a_phase;
}

void sigma_solver::set_sigma(const EBAMRIVData& a_sigma){
  CH_TIME("sigma_solver::set_sigma(ebamrivdata)");
  if(m_verbosity > 5){
    pout() << "sigma_solver::set_sigma(ebamrivdata)" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_sigma[lvl]->localCopyTo(*m_state[lvl]);
  }

  this->reset_cells(m_state);
}

void sigma_solver::set_sigma(const Real a_sigma){
  CH_TIME("sigma_solver::set_sigma(constant)");
  if(m_verbosity > 5){
    pout() << "sigma_solver::set_sigma(constant)" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::set_value(*m_state[lvl], a_sigma);
  }

  this->reset_cells(m_state);
}

void sigma_solver::set_verbosity(const int a_verbosity){
  CH_TIME("sigma_solver::set_verbosity");
  m_verbosity = a_verbosity;
  
  if(m_verbosity > 5){
    pout() << "sigma_solver::set_verbosity" << endl;
  }
  

}

void sigma_solver::set_time(const int a_step, const Real a_time, const Real a_dt){
  CH_TIME("sigma_solver::set_time");
  if(m_verbosity > 5){
    pout() << "sigma_solver::set_time" << endl;
  }
  
  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void sigma_solver::write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("sigma_solver::write_checkpoint_level");
  if(m_verbosity > 5){
    pout() << "sigma_solver::write_checkpoint_level" << endl;
  }

  EBCellFactory fact(m_amr->get_ebisl(phase::gas)[a_level]);
  LevelData<EBCellFAB> scratch(m_amr->get_grids()[a_level], 1, 3*IntVect::Unit, fact);
  data_ops::set_value(scratch, 0.0);
  data_ops::incr(scratch, *m_state[a_level], 1.0);

  // Write state vector
  write(a_handle, scratch, "sigma");
}

void sigma_solver::read_checkpoint_level(HDF5Handle& a_handle, const int a_level){
  CH_TIME("sigma_solver::read_checkpoint_level");
  if(m_verbosity > 5){
    pout() << "sigma_solver::read_checkpoint_level" << endl;
  }

  const EBISLayout& ebisl = m_amr->get_ebisl(phase::gas)[a_level];
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_level];
  
  EBCellFactory fact(ebisl);
  LevelData<EBCellFAB> scratch(dbl, 1, 3*IntVect::Unit, fact);
  data_ops::set_value(scratch, 0.0);
  read<EBCellFAB>(a_handle, scratch, "sigma", dbl, Interval(0,0), false);

		     
  data_ops::set_value(*m_state[a_level], 0.0);
  data_ops::incr(*m_state[a_level], scratch, 1.0);
}

void sigma_solver::write_plot_data(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("sigma_solver::write_plot_data");
  if(m_verbosity > 5){
    pout() << "sigma_solver::write_plot_data" << endl;
  }


  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_phase, 1);


  // Write sigma
  data_ops::set_value(scratch, 0.0);
  data_ops::incr(scratch, m_state, 1.0);
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    scratch[lvl]->localCopyTo(Interval(0,0), *a_output[lvl], Interval(a_comp, a_comp));
  }
  a_comp++;

  // Write flux
  data_ops::set_value(scratch, 0.0);
  data_ops::incr(scratch, m_flux, 1.0);
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    scratch[lvl]->localCopyTo(Interval(0,0), *a_output[lvl], Interval(a_comp, a_comp));
  }
  a_comp++;
}

void sigma_solver::set_plot_variables(){
  CH_TIME("sigma_solver::set_plot_variables");
  if(m_verbosity > 5){
    pout() << "sigma_solver::set_plot_variables" << endl;
  }

  m_plot_phi = true;
}

int sigma_solver::get_num_plotvars(){
  CH_TIME("sigma_solver::get_num_plotvars");
  if(m_verbosity > 5){
    pout() << "sigma_solver::get_num_plotvars" << endl;
  }
  
  return 2;
}
  
  
Vector<std::string> sigma_solver::get_plotvar_names() const{
  CH_TIME("sigma_solver::get_num_plotvars");
  if(m_verbosity > 5){
    pout() << "sigma_solver::get_num_plotvars" << endl;
  }
  Vector<std::string> ret(2);

  ret[0] = "surface charge density";
  ret[1] = "surface charge flux";
  
  return ret;
}


Real sigma_solver::compute_charge(){
  CH_TIME("sigma_solver::compute_charge");
  if(m_verbosity > 5){
    pout() << "sigma_solver::compute_charge" << endl;
  }

  m_amr->average_down(m_state, m_phase);

  Real charge = 0.0;

  const int comp               = 0;
  const DisjointBoxLayout& dbl = m_amr->get_grids()[0];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[0];
  const Real dx                = m_amr->get_dx()[0];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const Box box          = dbl.get(dit());
    const IntVectSet irreg = ebisbox.getIrregIVS(box);

    for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real& area    = ebisbox.bndryArea(vof);

      charge += area*(*m_state[0])[dit()](vof, comp);
    }
  }

  data_ops::sum(charge); // Parallell sum
  
  return charge*dx;
}

EBAMRIVData& sigma_solver::get_state(){
  return m_state;
}

EBAMRIVData& sigma_solver::get_flux(){
  return m_flux;
}
