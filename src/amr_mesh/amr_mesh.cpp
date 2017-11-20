/*!
  @file amr_mesh.cpp
  @brief Implementation of amr_mesh.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "amr_mesh.H"

amr_mesh::amr_mesh(){
}

amr_mesh::~amr_mesh(){
}

void amr_mesh::set_mfis(const RefCountedPtr<MFIndexSpace>& a_mfis){
  CH_TIME("amr_mesh::set_mfis");
  if(m_verbosity > 3){
    pout() << "amr_mesh::set_mfis" << endl;
  }

  m_mfis = a_mfis;
}

void amr_mesh::regrid(const Vector<IntVectSet>& a_tags){
  CH_TIME("amr_mesh::regrid");
  if(m_verbosity > 3){
    pout() << "amr_mesh::regrid" << endl;
  }

  this->define_eblevelgrid();
  this->define_eb_coar_ave();
  this->define_eb_quad_cfi();
  this->define_redist_oper();
  this->define_flux_reg();
}

void amr_mesh::average_down(EBAMRCellData& a_data, Phase::WhichPhase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(cell)" << endl;
  }

  //
  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);

    a_data[lvl]->exchange(interv);
  }
}

void amr_mesh::average_down(EBAMRFluxData& a_data, Phase::WhichPhase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(face)" << endl;
  }

  //
  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);

    a_data[lvl]->exchange(interv);
  }
}

void amr_mesh::average_down(EBAMRIVData& a_data, Phase::WhichPhase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(iv)" << endl;
  }

  //
  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);

    a_data[lvl]->exchange(interv);
  }
}

void amr_mesh::interp_ghost(EBAMRCellData& a_data, Phase::WhichPhase a_phase){
  CH_TIME("amr_mesh::interp_ghost");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost" << endl;
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv(0, ncomps -1);

    m_quadcfi[a_phase][lvl]->coarseFineInterp(*a_data[lvl], *a_data[lvl-1], 0, 0, ncomps);

    a_data[lvl]->exchange(interv);
  }
}

void amr_mesh::set_verbosity(const int a_verbosity){
  CH_TIME("amr_mesh::set_verbosity");

  m_verbosity = a_verbosity;
}

void amr_mesh::set_max_box_size(const int a_max_box_size){
  CH_TIME("amr_mesh::set_max_box_size");
  CH_assert(a_max_box_size >= 4);
  
  m_max_box_size = a_max_box_size;
}

void amr_mesh::set_blocking_factor(const int a_blocking_factor){
  CH_TIME("amr_mesh::set_blocking_factor");
  CH_assert(a_blocking_factor >= 4);
  
  m_blocking_factor = a_blocking_factor;
}

void amr_mesh::set_eb_ghost(const int a_ebghost){
  CH_TIME("amr_mesh:.set_eb_ghost");
  CH_assert(a_ebghost > 0);
  m_ebghost = a_ebghost;
}

void amr_mesh::set_num_ghost(const int a_num_ghost){
  CH_TIME("amr_mesh::set_num_ghost");
  CH_assert(a_num_ghost > 0);
  m_num_ghost = a_num_ghost;
}

void amr_mesh::set_redist_rad(const int a_redist_rad){
  CH_TIME("amr_mesh::set_redist_rads");
  CH_assert(a_redist_rad > 0);
  m_redist_rad = a_redist_rad;
}

int amr_mesh::get_finest_level(){
  return m_finest_level;
}

int amr_mesh::get_max_amr_depth(){
  return m_max_amr_depth;
}

Vector<Real>& amr_mesh::get_dx(){
  return m_dx;
}

Vector<int>& amr_mesh::get_ref_rat(){
  return m_refinement_ratios;
}

Vector<DisjointBoxLayout>& amr_mesh::get_grids(){
  return m_grids;
}

Vector<ProblemDomain>& amr_mesh::get_domains(){
  return m_domains;
}

Vector<RefCountedPtr<EBLevelGrid> >& amr_mesh::get_eblg(Phase::WhichPhase a_phase){
  return m_eblg[a_phase];
}

Vector<EBISLayout>& amr_mesh::get_ebisl(Phase::WhichPhase a_phase){
  return m_ebisl[a_phase];
}

Vector<RefCountedPtr<EBCoarseAverage> >& amr_mesh::get_coarave(Phase::WhichPhase a_phase){
  return m_coarave[a_phase];
}

Vector<RefCountedPtr<nwoebquadcfinterp> >& amr_mesh::get_quadcfi(Phase::WhichPhase a_phase){
  return m_quadcfi[a_phase];
}

Vector<RefCountedPtr<EBFastFR> >&  amr_mesh::get_flux_reg(Phase::WhichPhase a_phase){
  CH_assert(a_phase == Phase::Gas); // This is disabled since we only solve cdr in gas phase. 
  return m_flux_reg[a_phase];
}

Vector<RefCountedPtr<EBLevelRedist> >& amr_mesh::get_level_redist(Phase::WhichPhase a_phase){
  CH_assert(a_phase == Phase::Gas); // This is disabled since we only solve cdr in gas phase. 
  return m_level_redist[a_phase];
}

Vector<RefCountedPtr<EBCoarToFineRedist> >&  amr_mesh::get_coar_to_fine_redist(Phase::WhichPhase a_phase){
  CH_assert(a_phase == Phase::Gas); // This is disabled since we only solve cdr in gas phase.

  return m_coar_to_fine_redist[a_phase];;
}

Vector<RefCountedPtr<EBCoarToCoarRedist> >&  amr_mesh::get_coar_to_coar_redist(Phase::WhichPhase a_phase){
  CH_assert(a_phase == Phase::Gas); // This is disabled since we only solve cdr in gas phase.

  return m_coar_to_coar_redist[a_phase];;
}

Vector<RefCountedPtr<EBFineToCoarRedist> >&  amr_mesh::get_fine_to_coar_redist(Phase::WhichPhase a_phase){
  CH_assert(a_phase == Phase::Gas); // This is disabled since we only solve cdr in gas phase.

  return m_fine_to_coar_redist[a_phase];;
}

void amr_mesh::define_eblevelgrid(){
  CH_TIME("amr_mesh::define_eblevelgrid");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_levelgrid" << endl;
  }

  const EBIndexSpace* const ebis_gas = m_mfis->EBIS(Phase::Gas);
  const EBIndexSpace* const ebis_sol = m_mfis->EBIS(Phase::Solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    if(ebis_gas != NULL){
      m_eblg[Phase::Gas][lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl], m_domains[lvl], m_ebghost, ebis_gas));
      m_ebisl[Phase::Gas][lvl] = m_eblg[Phase::Gas][lvl]->getEBISL();
    }
    if(ebis_sol != NULL){
      m_eblg[Phase::Solid][lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl], m_domains[lvl], m_ebghost, ebis_sol));
      m_ebisl[Phase::Solid][lvl] = m_eblg[Phase::Solid][lvl]->getEBISL();
    }
  }
}

void amr_mesh::define_eb_coar_ave(){
  CH_TIME("amr_mesh::define_eb_coar_ave");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_eb_coar_ave" << endl;
  }

  const int comps = SpaceDim;

  const EBIndexSpace* const ebis_gas = m_mfis->EBIS(Phase::Gas);
  const EBIndexSpace* const ebis_sol = m_mfis->EBIS(Phase::Solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(ebis_gas != NULL){
	m_coarave[Phase::Gas][lvl] = RefCountedPtr<EBCoarseAverage> (new EBCoarseAverage(m_grids[lvl],
											 m_grids[lvl-1],
											 m_ebisl[Phase::Gas][lvl],
											 m_ebisl[Phase::Gas][lvl-1],
											 m_domains[lvl-1],
											 m_refinement_ratios[lvl-1],
											 comps,
											 ebis_gas));
      }
      if(ebis_sol != NULL){
	m_coarave[Phase::Solid][lvl] = RefCountedPtr<EBCoarseAverage> (new EBCoarseAverage(m_grids[lvl],
											   m_grids[lvl-1],
											   m_ebisl[Phase::Solid][lvl],
											   m_ebisl[Phase::Solid][lvl-1],
											   m_domains[lvl-1],
											   m_refinement_ratios[lvl-1],
											   comps,
											   ebis_sol));
      }
    }
  }
}

void amr_mesh::define_eb_quad_cfi(){
  CH_TIME("amr_mesh::define_eb_quad_cfi");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_eb_quad_cfi" << endl;
  }

  const int comps = SpaceDim;

  const EBIndexSpace* const ebis_gas = m_mfis->EBIS(Phase::Gas);
  const EBIndexSpace* const ebis_sol = m_mfis->EBIS(Phase::Solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(ebis_gas != NULL){
	LayoutData<IntVectSet>& cfivs = *(m_eblg[Phase::Gas][lvl]->getCFIVS());
	m_quadcfi[Phase::Gas][lvl] = RefCountedPtr<nwoebquadcfinterp> (new nwoebquadcfinterp(m_grids[lvl],
											     m_grids[lvl-1],
											     m_ebisl[Phase::Gas][lvl],
											     m_ebisl[Phase::Gas][lvl-1],
											     m_domains[lvl-1],
											     m_refinement_ratios[lvl-1],
											     comps,
											     m_dx[lvl],
											     m_num_ghost*IntVect::Unit,
											     cfivs,
											     ebis_gas));
      }
      if(ebis_sol != NULL){
	LayoutData<IntVectSet>& cfivs = *(m_eblg[Phase::Solid][lvl]->getCFIVS());
	m_quadcfi[Phase::Solid][lvl] = RefCountedPtr<nwoebquadcfinterp> (new nwoebquadcfinterp(m_grids[lvl],
											       m_grids[lvl-1],
											       m_ebisl[Phase::Solid][lvl],
											       m_ebisl[Phase::Solid][lvl-1],
											       m_domains[lvl-1],
											       m_refinement_ratios[lvl-1],
											       comps,
											       m_dx[lvl],
											       m_num_ghost*IntVect::Unit,
											       cfivs,
											       ebis_sol));
      }
    }
  }
}

void amr_mesh::define_flux_reg(){
  CH_TIME("amr_mesh::define_flux_reg");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_flux_reg" << endl;
  }

  const int comps = 1;

  const EBIndexSpace* const ebis_gas = m_mfis->EBIS(Phase::Gas);
  const EBIndexSpace* const ebis_sol = m_mfis->EBIS(Phase::Solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const bool has_fine = lvl < m_finest_level;

    if(has_fine){
      if(ebis_gas != NULL){
	m_flux_reg[Phase::Solid][lvl] = RefCountedPtr<EBFastFR> (new EBFastFR(*m_eblg[Phase::Gas][lvl+1],
									      *m_eblg[Phase::Gas][lvl],
									      m_refinement_ratios[lvl],
									      comps,
									      !m_ebcf));
      }
    }
  }
}

void amr_mesh::define_redist_oper(){
  CH_TIME("amr_mesh::define_redist_oper");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_redist_oper" << endl;
  }

  const int comps = 1;

  const EBIndexSpace* const ebis_gas = m_mfis->EBIS(Phase::Gas);
  const EBIndexSpace* const ebis_sol = m_mfis->EBIS(Phase::Solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl > m_finest_level;

    m_level_redist[Phase::Gas][lvl] = RefCountedPtr<EBLevelRedist> (new EBLevelRedist(m_grids[lvl],
										      m_ebisl[Phase::Gas][lvl],
										      m_domains[lvl],
										      comps,
										      m_redist_rad));

    
    if(m_ebcf){
      if(has_coar){
	m_fine_to_coar_redist[Phase::Gas][lvl] = RefCountedPtr<EBFineToCoarRedist> (new EBFineToCoarRedist());
	m_fine_to_coar_redist[Phase::Gas][lvl]->define(m_grids[lvl],
						       m_grids[lvl-1],
						       m_ebisl[Phase::Gas][lvl],
						       m_ebisl[Phase::Gas][lvl-1],
						       m_domains[lvl-1].domainBox(),
						       m_refinement_ratios[lvl-1],
						       comps,
						       m_redist_rad,
						       ebis_gas);

	// Set register to zero
	m_fine_to_coar_redist[Phase::Gas][lvl]->setToZero();
      }

      if(has_fine){

	m_coar_to_fine_redist[Phase::Gas][lvl] = RefCountedPtr<EBCoarToFineRedist> (new EBCoarToFineRedist());
	m_coar_to_fine_redist[Phase::Gas][lvl]->define(m_grids[lvl+1],
						       m_grids[lvl],
						       m_ebisl[Phase::Gas][lvl],
						       m_domains[lvl].domainBox(),
						       m_refinement_ratios[lvl],
						       comps,
						       m_redist_rad,
						       ebis_gas);

	// Coarse to coarse redistribution
	m_coar_to_coar_redist[Phase::Gas][lvl] = RefCountedPtr<EBCoarToCoarRedist> (new EBCoarToCoarRedist());
	m_coar_to_coar_redist[Phase::Gas][lvl]->define(*m_eblg[Phase::Gas][lvl+1],
						       *m_eblg[Phase::Gas][lvl],
						       m_refinement_ratios[lvl],
						       comps,
						       m_redist_rad);


	// Set registers to zero
	m_coar_to_fine_redist[Phase::Gas][lvl]->setToZero();
	m_coar_to_coar_redist[Phase::Gas][lvl]->setToZero();
      }
    }
  }
}


