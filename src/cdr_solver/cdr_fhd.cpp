/*!
  @file   cdr_fhd.cpp
  @brief  Implementation of cdr_fhd.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "cdr_fhd.H"
#include "cdr_fhdF_F.H"

#include <ParmParse.H>

cdr_fhd::cdr_fhd() : cdr_gdnv() {
  m_name     = "cdr_fhd";
  
  std::string str;
  ParmParse pp("cdr_fhd");
  
  pp.get("stochastic_diffusion", str); m_stochastic_diffusion = (str == "true") ? true : false;
  pp.get("stochastic_advection", str); m_stochastic_advection = (str == "true") ? true : false;
  pp.get("limit_slopes", str);         m_slopelim             = (str == "true") ? true : false;

  this->set_divF_nc(1); // Use covered face stuff for nonconservative divergence
}

cdr_fhd::~cdr_fhd(){
  this->delete_covered();
}


void cdr_fhd::advance_euler(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt){
  CH_TIME("cdr_fhd::advance_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_euler (no source)" << endl;
  }
  
  if(!m_stochastic_diffusion){
    cdr_tga::advance_euler(a_new_state, a_old_state, a_dt);
  }
  else{
    MayDay::Abort("cdr_fhd::advance_euler - stochastic not implemented");
  }
}

void cdr_fhd::advance_euler(EBAMRCellData&       a_new_state,
			    const EBAMRCellData& a_old_state,
			    const EBAMRCellData& a_source,
			    const Real           a_dt){
  CH_TIME("cdr_fhd::advance_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_euler (with source)" << endl;
  }
  
  if(!m_stochastic_diffusion){
    cdr_tga::advance_euler(a_new_state, a_old_state, a_source, a_dt);
  }
  else{


    // Face states subject to random flux
    EBAMRFluxData face_states;
    m_amr->allocate(face_states, m_phase, 1);
    this->smooth_heaviside_faces(face_states, a_old_state);

    
    MayDay::Abort("cdr_fhd::advance_euler2 - stochastic not implemented");
  }
}

void cdr_fhd::advance_tga(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt){
  CH_TIME("cdr_fhd::advance_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_tga (no source)" << endl;
  }
  
  if(!m_stochastic_diffusion){
    cdr_tga::advance_tga(a_new_state, a_old_state, a_dt);
  }
  else{
    MayDay::Abort("cdr_fhd::advance_tga - stochastic not implemented");
  }
}

void cdr_fhd::advance_tga(EBAMRCellData&       a_new_state,
			  const EBAMRCellData& a_old_state,
			  const EBAMRCellData& a_source,
			  const Real           a_dt){
  CH_TIME("cdr_fhd::advance_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_tga (with source)" << endl;
  }
  
  if(!m_stochastic_diffusion){
    cdr_tga::advance_tga(a_new_state, a_old_state, a_source, a_dt);
  }
  else{
    MayDay::Abort("cdr_fhd::advance_tga2 - stochastic not implemented");
  }
}

void cdr_fhd::smooth_heaviside_faces(EBAMRFluxData& a_face_states, const EBAMRCellData& a_cell_states){
  CH_TIME("cdr_fhd::smooth_heaviside_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::smooth_heaviside_faces" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();

      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& face_states       = (*a_face_states[lvl])[dit()][dir];
	const EBCellFAB& cell_states = (*a_cell_states[lvl])[dit()];

	BaseFab<Real>& reg_face       = face_states.getSingleValuedFAB();
	const BaseFab<Real>& reg_cell = face_states.getSingleValuedFAB();

	Box facebox = box;
	facebox.surroundingNodes(dir);
	
	FORT_HEAVISIDE_MEAN(CHF_FRA1(reg_face, comp),
			    CHF_CONST_FRA1(reg_cell, comp),
			    CHF_BOX(facebox));
      }
    }

    // Covered is bogus
    EBLevelDataOps::setCoveredVal(*a_face_states[lvl], 0.0);
  }
}
