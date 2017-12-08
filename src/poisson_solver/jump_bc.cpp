/*!
  @file jump_bc.cpp
  @brief Implementation of jump_bc.H
  @author Robert Marskar
  @date Dec. 2017
*/

#include "jump_bc.H"

jump_bc::jump_bc(){
  
}

jump_bc::~jump_bc(){

}

void jump_bc::define(const MFLevelGrid& a_mflg, const Real& a_dx, const int a_order){
  CH_TIME("jump_bc::define");
  m_mflg   = a_mflg;
  m_domain = m_mflg.get_domain();
  m_mfis   = m_mflg.get_mfis();
  m_grids  = m_mflg.get_grids();
  m_order  = a_order;

  m_bco.define(m_grids);
  m_soln.define(m_grids);
  m_weights.define(m_grids);
  m_stencils.define(m_grids);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    MFInterfaceFAB<Real>& bco         = m_bco[dit()];
    MFInterfaceFAB<Real>& soln        = m_soln[dit()];
    MFInterfaceFAB<Real>& weights     = m_weights[dit()];
    MFInterfaceFAB<VoFStencil>& stens = m_stencils[dit()];

    bco.define(m_mflg,     dit());
    soln.define(m_mflg,    dit());
    weights.define(m_mflg, dit());
    stens.define(m_mflg,   dit());
  }
}

void jump_bc::build_stencils(){
  CH_TIME("jump_bc::build_stencils");

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    for (int iphase = 0; iphase < m_mfis->num_phases(); iphase ++){

      BaseIVFAB<Real>& bco            = m_bco[dit()].get_ivfab(iphase);
      BaseIVFAB<Real>& weights        = m_weights[dit()].get_ivfab(iphase);
      BaseIVFAB<VoFStencil>& stencils = m_stencils[dit()].get_ivfab(iphase);
    }
  }
}


bool jump_bc::get_second_order_sten(Real& a_weight, VoFStencil& a_stencil, const VolIndex& a_vof, const EBISBox& a_ebisbox){
  CH_TIME("jump_bc::get_second_order_sten");
}

bool jump_bc::get_first_order_sten(Real& a_weight, VoFStencil& a_stencil, const VolIndex& a_vof, const EBISBox& a_ebisbox){
  CH_TIME("jump_bc::get_first_order_sten");
}
