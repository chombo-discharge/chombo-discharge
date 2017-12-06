/*!
  @file mf_helmholtz_opfactory.cpp
  @brief Implementation of mf_helmholtz_opfactory.H
  @author Robert Marskar
*/

#include "mf_helmholtz_opfactory.H"

int mf_helmholtz_opfactory::s_max_box_size = 32;
int mf_helmholtz_opfactory::s_test_ref     = 4;
int mf_helmholtz_opfactory::s_relax_type   = 2;

mf_helmholtz_opfactory::~mf_helmholtz_opfactory(){

}

void mf_helmholtz_opfactory::set_bottom_drop(const int a_bottom_drop){
  s_test_ref = a_bottom_drop;
}

void mf_helmholtz_opfactory::set_relax_type(int a_relax_type){
  s_relax_type = a_relax_type;
}

void mf_helmholtz_opfactory::reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim){
  delete a_reclaim;
}

void mf_helmholtz_opfactory::AMRreclaim(mf_helmholtz_op* a_reclaim){
  delete a_reclaim;
}

void mf_helmholtz_opfactory::average_down_amr(){
  CH_TIME("mf_helmholtz_opfactory::average_down_amr");
    
  const int ncomp        = 0;
  const Interval interv  = Interval(0, ncomp -1);
  const int finest_level = m_num_levels - 1;

  for (int lvl = finest_level; lvl > 0; lvl--){ // Average down AMR levels
    m_aveop[lvl]->average(*m_jump[lvl-1], *m_jump[lvl], interv);
  }
}

void mf_helmholtz_opfactory::set_jump(const Real& a_sigma, const Real& a_scale){
  CH_TIME("mf_helmholtz_opfactory::set_jump(scalar)");
  for (int lvl = 0; lvl < m_num_levels; lvl++){
    EBLevelDataOps::setVal(*m_jump[lvl], a_sigma);
    data_ops::scale(*m_jump[lvl], a_sigma);
  }

  this->average_down_amr();
}

void mf_helmholtz_opfactory::set_jump(const EBAMRIVData& a_sigma, const Real& a_scale){
  CH_TIME("mf_helmholtz_opfactory::set_jump(data based)");

  for (int lvl = 0; lvl < m_num_levels; lvl++){
    a_sigma[lvl]->copyTo(*m_jump[lvl]);
    data_ops::scale(*m_jump[lvl], a_scale);
  }

  this->average_down_amr();
}

int mf_helmholtz_opfactory::refToFiner(const ProblemDomain& a_domain) const{
  int retval = -1;
  bool found = false;

  for (int lvl = 0; lvl < m_domains.size(); lvl++){
    if(m_domains[lvl] == a_domain){
      retval = m_ref_rat[lvl];
      found  = true;
    }
  }
  
  if(!found){
    MayDay::Error("mf_helmholtz_opfactory::refToFiner - domain not found in AMR hierarchy");
  }
  
  return retval;
}

AMRLevelOp<LevelData<MFCellFAB> >* mf_helmholtz_opfactory::AMRnewOp(const ProblemDomain& a_fineindexspace){
  MayDay::Abort("mf_helmholtz_opfactory::AMRnewOp - not implemented");
  return static_cast<AMRLevelOp<LevelData<MFCellFAB> >* > (NULL);
}
