/*!
  @file mf_helmholtz_opfactory.cpp
  @brief Implementation of mf_helmholtz_opfactory.H
  @author Robert Marskar
*/

#include "mf_helmholtz_opfactory.H"

int mf_helmholtz_opfactory::s_max_box_size = 32;
int mf_helmholtz_opfactory::s_test_ref     = 32;
int mf_helmholtz_opfactory::s_relax_type   = 2;

mf_helmholtz_opfactory::~mf_helmholtz_opfactory(){

}

void mf_helmholtz_opfactory::set_bottom_drop(const int a_bottom_drop){
  s_test_ref = a_bottom_drop;
}

void mf_helmholtz_opfactory::set_relax_type(int a_relax_type){
  s_relax_type = a_relax_type;
}

void mf_helmholtz_opfactory::set_jump(const EBAMRIVData& a_ajump,
				      const EBAMRIVData& a_bjump,
				      const EBAMRIVData& a_cjump,
				      const EBAMRIVData& a_sigma){
  m_ajump   = a_ajump;
  m_bjump   = a_bjump;
  m_cjump   = a_cjump;
  m_sigma   = a_sigma;
  m_jumpset = true;
}

MGLevelOp<LevelData<MFCellFAB> >* mf_helmholtz_opfactory::MGnewOp(const ProblemDomain& a_fine_ebis, int a_depth, bool a_homo_only){
  MayDay::Abort("mf_helmholtz_opfactory::MGnewOp - not implemented");
  return static_cast<MGLevelOp<LevelData<MFCellFAB> >* > (NULL);
}

AMRLevelOp<LevelData<MFCellFAB> >* mf_helmholtz_opfactory::AMRnewOp(const ProblemDomain& a_fineindexspace){
  MayDay::Abort("mf_helmholtz_opfactory::AMRnewOp - not implemented");
  return static_cast<AMRLevelOp<LevelData<MFCellFAB> >* > (NULL);
}

void mf_helmholtz_opfactory::reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim){
  delete a_reclaim;
}

void mf_helmholtz_opfactory::AMRreclaim(mf_helmholtz_op* a_reclaim){
  delete a_reclaim;
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
