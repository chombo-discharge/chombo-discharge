/*!
  @file mf_helmholtz_opfactory.cpp
  @brief Implementation of mf_helmholtz_opfactory.H
  @author Robert Marskar
*/

#include "mf_helmholtz_opfactory.H"

int mf_helmholtz_opfactory::s_drop_bottom = 2;
int mf_helmholtz_opfactory::s_relax_type  = 2;


mf_helmholtz_opfactory::~mf_helmholtz_opfactory(){

}

void mf_helmholtz_opfactory::set_bottom_drop(const int a_bottom_drop){
  s_drop_bottom = a_bottom_drop;
}

void mf_helmholtz_opfactory::set_relax_type(int a_relax_type){
  s_relax_type = a_relax_type;
}

void mf_helmholtz_opfactory::set_jump(const Real a_gas_perm, const EBAMRIVData& a_solid_perm, const EBAMRIVData& a_jump){
  MayDay::Abort("mf_helmholtz_opfactory::not implemented");
}

MGLevelOp<LevelData<MFCellFAB> >* mf_helmholtz_opfactory::MGnewOp(const ProblemDomain& a_fine_ebis, int a_depth, bool a_homo_only){
  MayDay::Abort("mf_helmholtz_opfactory::MGnewOp - not implemented");
  return static_cast<MGLevelOp<LevelData<MFCellFAB> >* > (NULL);
}

AMRLevelOp<LevelData<MFCellFAB> >* mf_helmholtz_opfactory::AMRnewOp(const ProblemDomain& a_fineindexspace){
  MayDay::Abort("mf_helmholtz_opfactory::AMRnewOp - not implemented");
  return static_cast<AMRLevelOp<LevelData<MFCellFAB> >* > (NULL);
}

mf_helmholtz_op* mf_helmholtz_opfactory::createOperator(const DisjointBoxLayout&       a_dilboMGLevel,
							const DisjointBoxLayout&       a_dilboCoarMG,
							const ProblemDomain&           a_domainMGLevel,
							const bool&                    a_hasMGObjects,
							const bool&                    a_layoutChanged,
							const RealVect&                a_dxMGLevel,
							const RealVect&                a_dxCoar,
							const int&                     a_whichLevel,
							const int&                     a_mgLevel){
  MayDay::Abort("mf_helmholtz_opfactory::createOperator - not implemented");
  return static_cast<mf_helmholtz_op*> (NULL);
}

void mf_helmholtz_opfactory::reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim){
  MayDay::Abort("mf_helmholtz_opfactory::reclaim - not implemented");
}

void mf_helmholtz_opfactory::AMRreclaim(mf_helmholtz_op* a_reclaim){
  MayDay::Abort("mf_helmholtz_opfactory::AMRreclaim - not implemented");
}

int mf_helmholtz_opfactory::refToFiner(const ProblemDomain& a_domain) const{
  MayDay::Abort("mf_helmholtz_opfactory::refToFiner - not implemented");
}
