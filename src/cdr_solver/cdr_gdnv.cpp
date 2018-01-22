/*!
  @file cdr_gdnv.cpp
  @brief Implementation of cdr_gdnv.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "cdr_gdnv.H"

#include <ExtrapAdvectBC.H>

cdr_gdnv::cdr_gdnv() : cdr_solver() {
  m_name = "cdr_gdnv";
}


cdr_gdnv::~cdr_gdnv(){

}

int cdr_gdnv::query_ghost() const {
  return 3;
}
  

void cdr_gdnv::extrapolate_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state){
  CH_TIME("cdr_gdnv::extrapolate_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::extrapolate_to_faces" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  const int comp = 0;
  
  EBAdvectPatchIntegrator::setCurComp(0);
  EBAdvectPatchIntegrator::setDoingVel(0);

  RefCountedPtr<ExtrapAdvectBCFactory> bcfact = RefCountedPtr<ExtrapAdvectBCFactory>
    (new ExtrapAdvectBCFactory());

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const RefCountedPtr<EBAdvectLevelIntegrator>& leveladvect = m_amr->get_level_advect(m_phase)[lvl];
    leveladvect->resetBCs(bcfact);
    CH_assert(!leveladvect.isNull());

    const bool has_coar = lvl > 0;
    LevelData<EBCellFAB>* coarstate_old = NULL;
    LevelData<EBCellFAB>* coarstate_new = NULL;
    LevelData<EBCellFAB>* coarvel_old   = NULL;
    LevelData<EBCellFAB>* coarvel_new   = NULL;
    LevelData<EBCellFAB>* coarsrc_old   = NULL;
    LevelData<EBCellFAB>* coarsrc_new   = NULL;

    
    if(has_coar){
      coarstate_old = a_state[lvl-1];
      coarstate_new = a_state[lvl-1];
      coarvel_old   = m_velo_cell[lvl-1];
      coarvel_new   = m_velo_cell[lvl-1];
      coarsrc_old   = m_source[lvl-1];
      coarsrc_new   = m_source[lvl-1];
    }
    leveladvect->advectToFacesBCG(*a_face_state[lvl],
				  *a_state[lvl],
				  *m_velo_cell[lvl],
				  *m_velo_face[lvl],
				  coarstate_old,
				  coarstate_new,
				  coarvel_old,
				  coarvel_new,
				  m_time,
				  m_time,
				  m_time,
				  m_dt,
				  m_source[lvl], 
				  coarsrc_old,
				  coarsrc_new);
  }
}
