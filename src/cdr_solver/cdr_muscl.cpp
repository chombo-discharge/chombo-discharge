/*!
  @file   cdr_muscl.cpp
  @brief  Implementation of cdr_muscl.H
  @author Robert Marskar
  @date   Feb. 2018
  @todo   The slope computation and the upwinding should definitely be done in Fortran....
*/

#include "cdr_muscl.H" 
#include "slope_limiters.H"
#include "data_ops.H"

#include <ParmParse.H>

cdr_muscl::cdr_muscl(){

  m_slope_func = &(slope_limiters::superbee);

  { // Select slope from input
    ParmParse pp("cdr_muscl");
    if(pp.contains("slope_limiter")){
      std::string str;
      pp.get("slope_limiter", str);
      if(str == "koren"){
	m_slope_func = &(slope_limiters::koren);
      }
      else if(str == "minmod"){
	m_slope_func = &(slope_limiters::minmod);
      }
      else if(str == "superbee"){
	m_slope_func = &(slope_limiters::superbee);
      }
      else if(str == "van_leer"){
	m_slope_func = &(slope_limiters::van_leer);
      }
      else {
	MayDay::Abort("cdr_muscl::cdr_muscl - unsupported slope limiter requested"); 
      }
    }
  }
}

cdr_muscl::~cdr_muscl(){

}

int cdr_muscl::query_ghost() const {
  return 3;
}

void cdr_muscl::advect_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state, const Real a_extrap_dt){
  CH_TIME("cdr_muscl::advect_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::advect_to_faces" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  
  EBAMRCellData copy_state; // Need to overwrite ghost cells, so I do that here.
  EBAMRCellData slopes;     // Cell-centered slopes

  m_amr->allocate(copy_state, m_phase, ncomp);
  m_amr->allocate(slopes,     m_phase, SpaceDim);

  data_ops::set_value(copy_state, 0.0);
  data_ops::set_value(slopes,     0.0);
  data_ops::incr(copy_state, a_state, 1.0);

  m_amr->interp_ghost_pwl(copy_state, m_phase);


  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      EBFluxFAB& face_state  = (*a_face_state[lvl])[dit()];
      EBCellFAB& slopes_fab  = (*slopes[lvl])[dit()];
      const Box box          = dbl.get(dit());
      const EBCellFAB& state = (*copy_state[lvl])[dit()];

      // Here's the strategy:
      // 2. Compute cell-centered slopes. For boundary faces, set covered state to zero when computing the slope. 
      // 3. Iterate over faces, compute left and right eigenstates. Select upwind state.

      this->compute_slopes(slopes_fab, state, box);
      this->upwind(face_state, slopes_fab, state, box);
    }
  }
}

void cdr_muscl::compute_slopes(EBCellFAB& a_slopes, const EBCellFAB& a_state, const Box& a_box){
  CH_TIME("cdr_muscl::compute_slopes");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_slopes" << endl;
  }

  const EBISBox& ebisbox = a_state.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const IntVectSet ivs(a_box);

  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();

    for (int dir = 0; dir < SpaceDim; dir++){

      // Strategy: Get left and right states, compute left and right slopes.
      // PS: I really need to get this right wrt. irregular cells (and also multicells)
    }
  }

  MayDay::Abort("cdr_muscl::compute_slopes - error");
  
}

void cdr_muscl::upwind(EBFluxFAB& a_face_states, const EBCellFAB& a_slopes, const EBCellFAB& a_state, const Box& a_box){
  CH_TIME("cdr_muscl::upwind");
  if(m_verbosity > 5){
    pout() << m_name + "::upwind" << endl;
  }

  MayDay::Abort("cdr_muscl::upwind - error");
}
