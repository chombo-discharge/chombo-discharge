/*!
  @file   advection_diffusion_tagger.cpp
  @brief  Implementation of advection_diffusion_tagger.H
  @author Robert Marskar
  @date   Nov. 2017
  @todo   Rename to advection_diffusion_tagger. 
*/

#include "advection_diffusion_tagger.H"
#include "data_ops.H"

#include <ParmParse.H> 

using namespace physics::advection_diffusion;

advection_diffusion_tagger::advection_diffusion_tagger(RefCountedPtr<cdr_solver>& a_solver,
						       RefCountedPtr<amr_mesh>&   a_amr){
  m_solver = a_solver;
  m_amr    = a_amr;
  m_name   = "advection_diffusion";
}

advection_diffusion_tagger::~advection_diffusion_tagger(){

}

void advection_diffusion_tagger::regrid(){
  pout() << "regridding cell tagger" << endl;
}

void advection_diffusion_tagger::parse_options(){
  ParmParse pp(m_name.c_str());
  pp.get("refine_curv", m_refi_curv);
  pp.get("refine_magn", m_refi_magn);
  
  parse_buffer(); // Derived from cell
}

bool advection_diffusion_tagger::tag_cells(EBAMRTags& a_tags){
  EBAMRCellData sca;
  EBAMRCellData vec;

  m_amr->allocate(sca, phase::gas, 1);
  m_amr->allocate(vec, phase::gas, SpaceDim);

  const EBAMRCellData& state = m_solver->get_state();

  // Compute the gradient, vec = grad(phi)
  m_amr->compute_gradient(vec, state, phase::gas); // vec = grad(phi)
  data_ops::vector_length(sca, vec);               // sca = |grad(phi)|
  data_ops::set_covered_value(sca, 0, 0.0);        // covered cell values are set to 0.0

  bool found_tags = false;

  // Never tag on max_amr_depth
  const int finest_level     = m_amr->get_finest_level();
  const int max_depth        = m_amr->get_max_amr_depth();
  const int finest_tag_level = (finest_level == max_depth) ? max_depth - 1 : finest_level; // Never tag on max_amr_depth

  for (int lvl = 0; lvl <= finest_tag_level; lvl++){
    data_ops::scale(*sca[lvl], m_amr->get_dx()[lvl]); // sca = |grad(phi)|*dx

    const Real SAFETY = 1.E-6;
    
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBCellFAB& c   = (*sca[lvl])[dit()];
      const EBCellFAB& phi = (*state[lvl])[dit()];

      const BaseFab<Real>& cReg   = c.getSingleValuedFAB();
      const BaseFab<Real>& phiReg = phi.getSingleValuedFAB();

      const Box box = dbl.get(dit());


      // These are the tags
      DenseIntVectSet& tags = (*a_tags[lvl])[dit()];

      // Do regular cells
      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv = bit();

	const Real crit = Abs(cReg(iv,0))/(SAFETY + Abs(phiReg(iv,0)));
	if(crit > m_refi_curv && Abs(phiReg(iv,0)) > m_refi_magn){
	  tags |= iv;
	  found_tags = true;
	}
      }

      tags &= box;

      // Do irregular cells
      const EBISBox& ebisbox = c.getEBISBox();
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      for (VoFIterator vofit(ebisbox.getIrregIVS(box), ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
      }
    }
  }

  // Some ranks may have gotten new tags while others have not. This little code snippet
  // sets got_new_tags = true for all ranks if any rank originally had got_new_tags = true
#ifdef CH_MPI
  int glo = 1;
  int loc = found_tags ? 1 : 0;

  const int result = MPI_Allreduce(&loc, &glo, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);

  found_tags = (glo == 1) ? true : false;
#endif
    
  return true;
}
