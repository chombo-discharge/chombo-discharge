/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_AdvectionDiffusionTagger.cpp
  @brief  Implementation of CD_AdvectionDiffusionTagger.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_AdvectionDiffusionTagger.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::AdvectionDiffusion;

AdvectionDiffusionTagger::AdvectionDiffusionTagger(RefCountedPtr<CdrSolver>& a_solver,
						   RefCountedPtr<AmrMesh>&   a_amr){
  m_solver    = a_solver;
  m_amr       = a_amr;
  m_name      = "AdvectionDiffusion";
  m_verbosity = -1;
  m_realm     = a_solver->getRealm();
}


AdvectionDiffusionTagger::~AdvectionDiffusionTagger(){

}

void AdvectionDiffusionTagger::regrid(){

}

void AdvectionDiffusionTagger::parseOptions(){
  ParmParse pp(m_name.c_str());
  pp.get("refine_curv", m_refi_curv);
  pp.get("refine_magn", m_refi_magn);
  
  parseBuffer(); // Derived from cell
}

bool AdvectionDiffusionTagger::tagCells(EBAMRTags& a_tags){
  EBAMRCellData sca;
  EBAMRCellData vec;

  m_amr->allocate(sca, m_realm, phase::gas, 1);
  m_amr->allocate(vec, m_realm, phase::gas, SpaceDim);

  const EBAMRCellData& state = m_solver->getPhi();

  // Compute the gradient, vec = grad(phi)
  m_amr->computeGradient(vec, state, m_realm, phase::gas); // vec = grad(phi)
  DataOps::vectorLength(sca, vec);                        // sca = |grad(phi)|
  DataOps::set_covered_value(sca, 0, 0.0);                 // covered cell values are set to 0.0

  bool found_tags = false;

  // Never tag on max_amr_depth
  const int finest_level     = m_amr->getFinestLevel();
  const int max_depth        = m_amr->getMaxAmrDepth();
  const int finest_tag_level = (finest_level == max_depth) ? max_depth - 1 : finest_level; // Never tag on max_amr_depth

  for (int lvl = 0; lvl <= finest_tag_level; lvl++){
    DataOps::scale(*sca[lvl], m_amr->getDx()[lvl]); // sca = |grad(phi)|*dx

    const Real SAFETY = 1.E-6;
    
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBCellFAB& c   = (*sca[lvl])[dit()];
      const EBCellFAB& phi = (*state[lvl])[dit()];

      const BaseFab<Real>& cReg   = c.getSingleValuedFAB();
      const BaseFab<Real>& phiReg = phi.getSingleValuedFAB();

      const Box box = dbl.get(dit());

      // These are the tags
      DenseIntVectSet& tags = (*a_tags[lvl])[dit()];
      tags.makeEmptyBits(); // Clear all previous tags

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
  return found_tags;
}

#include <CD_NamespaceFooter.H>
