/*!
  @file   ito_plasma_tagger.cpp
  @brief  Implementation of ito_plasma_tagger.H
  @author Robert marskar
  @date   May. 2018
*/

#include "ito_plasma_tagger.H"
#include "ito_iterator.H"
#include "rte_iterator.H"
#include "data_ops.H"

#include <EBArith.H>
#include <ParmParse.H>

#include "CD_NamespaceHeader.H"
using namespace physics::ito_plasma;

ito_plasma_tagger::ito_plasma_tagger(){
  CH_TIME("ito_plasma_tagger::ito_plasma_tagger");
  m_verbosity = 10;
  if(m_verbosity > 5){
    pout() << "ito_plasma_tagger::ito_plasma_tagger" << endl;
  }

  m_name  = "ito_plasma_tagger";
  m_phase = phase::gas;
  m_Realm = Realm::Primal;
}

ito_plasma_tagger::ito_plasma_tagger(const RefCountedPtr<ito_plasma_physics>&     a_physics,
				     const RefCountedPtr<ito_plasma_stepper>&     a_timeStepper,
				     const RefCountedPtr<AmrMesh>&               a_amr,
				     const RefCountedPtr<computational_geometry>& a_computationalGeometry) : ito_plasma_tagger() {
  this->define(a_physics, a_timeStepper, a_amr, a_computationalGeometry);
}

ito_plasma_tagger::~ito_plasma_tagger(){

}

void ito_plasma_tagger::define(const RefCountedPtr<ito_plasma_physics>&     a_physics,
			       const RefCountedPtr<ito_plasma_stepper>&     a_timeStepper,
			       const RefCountedPtr<AmrMesh>&               a_amr,
			       const RefCountedPtr<computational_geometry>& a_computationalGeometry){
  CH_TIME("ito_plasma_tagger::define");
  if(m_verbosity > 5){
    pout() << m_name + "::define" << endl;
  }

  m_physics     = a_physics;
  m_timeStepper = a_timeStepper;
  m_amr         = a_amr;
  m_computationalGeometry    = a_computationalGeometry;
}

void ito_plasma_tagger::regrid(){
  CH_TIME("ito_plasma_tagger::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }
  
  if(m_num_tracers > 0){
    m_tracer.resize(m_num_tracers);
    m_grad_tracer.resize(m_num_tracers);
    for (int i = 0; i < m_num_tracers; i++){
      m_amr->allocate(m_tracer[i],      m_Realm, m_phase, 1);
      m_amr->allocate(m_grad_tracer[i], m_Realm, m_phase, SpaceDim);
    }
  }
}

void ito_plasma_tagger::set_phase(const phase::which_phase a_phase){
  CH_TIME("ito_plasma_tagger::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::set_phase" << endl;
  }
  
  m_phase = a_phase;
}

int ito_plasma_tagger::getNumberOfPlotVariables(){
  CH_TIME("ito_plasma_tagger::getNumberOfPlotVariables_cells");
  if(m_verbosity > 5){
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }
  
  return m_num_tracers;
}

Vector<EBAMRCellData>& ito_plasma_tagger::get_tracer_fields() {
  return m_tracer;
}

void ito_plasma_tagger::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) {
  CH_TIME("ito_plasma_tagger::writePlotData");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotData" << endl;
  }

  this->compute_tracers();
  for (int i = 0; i < m_num_tracers; i++){
    std::string one = "Tracer field-";
    long int j = i;
    char s[2]; 
    sprintf(s,"%ld", j);
    std::string two(s);
      
    const EBAMRCellData& tracer = m_tracer[i];
    
    const Interval src_interv(0, 0);
    const Interval dst_interv(a_icomp, a_icomp);
    
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      tracer[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
      data_ops::set_covered_value(*a_output[lvl], a_icomp, 0.0);
    }

    // Add component and name
    a_plotVariableNames.push_back(one+two);
    a_icomp++;
  }
}

bool ito_plasma_tagger::tagCells(EBAMRTags& a_tags){
  CH_TIME("ito_plasma_tagger::tagCells");
  if(m_verbosity > 5){
    pout() << m_name + "::tagCells" << endl;
  }

  bool got_new_tags = false;

  const RealVect origin      = m_amr->getProbLo();
  const Real time            = m_timeStepper->get_time();
  const int finest_level     = m_amr->getFinestLevel();
  const int max_depth        = m_amr->getMaxAmrDepth();
  const int finest_tag_level = (finest_level == max_depth) ? max_depth - 1 : finest_level; // Never tag on max_amr_depth

  if(m_num_tracers > 0){

    // Compute tracer fields. This is an overriden pure function
    compute_tracers();
    
    for (int lvl = 0; lvl <= finest_tag_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_Realm)[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_Realm, m_phase)[lvl];
      const Real dx                = m_amr->getDx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const Box box          = dbl.get(dit());
	const EBISBox& ebisbox = ebisl[dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();

	const IntVectSet irreg_ivs = ebisbox.getIrregIVS(box);
	const IntVectSet prev_tags = IntVectSet((*a_tags[lvl])[dit()]);

	DenseIntVectSet coarsen_tags(box, false); // Cells that will be coarsened
	DenseIntVectSet refine_tags(box, false);  // Cells that will be refined

	Vector<EBCellFAB*> tracers;
	Vector<EBCellFAB*> gtracers;

	for (int i = 0; i < m_num_tracers; i++){
	  tracers.push_back(&((*m_tracer[i][lvl])[dit()]));
	  gtracers.push_back(&((*m_grad_tracer[i][lvl])[dit()]));
	}

	DenseIntVectSet& tags = (*a_tags[lvl])[dit()];
	
	// Refinement and coarsening
	refine_cells_box(refine_tags, tracers, gtracers, lvl, box, ebisbox, time, dx, origin);
	coarsen_cells_box(coarsen_tags, tracers, gtracers, lvl, box, ebisbox, time, dx, origin);

	// Check if we got any new tags, or we are just recycling old tags.
	// Basically we will check if (current_tags + refined_tags - coarsen_tags) == current_tags
	DenseIntVectSet cpy1 = tags;
	tags -= coarsen_tags;
	tags |= refine_tags;
	DenseIntVectSet cpy2 = tags;

	cpy2 -= cpy1; // = new tags minus old tags. If nonzero, we got some new tags. 
	cpy1 -= tags; // = old_tags minus new tags. If nonzero, we got some new tags
	if(cpy1.numPts() != 0 || cpy2.numPts() != 0){
	  got_new_tags = true;
	}

	tags &= box;
      }
    }
  }

  // Some ranks may have gotten new tags while others have not. This little code snippet
  // sets got_new_tags = true for all ranks if any rank originally had got_new_tags = true
#ifdef CH_MPI
  int glo = 1;
  int loc = got_new_tags ? 1 : 0;

  const int result = MPI_Allreduce(&loc, &glo, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);

  got_new_tags = (glo == 1) ? true : false;
#endif

  return got_new_tags;
}

void ito_plasma_tagger::refine_cells_box(DenseIntVectSet&          a_refined_tags,
					 const Vector<EBCellFAB*>& a_tracers,
					 const Vector<EBCellFAB*>& a_grad_tracers,
					 const int                 a_lvl,
					 const Box                 a_box,
					 const EBISBox&            a_ebisbox,
					 const Real                a_time,
					 const Real                a_dx,
					 const RealVect            a_origin){


  
  Vector<BaseFab<Real>* > reg_tracers;
  Vector<BaseFab<Real>* > reg_gtracer;

  for (int i = 0; i < m_num_tracers; i++){
    reg_tracers.push_back(&(a_tracers[i]->getSingleValuedFAB()));
    reg_gtracer.push_back(&(a_grad_tracers[i]->getSingleValuedFAB()));
  }

  // Regular box loop
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv   = bit();
    const RealVect pos = a_origin + a_dx*RealVect(iv) + 0.5*a_dx*RealVect::Unit;
    
    // If position is inside any of the tagging boxes, we can refine
    if(insideTagBox(pos) && a_ebisbox.isRegular(iv)){
      
      // Build point-wise tracer fields
      Vector<Real>     tr(m_num_tracers); 
      Vector<RealVect> gt(m_num_tracers);
      for (int i = 0; i < m_num_tracers; i++){
	tr[i] = (*reg_tracers[i])(iv, 0);
	gt[i] = RealVect(D_DECL((*reg_gtracer[i])(iv, 0),(*reg_gtracer[i])(iv, 1),(*reg_gtracer[i])(iv, 2)));
      }

      // Check if this cell should be refined
      const bool refine = refine_cell(pos, a_time, a_dx, a_lvl, tr, gt);

      // If we refine, grow with buffer and increment to a_refined_tags
      if(refine){
	a_refined_tags |= iv;
      }
    }
  }

  // Irregular box loop
  const IntVectSet& irreg = a_ebisbox.getIrregIVS(a_box);
  const EBGraph& ebgraph  = a_ebisbox.getEBGraph();
  for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect pos  = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, a_origin);

    // If position is inside any of the tagging boxes, we can refine
    if(insideTagBox(pos)){
      
      // Build point-wise tracer fields
      Vector<Real>     tr(m_num_tracers); 
      Vector<RealVect> gt(m_num_tracers);
      for (int i = 0; i < m_num_tracers; i++){
	tr[i] = (*a_tracers[i])(vof, 0);
	gt[i] = RealVect(D_DECL((*a_grad_tracers[i])(vof, 0),(*a_grad_tracers[i])(vof, 1),(*a_grad_tracers[i])(vof, 2)));
      }

      // Check if this cell should be refined
      const bool refine = refine_cell(pos, a_time, a_dx, a_lvl, tr, gt);

      if(refine){
	a_refined_tags |= vof.gridIndex();
      }
    }
  }
}

void ito_plasma_tagger::coarsen_cells_box(DenseIntVectSet&         a_coarsened_tags,
					  const Vector<EBCellFAB*>& a_tracers,
					  const Vector<EBCellFAB*>& a_grad_tracers,
					  const int                 a_lvl,
					  const Box                 a_box,
					  const EBISBox&            a_ebisbox,
					  const Real                a_time,
					  const Real                a_dx,
					  const RealVect            a_origin){


  
  Vector<BaseFab<Real>* > reg_tracers;
  Vector<BaseFab<Real>* > reg_gtracer;

  for (int i = 0; i < m_num_tracers; i++){
    reg_tracers.push_back(&(a_tracers[i]->getSingleValuedFAB()));
    reg_gtracer.push_back(&(a_grad_tracers[i]->getSingleValuedFAB()));
  }

  // Regular box loop
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv   = bit();
    const RealVect pos = a_origin + a_dx*RealVect(iv);
    
    // If position is inside any of the tagging boxes, we can refine
    const bool inside = insideTagBox(pos);
    if(!inside){
      a_coarsened_tags |= iv;
    }
    else if(inside && a_ebisbox.isRegular(iv)){
      
      // Build point-wise tracer fields
      Vector<Real>     tr(m_num_tracers); 
      Vector<RealVect> gt(m_num_tracers);
      for (int i = 0; i < m_num_tracers; i++){
	tr[i] = (*reg_tracers[i])(iv, 0);
	gt[i] = RealVect(D_DECL((*reg_gtracer[i])(iv, 0),(*reg_gtracer[i])(iv, 1),(*reg_gtracer[i])(iv, 2)));
      }

      // Check if this cell should be refined
      const bool coarsen = coarsen_cell(pos, a_time, a_dx, a_lvl, tr, gt);

      // If we refine, grow with buffer and increment to a_refined_tags
      if(coarsen){
	a_coarsened_tags |= iv;
      }
    }
  }

  // Irregular box loop
  const IntVectSet& irreg = a_ebisbox.getIrregIVS(a_box);
  const EBGraph& ebgraph  = a_ebisbox.getEBGraph();
  for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect pos  = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, a_origin);

    // If position is inside any of the tagging boxes, we can refine
    const bool inside = insideTagBox(pos);
    if(!inside){
      a_coarsened_tags |= vof.gridIndex();
    }
    else if(inside){
      
      // Build point-wise tracer fields
      Vector<Real>     tr(m_num_tracers); 
      Vector<RealVect> gt(m_num_tracers);
      for (int i = 0; i < m_num_tracers; i++){
	tr[i] = (*a_tracers[i])(vof, 0);
	gt[i] = RealVect(D_DECL((*a_grad_tracers[i])(vof, 0),(*a_grad_tracers[i])(vof, 1),(*a_grad_tracers[i])(vof, 2)));
      }

      // Check if this cell should be refined
      const bool coarsen = coarsen_cell(pos, a_time, a_dx, a_lvl, tr, gt);

      if(coarsen){
	a_coarsened_tags |= vof.gridIndex();
      }
    }
  }
}
#include "CD_NamespaceFooter.H"
