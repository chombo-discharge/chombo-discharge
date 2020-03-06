/*!
  @file   cell_tagger.cpp
  @brief  Implementation of cell_tagger.H
  @author Robert marskar
  @date   May. 2018
*/

#include "cell_tagger.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"
#include "data_ops.H"

#include <EBArith.H>
#include <ParmParse.H>

cell_tagger::cell_tagger(){
  CH_TIME("cell_tagger::cell_tagger");
  m_verbosity = 10;
  if(m_verbosity > 5){
    pout() << "cell_tagger::cell_tagger" << endl;
  }

  m_name  = "cell_tagger";
  m_phase = phase::gas;

  // Parse class options directly
#if 0 // Should be moved
  parse_verbosity();
  parse_buffer();
  parse_boxes();
#endif
}

cell_tagger::~cell_tagger(){

}

int cell_tagger::get_num_tracers(){
  return m_num_tracers;
}

int cell_tagger::get_buffer(){
  return m_buffer;
}

void cell_tagger::parse_boxes(){
  
  ParmParse pp(m_name.c_str());

  int num_boxes = 0;
  pp.get("num_boxes", num_boxes);

  m_tagboxes.resize(0);
  if(num_boxes > 0){
    m_tagboxes.resize(num_boxes);

    const int ndigits = (int) log10((double) num_boxes) + 1;
      
    for (int ibox = 0; ibox < num_boxes; ibox++){
      char* cstr = new char[ndigits];
      sprintf(cstr, "%d", 1+ibox);

      std::string str1 = "box" + std::string(cstr) + "_lo";
      std::string str2 = "box" + std::string(cstr) + "_hi";

      Vector<Real> corner_lo(SpaceDim);
      Vector<Real> corner_hi(SpaceDim);

      pp.getarr(str1.c_str(), corner_lo, 0, SpaceDim);
      pp.getarr(str2.c_str(), corner_hi, 0, SpaceDim);

      const RealVect c1 = RealVect(D_DECL(corner_lo[0], corner_lo[1], corner_lo[2]));
      const RealVect c2 = RealVect(D_DECL(corner_hi[0], corner_hi[1], corner_hi[2]));

      m_tagboxes[ibox] = real_box(c1,c2);
	
      delete cstr;
    }
  }
}

void cell_tagger::parse_buffer(){

  ParmParse pp(m_name.c_str());
  pp.get("buffer", m_buffer);
  m_buffer = Max(0, m_buffer);
}

void cell_tagger::define(const RefCountedPtr<plasma_kinetics>&        a_plaskin,
			 const RefCountedPtr<time_stepper>&           a_timestepper,
			 const RefCountedPtr<amr_mesh>&               a_amr,
			 const RefCountedPtr<computational_geometry>& a_compgeom,
			 const RefCountedPtr<physical_domain>&        a_physdom){
  CH_TIME("cell_tagger::define");
  if(m_verbosity > 5){
    pout() << m_name + "::define" << endl;
  }

  m_plaskin     = a_plaskin;
  m_timestepper = a_timestepper;
  m_amr         = a_amr;
  m_compgeom    = a_compgeom;
  m_physdom     = a_physdom;
}

void cell_tagger::regrid(){
  CH_TIME("cell_tagger::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }
  
  if(m_num_tracers > 0){
    m_tracer.resize(m_num_tracers);
    m_grad_tracer.resize(m_num_tracers);
    for (int i = 0; i < m_num_tracers; i++){
      m_amr->allocate(m_tracer[i],      m_phase, 1);
      m_amr->allocate(m_grad_tracer[i], m_phase, SpaceDim);
    }
  }
}

void cell_tagger::set_phase(const phase::which_phase a_phase){
  CH_TIME("cell_tagger::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::set_phase" << endl;
  }
  
  m_phase = a_phase;
}

void cell_tagger::parse_verbosity(){
  CH_TIME("cell_tagger::parse_verbosity");

  ParmParse pp(m_name.c_str());
  pp.get("verbosity", m_verbosity);
}

Vector<EBAMRCellData>& cell_tagger::get_tracer_fields() {
  return m_tracer;
}

bool cell_tagger::tag_cells(EBAMRTags& a_tags){
  CH_TIME("cell_tagger::tag_cells");
  if(m_verbosity > 5){
    pout() << m_name + "::tag_cells" << endl;
  }

  bool got_new_tags = false;

  const RealVect origin      = m_physdom->get_prob_lo();
  const Real time            = m_timestepper->get_time();
  const int finest_level     = m_amr->get_finest_level();
  const int max_depth        = m_amr->get_max_amr_depth();
  const int finest_tag_level = (finest_level == max_depth) ? max_depth - 1 : finest_level; // Never tag on max_amr_depth

  if(m_num_tracers > 0){

    // Compute tracer fields. This is a virtual function. 
    compute_tracers();
    
    for (int lvl = 0; lvl <= finest_tag_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
      const Real dx                = m_amr->get_dx()[lvl];

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


void cell_tagger::refine_cells_box(DenseIntVectSet&          a_refined_tags,
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
    if(inside_tag_box(pos) && a_ebisbox.isRegular(iv)){
      
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
    if(inside_tag_box(pos)){
      
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

void cell_tagger::coarsen_cells_box(DenseIntVectSet&         a_coarsened_tags,
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
    const bool inside = inside_tag_box(pos);
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
    const bool inside = inside_tag_box(pos);
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


bool cell_tagger::inside_tag_box(const RealVect a_pos){
  bool do_this_refine = (m_tagboxes.size() > 0) ? false : true;
  for (int ibox = 0; ibox < m_tagboxes.size(); ibox++){
    const RealVect lo = m_tagboxes[ibox].get_lo();
    const RealVect hi = m_tagboxes[ibox].get_hi();

    if(a_pos >= lo && a_pos <= hi){
      do_this_refine = true;
    }
  }

  return do_this_refine;
}
