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

cell_tagger::cell_tagger(const int a_num_tracers){
  CH_TIME("cell_tagger::cell_tagger");
  this->set_verbosity(-1);
  if(m_verbosity > 5){
    pout() << "cell_tagger::cell_tagger" << endl;
  }

  m_tagboxes.resize(0);
  m_num_tracers = a_num_tracers;
  m_name        = "cell_tagger";
  m_buffer      = 0;


  this->set_phase(phase::gas);

  { // Get options from input script
    ParmParse pp("cell_tagger");

    int num_boxes = 0;
    pp.query("num_boxes", num_boxes);
    pp.query("buffer",    m_buffer);

    m_buffer = Max(0, m_buffer);

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
  

}

cell_tagger::~cell_tagger(){

}

int cell_tagger::get_num_tracers(){
  return m_num_tracers;
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

void cell_tagger::set_verbosity(const int a_verbosity){
  CH_TIME("cell_tagger::set_verbosity");
  m_verbosity = a_verbosity;
  
  if(m_verbosity > 5){
    pout() << m_name + "::set_verbosity" << endl;
  }
}

bool cell_tagger::tag_cells(EBAMRTags& a_tags){
  CH_TIME("cell_tagger::tag_cells");
  if(m_verbosity > 5){
    pout() << m_name + "::tag_cells" << endl;
  }

  bool got_new_tags = false;

  if(m_num_tracers > 0){
    
    const RealVect origin      = m_physdom->get_prob_lo();
    const Real time            = m_timestepper->get_time();
    const int finest_level     = m_amr->get_finest_level();
    const int max_depth        = m_amr->get_max_amr_depth();
    const int finest_tag_level = (finest_level == max_depth) ? max_depth - 1 : finest_level; // Never tag on max_amr_depth

    this->compute_tracers();

    for (int lvl = 0; lvl <= finest_tag_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
      const Real dx                = m_amr->get_dx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const Box box          = dbl.get(dit());
	const EBISBox& ebisbox = ebisl[dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();

	const IntVectSet irreg_ivs = ebisbox.getIrregIVS(box);
	const IntVectSet prev_tags = IntVectSet((*a_tags[lvl])[dit()].get_ivs());

	DenseIntVectSet coarsen_tags(box, false);
	DenseIntVectSet refine_tags(box, false);

	// Coarsening loop - do not coarsen irregular cells that have been tagged previously (we consider them to be too important)
	const IntVectSet coarsen_ivs = (prev_tags - irreg_ivs); 
	for (VoFIterator vofit(coarsen_ivs, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	  const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin);

	  Vector<Real> tracers(m_num_tracers);
	  Vector<RealVect> grad_tracers(m_num_tracers);
	  
	  for (int itracer = 0; itracer < m_num_tracers; itracer++){
	    tracers[itracer]     = (*m_tracer[itracer][lvl])[dit()](vof, 0);
	    grad_tracers[itracer] = RealVect(D_DECL((*m_grad_tracer[itracer][lvl])[dit()](vof, 0),
						    (*m_grad_tracer[itracer][lvl])[dit()](vof, 1),
						    (*m_grad_tracer[itracer][lvl])[dit()](vof, 2)));
	  }
	  const bool coarsen = this->coarsen_cell(pos,
						  time,
						  dx,
						  lvl,
						  tracers,
						  grad_tracers);

	  if(coarsen){
	    coarsen_tags |= vof.gridIndex();
	  }
	}

	// Refinement loop
	const IntVectSet refine_ivs(box);
	for (VoFIterator vofit(refine_ivs, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	  const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, origin);

	  Vector<Real> tracers(m_num_tracers);
	  Vector<RealVect> grad_tracers(m_num_tracers);
	  
	  for (int itracer = 0; itracer < m_num_tracers; itracer++){
	    tracers[itracer]     = (*m_tracer[itracer][lvl])[dit()](vof, 0);
	    grad_tracers[itracer] = RealVect(D_DECL((*m_grad_tracer[itracer][lvl])[dit()](vof, 0),
						    (*m_grad_tracer[itracer][lvl])[dit()](vof, 1),
						    (*m_grad_tracer[itracer][lvl])[dit()](vof, 2)));
	  }
	  const bool refine = this->refine_cell(pos,
						time,
						dx,
						lvl,
						tracers,
						grad_tracers);

	  bool do_this_refine = (m_tagboxes.size() > 0) ? false : true;
	  for (int ibox = 0; ibox < m_tagboxes.size(); ibox++){
	    const RealVect lo = m_tagboxes[ibox].get_lo();
	    const RealVect hi = m_tagboxes[ibox].get_hi();

	    if(pos >= lo && pos <= hi){
	      do_this_refine = true;
	    }
	  }

	  if(refine && do_this_refine){
	    IntVectSet buf(vof.gridIndex());
	      buf.grow(m_buffer);
	      buf &= box;
	    for (IVSIterator ivs_it(buf); ivs_it.ok(); ++ivs_it){
	      refine_tags |= ivs_it();
	    }
	  }
	}


	DenseIntVectSet& tags    = (*a_tags[lvl])[dit()].get_ivs();
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

#ifdef CH_MPI
  int glo = 1;
  int loc = got_new_tags ? 1 : 0;

  const int result = MPI_Allreduce(&loc, &glo, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);

  got_new_tags = (glo == 1) ? true : false;
#endif

  return got_new_tags;
}

Vector<EBAMRCellData>& cell_tagger::get_tracer_fields() {
  return m_tracer;
}
