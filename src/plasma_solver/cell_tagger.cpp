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

cell_tagger::cell_tagger(const int a_num_tracers){
  CH_TIME("cell_tagger::cell_tagger");
  this->set_verbosity(-1);
  if(m_verbosity > 5){
    pout() << "cell_tagger::cell_tagger" << endl;
  }
  
  m_num_tracers = a_num_tracers;
  m_name        = "cell_tagger";


  this->set_phase(phase::gas);
}

cell_tagger::~cell_tagger(){

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

void cell_tagger::tag_cells(EBAMRTags& a_tags){
  CH_TIME("cell_tagger::tag_cells");
  if(m_verbosity > 5){
    pout() << m_name + "::tag_cells" << endl;
  }

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

	const IntVectSet& irreg_ivs = ebisbox.getIrregIVS(box);
	const IntVectSet prev_tags  = IntVectSet((*a_tags[lvl])[dit()].get_ivs());

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

	  if(refine){
	    refine_tags |= vof.gridIndex();
	  }
	}

	DenseIntVectSet& tags = (*a_tags[lvl])[dit()].get_ivs();
	tags -= coarsen_tags;
	tags |= refine_tags;
	tags &= box;
      }
    }
  }
}

Vector<EBAMRCellData>& cell_tagger::get_tracer_fields() {
  return m_tracer;
}
