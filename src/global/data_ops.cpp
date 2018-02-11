/*!
  @file data_ops.cpp
  @brief Implementation of data_ops.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "data_ops.H"
#include "EBLevelDataOps.H"
#include "MFLevelDataOps.H"

void data_ops::average_cell_to_face(EBAMRFluxData&               a_facedata,
				    const EBAMRCellData&         a_celldata,
				    const Vector<ProblemDomain>& a_domains){
  for (int lvl = 0; lvl < a_facedata.size(); lvl++){
    data_ops::average_cell_to_face(*a_facedata[lvl], *a_celldata[lvl], a_domains[lvl]);
  }
}

void data_ops::average_cell_to_face(LevelData<EBFluxFAB>&       a_facedata,
				    const LevelData<EBCellFAB>& a_celldata,
				    const ProblemDomain&        a_domain){

  CH_assert(a_facedata.nComp() == 1);
  CH_assert(a_celldata.nComp() == SpaceDim);
  
  for (DataIterator dit = a_facedata.dataIterator(); dit.ok(); ++dit){
    EBFluxFAB& flux_vel       = a_facedata[dit()];
    const EBCellFAB& cell_vel = a_celldata[dit()];
    const EBISBox& ebisbox    = cell_vel.getEBISBox();
    const EBGraph& ebgraph    = ebisbox.getEBGraph();
    const Box& box            = a_celldata.disjointBoxLayout().get(dit());
    
    for (int dir = 0; dir < SpaceDim; dir++){
      EBLevelDataOps::averageCellToFace(flux_vel[dir], cell_vel, ebgraph, box, dir, dir, a_domain, dir, 0);
    }
  }
}

void data_ops::average_cell_to_face_allcomps(LevelData<EBFluxFAB>&       a_facedata,
					     const LevelData<EBCellFAB>& a_celldata,
					     const ProblemDomain&        a_domain){

  CH_assert(a_facedata.nComp() == a_celldata.nComp());
  
  const int ncomp = a_facedata.nComp();

  for (DataIterator dit = a_facedata.dataIterator(); dit.ok(); ++dit){
    EBFluxFAB& facedata       = a_facedata[dit()];
    const EBCellFAB& celldata = a_celldata[dit()];
    const EBISBox& ebisbox    = celldata.getEBISBox();
    const EBGraph& ebgraph    = ebisbox.getEBGraph();
    const Box& box            = a_celldata.disjointBoxLayout().get(dit());
    const IntVectSet ivs(box);

    // Interior faces
    FaceStop::WhichFaces stop_crit;
    for (int dir = 0; dir < SpaceDim; dir++){

      stop_crit = FaceStop::SurroundingWithBoundary;
      for (FaceIterator faceit(ivs, ebgraph, dir, stop_crit); faceit.ok(); ++faceit){
	const FaceIndex& face  = faceit();
	const VolIndex& vof_lo = face.getVoF(Side::Lo);
	const VolIndex& vof_hi = face.getVoF(Side::Hi);

	if(!face.isBoundary()){
	  for (int comp = 0; comp < ncomp; comp++){
	    facedata[dir](face, comp) = 0.5*(celldata(vof_lo, comp) + celldata(vof_hi, comp));
	  }
	}
	else{
	  VolIndex vof;
	  if(a_domain.contains(face.gridIndex(Side::Lo))){
	    vof = face.getVoF(Side::Lo);
	  }
	  else if(a_domain.contains(face.gridIndex(Side::Hi))){
	    vof = face.getVoF(Side::Hi);
	  }
	  else{
	    MayDay::Error("data_ops::average_cell_to_face - logic bust in average cells to faces");
	  }
	  for (int comp = 0; comp < ncomp; comp++){
	    facedata[dir](face, comp) = celldata(vof, comp);
	  }
	}
      }
    }
  }
}

void data_ops::incr(MFAMRCellData& a_lhs, const MFAMRCellData& a_rhs, const Real a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}
  
void data_ops::incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const Real a_scale){
  MFLevelDataOps::incr(a_lhs, a_rhs, a_scale);
}

void data_ops::incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void data_ops::incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real& a_scale){
  EBLevelDataOps::incr(a_lhs, a_rhs, a_scale);
}

void data_ops::incr(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void data_ops::incr(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs, const Real& a_scale){
  CH_assert(a_lhs.nComp() == a_rhs.nComp());

  const int ncomp = a_lhs.nComp();
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    BaseIVFAB<Real>& lhs       = a_lhs[dit()];
    const BaseIVFAB<Real>& rhs = a_rhs[dit()];

    for (VoFIterator vofit(lhs.getIVS(), lhs.getEBGraph()); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      for (int comp = 0; comp < ncomp; comp++){
	lhs(vof, comp) += rhs(vof, comp)*a_scale;
      }
    }
  }
}

void data_ops::incr(EBAMRCellData& a_lhs, const EBAMRIVData& a_rhs, const Real a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void data_ops::incr(LevelData<EBCellFAB>& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs, const Real a_scale){
  CH_assert(a_lhs.nComp() == a_rhs.nComp());

  const int ncomp = a_lhs.nComp();

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs = a_lhs[dit()];
    const BaseIVFAB<Real>& rhs = a_rhs[dit()];
    const EBGraph& ebgraph     = rhs.getEBGraph();
    const IntVectSet& ivs      = rhs.getIVS();

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      for (int comp = 0; comp < ncomp; comp++){
	lhs(vof, comp) += rhs(vof, comp)*a_scale;
      }
    }
  }
}

void data_ops::incr(EBAMRIVData& a_lhs, const EBAMRCellData& a_rhs, const Real a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void data_ops::incr(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale){
  const int ncomp = a_lhs.nComp();
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    BaseIVFAB<Real>& lhs = a_lhs[dit()];
    const EBCellFAB& rhs = a_rhs[dit()];

    for (VoFIterator vofit(lhs.getIVS(), lhs.getEBGraph()); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      for (int comp = 0; comp < ncomp; comp++){
	lhs(vof, comp) += rhs(vof, comp)*a_scale;
      }
    }
  }
}

void data_ops::copy(EBAMRCellData& a_dst, const EBAMRCellData& a_src){
  for (int lvl = 0; lvl < a_dst.size(); lvl++){
    if(a_src[lvl] != NULL && a_dst[lvl] != NULL){
      a_src[lvl]->copyTo(*a_dst[lvl]);
    }
  }
}

void data_ops::exponentiate(EBAMRCellData& a_lhs, const Real a_factor){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::exponentiate(*a_lhs[lvl], a_factor);
  }
}

void data_ops::exponentiate(LevelData<EBCellFAB>& a_lhs, const Real a_factor){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs         = a_lhs[dit()];
    const Box box          = a_lhs.disjointBoxLayout().get(dit());
    const EBISBox& ebisbox = lhs.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs(box);

    
    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      for (int comp = 0; comp < lhs.nComp(); comp++){
	const Real value = lhs(vof, comp);
	lhs(vof,comp) = exp(a_factor*value);
      }
    }
  }
}

void data_ops::get_max_min(Real& a_max, Real& a_min, EBAMRCellData& a_E, const int a_comp){
  for (int lvl = 0; lvl < a_E.size(); lvl++){
    data_ops::get_max_min(a_max, a_min, *a_E[lvl], a_comp);
  }
}

void data_ops::get_max_min(Real& a_max, Real& a_min, LevelData<EBCellFAB>& a_E, const int a_comp){
  EBLevelDataOps::getMaxMin(a_max, a_min, a_E, a_comp);
}

void data_ops::get_max_min(Vector<Real>& a_max, Vector<Real>& a_min, Vector<EBAMRCellData>& a_data){
  a_max.resize(a_data.size(), -1.234567E89);
  a_min.resize(a_data.size(),  1.234567E89);

  const int comp  = 0;
  const int ncomp = 1;
  
  for (int i = 0; i < a_data.size(); i++){
    CH_assert(a_data[i][0]->nComp() == ncomp);
    data_ops::get_max_min(a_max[i], a_min[i], a_data[i], comp);
  }
}

void data_ops::get_max_min_norm(Real& a_max, Real& a_min, EBAMRCellData& a_data){
  a_max = -1.234567E89;
  a_min =  1.234567E89;
  
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    Real max, min;
    data_ops::get_max_min_norm(max, min, *a_data[lvl]);

    a_max = Max(a_max, max);
    a_min = Min(a_min, min);
			       
  }
}

void data_ops::get_max_min_norm(Real& a_max, Real& a_min, LevelData<EBCellFAB>& a_data){
  a_max = -1.234567E89;
  a_min =  1.234567E89;

  const int ncomp = a_data.nComp();
  
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit){
    const Box& box         = a_data.disjointBoxLayout().get(dit());
    const EBCellFAB& data  = a_data[dit()];
    const EBISBox& ebisbox = data.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs(box);
    
    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      Real cur = 0.0;
      for (int comp = 0; comp < ncomp; comp++){
	cur += data(vof, comp)*data(vof, comp);
      }
      cur = sqrt(cur);

      a_max = Max(a_max, cur);
      a_min = Min(a_min, cur);
    }
  }

  // Communicate result
#ifdef CH_MPI
  int result;
  Real tmp = 1.;
  
  result = MPI_Allreduce(&a_max, &tmp, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("data_ops::get_max_min_norm - communication error on norm");
  }
  a_max = tmp;
  
  result = MPI_Allreduce(&a_min, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("data_ops::get_max_min_norm - communication error on norm");
  }
  a_min = tmp;
#endif
}

void data_ops::scale(MFAMRCellData& a_lhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::scale(*a_lhs[lvl], a_scale);
  }
}

void data_ops::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale){
  MFLevelDataOps::scale(a_lhs, a_scale);
}

void data_ops::scale(EBAMRIVData& a_lhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::scale(*a_lhs[lvl], a_scale);
  }
}

void data_ops::scale(EBAMRCellData& a_lhs, const Real a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::scale(*a_lhs[lvl], a_scale);
  }
}

void data_ops::scale(LevelData<EBCellFAB>& a_lhs, const Real a_scale){
  EBLevelDataOps::scale(a_lhs, a_scale);
}

void data_ops::scale(EBAMRFluxData& a_lhs, const Real a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    EBLevelDataOps::scale(*a_lhs[lvl], a_scale);
  }
}

void data_ops::scale(LevelData<BaseIVFAB<Real> >& a_lhs, const Real& a_scale){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    BaseIVFAB<Real>& lhs = a_lhs[dit()];

    for (VoFIterator vofit(lhs.getIVS(), lhs.getEBGraph()); vofit.ok(); ++vofit){
      for (int comp = 0; comp < a_lhs.nComp(); comp++){
	lhs(vofit(), comp) *= a_scale;
      }
    }
  }
}

void data_ops::divide(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_lcomp, const int a_rcomp){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::divide(*a_lhs[lvl], *a_rhs[lvl], a_lcomp, a_rcomp);
  }
}

void data_ops::divide(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int a_lcomp, const int a_rcomp){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs       = a_lhs[dit()];
    const EBCellFAB& rhs = a_rhs[dit()];

    lhs.divide(rhs, a_rcomp, a_lcomp, 1);
  }
}

void data_ops::divide_scalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::divide_scalar(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void data_ops::divide_scalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs){
  const int lcomps = a_lhs.nComp();
  const int rcomps = a_rhs.nComp();
  
  CH_assert(a_rhs.nComp() == 1);
  CH_assert(a_lhs.nComp() >= 1);

  for (int comp = 0; comp < lcomps; comp++){
    data_ops::divide(a_lhs, a_rhs, comp, 0);
  }
}

void data_ops::floor(EBAMRCellData& a_lhs, const Real a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::floor(*a_lhs[lvl], a_value);
  }
}

void data_ops::floor(LevelData<EBCellFAB>& a_lhs, const Real a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs = a_lhs[dit()];
    const Box box          = a_lhs.disjointBoxLayout().get(dit());
    const EBISBox& ebisbox = lhs.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs(lhs.getRegion());

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      for (int comp = 0; comp < a_lhs.nComp(); comp++){
	const Real value = lhs(vof, comp);
	lhs(vof, comp) = Max(value, a_value);
      }
    }
  }
}

void data_ops::kappa_sum(Real& a_mass, const LevelData<EBCellFAB>& a_lhs){


  Real mass = 0.;

  CH_assert(a_lhs.nComp() == 1);
  
  const int comp  = 0;
  const int ncomp = 1;
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    const Box box          = a_lhs.disjointBoxLayout().get(dit());
    const EBCellFAB& lhs   = a_lhs[dit()];
    const EBISBox& ebisbox = lhs.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs(box);
    
    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      mass += ebisbox.volFrac(vof)*lhs(vof, comp);
    }
  }

  a_mass = EBLevelDataOps::parallelSum(mass);
}

void data_ops::kappa_scale(EBAMRCellData& a_data){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    data_ops::kappa_scale(*a_data[lvl]);
  }
}

void data_ops::kappa_scale(LevelData<EBCellFAB>& a_data){
  EBLevelDataOps::kappaWeight(a_data);
}

void data_ops::kappa_scale(MFAMRCellData& a_data){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    data_ops::kappa_scale(*a_data[lvl]);
  }
}

void data_ops::kappa_scale(LevelData<MFCellFAB>& a_data){
  MFLevelDataOps::kappaWeight(a_data);
}

void data_ops::ln(EBAMRCellData& a_lhs){
  for (int lvl = 0; lvl <= a_lhs.size(); lvl++){
    data_ops::ln(*a_lhs[lvl]);
  }
}

void data_ops::ln(LevelData<EBCellFAB>& a_lhs){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs         = a_lhs[dit()];
    const Box box          = a_lhs.disjointBoxLayout().get(dit());
    const EBISBox& ebisbox = lhs.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs(box);

    
    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      for (int comp = 0; comp < lhs.nComp(); comp++){
	const Real value = lhs(vof, comp);
	lhs(vof,comp) = log(1.E-20 + value);
      }
    }
  }
}

void data_ops::multiply(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::multiply(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void data_ops::multiply(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()] *= a_rhs[dit()];
  }
}

void data_ops::multiply_scalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs){
  CH_assert(a_rhs.nComp() == 1);

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs       = a_lhs[dit()];
    const EBCellFAB& rhs = a_rhs[dit()];

    for (int comp = 0; comp < lhs.nComp(); comp++){
      lhs.mult(rhs, 0, comp, 1);
    }
  }
}

void data_ops::set_covered_value(EBAMRCellData& a_lhs, const int a_comp, const Real a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::set_covered_value(*a_lhs[lvl], a_comp, a_value);
  }
}

void data_ops::set_covered_value(LevelData<EBCellFAB>& a_lhs, const int a_comp, const Real a_value){
  EBLevelDataOps::setCoveredVal(a_lhs, a_comp, a_value);
}

void data_ops::set_value(EBAMRCellData& a_data, const Real& a_value){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    EBLevelDataOps::setVal(*a_data[lvl], a_value);
  }
}

void data_ops::set_value(EBAMRCellData& a_lhs, const Real a_value, const int a_comp){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::set_value(*a_lhs[lvl], a_value, a_comp);
  }
}

void data_ops::set_value(LevelData<EBCellFAB>& a_lhs, const Real a_value, const int a_comp){
  EBLevelDataOps::setVal(a_lhs, a_value, a_comp);
}

void data_ops::set_value(LevelData<EBCellFAB>& a_lhs, const Real a_value){
  EBLevelDataOps::setVal(a_lhs, a_value);
}

void data_ops::set_value(LevelData<EBFluxFAB>& a_lhs, const Real a_value){
  EBLevelDataOps::setVal(a_lhs, a_value);
}
  
void data_ops::set_value(LevelData<BaseIVFAB<Real> >& a_lhs, const Real a_value){
  EBLevelDataOps::setVal(a_lhs, a_value);
}

void data_ops::set_value(EBAMRFluxData& a_data, const Real& a_value){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    EBLevelDataOps::setVal(*a_data[lvl], a_value);
  }
}

void data_ops::set_value(EBAMRIVData& a_data, const Real& a_value){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    EBLevelDataOps::setVal(*a_data[lvl], a_value);
  }
}

void data_ops::set_value(MFAMRCellData& a_lhs, const Real& a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::set_value(*a_lhs[lvl], a_value);
  }
}

void data_ops::set_value(LevelData<MFCellFAB>& a_lhs, const Real& a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    MFCellFAB& lhs = a_lhs[dit()];
    lhs.setVal(a_value);
  }
}

void data_ops::set_value(MFAMRFluxData& a_lhs, const Real& a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::set_value(*a_lhs[lvl] , a_value);
  }
}

void data_ops::set_value(LevelData<MFFluxFAB>& a_lhs, const Real& a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    MFFluxFAB& lhs = a_lhs[dit()];
    lhs.setVal(a_value);
  }
}

void data_ops::set_value(MFAMRIVData& a_lhs, const Real& a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::set_value(*a_lhs[lvl] , a_value);
  }
}

void data_ops::set_value(LevelData<MFBaseIVFAB>& a_lhs, const Real& a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()].setVal(a_value);
  }
}

void data_ops::vector_length(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::vector_length(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void data_ops::vector_length(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs){
  CH_assert(a_lhs.nComp() == 1);
  CH_assert(a_rhs.nComp() == SpaceDim);

  const int comp  = 0;
  const int ncomp = a_rhs.nComp();
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs             = a_lhs[dit()];
    const EBCellFAB& rhs       = a_rhs[dit()];
    const EBISBox& ebisbox     = lhs.getEBISBox();
    const EBGraph& ebgraph     = ebisbox.getEBGraph();

    const IntVectSet ivs(a_lhs.disjointBoxLayout()[dit()]);

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      lhs(vof, comp) = 0.;

      for (int i = 0; i < ncomp; i++){
	lhs(vof, comp) += rhs(vof, i)*rhs(vof, i);
      }

      lhs(vof, comp) = sqrt(lhs(vof, comp));
    }
  }
}
