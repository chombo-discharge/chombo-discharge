/*!
  @file data_ops.cpp
  @brief Implementation of data_ops.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "data_ops.H"
#include "data_opsF_F.H"
#include "EBArith.H"
#include "EBLevelDataOps.H"
#include "MFLevelDataOps.H"
#include "CCProjectorF_F.H"
#include "EBLevelDataOpsF_F.H"

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
      EBLevelDataOps::averageCellToFace(flux_vel[dir], cell_vel, ebgraph, box, 1, dir, a_domain, dir, 0);
    }
  }
}

void data_ops::average_cell_to_face_allcomps(EBAMRFluxData& a_face_data,
					     const EBAMRCellData& a_cell_data,
					     const Vector<ProblemDomain>& a_domains){
  for (int lvl = 0; lvl < a_face_data.size(); lvl++){
    data_ops::average_cell_to_face_allcomps(*a_face_data[lvl], *a_cell_data[lvl], a_domains[lvl]);
  }
}

void data_ops::average_cell_to_face_allcomps(LevelData<EBFluxFAB>&       a_facedata,
					     const LevelData<EBCellFAB>& a_celldata,
					     const ProblemDomain&        a_domain){

  CH_assert(a_facedata.nComp() == a_celldata.nComp());

#if 1 // New code, should perform 5x faster than the crap below. 
  const int ncomp = a_facedata.nComp();
  for (DataIterator dit = a_facedata.dataIterator(); dit.ok(); ++dit){
    const EBCellFAB& celldata = a_celldata[dit()];
    const EBISBox& ebisbox    = celldata.getEBISBox();
    const EBGraph& ebgraph    = ebisbox.getEBGraph();
    const Box box             = a_celldata.disjointBoxLayout().get(dit());
    
    for (int dir = 0; dir < SpaceDim; dir++){
      for (int icomp = 0; icomp < ncomp; icomp++){
	EBLevelDataOps::averageCellToFace(a_facedata[dit()][dir],
					  celldata,
					  ebgraph,
					  box,
					  1,
					  dir,
					  a_domain,
					  icomp,
					  icomp);

      }
    }
  }

#else // Original, non-performant code
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
#endif
}

void data_ops::average_face_to_cell(EBAMRCellData&               a_celldata,
				    const EBAMRFluxData&         a_facedata,
				    const Vector<ProblemDomain>& a_domains){
  for (int lvl = 0; lvl < a_facedata.size(); lvl++){
    data_ops::average_face_to_cell(*a_celldata[lvl], *a_facedata[lvl], a_domains[lvl]);
  }
}

void data_ops::average_face_to_cell(LevelData<EBCellFAB>&       a_celldata,
				    const LevelData<EBFluxFAB>& a_fluxdata,
				    const ProblemDomain&        a_domain){

  const int nc = a_celldata.nComp();
  CH_assert(a_fluxdata.nComp() == nc);

  for (DataIterator dit = a_celldata.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& celldata       = a_celldata[dit()];
    const EBFluxFAB& fluxdata = a_fluxdata[dit()];
    const Box& box            = a_celldata.disjointBoxLayout().get(dit());

    // Irregular cells
    const EBISBox& ebisbox = celldata.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);

    // Reset
    celldata.setVal(0.0);

    // Do regular cells; all components in each direction
    for (int dir = 0; dir < SpaceDim; dir++){
      BaseFab<Real>& cellreg       = celldata.getSingleValuedFAB();
      //      const BaseFab<Real>& facereg = fluxdata[dir].getSingleValuedFAB();
      auto& facereg = fluxdata[dir].getSingleValuedFAB();
      FORT_AVERAGE_FACE_TO_CELL(CHF_FRA(cellreg),
				CHF_CONST_FRA(facereg),
				CHF_CONST_INT(dir),
				CHF_CONST_INT(nc),
				CHF_BOX(box));
    }

    // Reset irregular cells
    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      for (int ic = 0; ic < nc; ic++){
	celldata(vof, ic) = 0.0;
      }
    }
  }
}

void data_ops::dot_prod(MFAMRCellData& a_result, const MFAMRCellData& a_data1, const MFAMRCellData& a_data2){
  for (int lvl = 0; lvl < a_result.size(); lvl++){
    data_ops::dot_prod(*a_result[lvl], *a_data1[lvl], *a_data2[lvl]);
  }
}

void data_ops::dot_prod(LevelData<MFCellFAB>& a_result, const LevelData<MFCellFAB>& a_data1, const LevelData<MFCellFAB>& a_data2){
  const int nc = a_data1.nComp();

  CH_assert(a_data2.nComp() == nc);
  CH_assert(a_result.nComp() == 1);

  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit){
    MFCellFAB& result      = a_result[dit()];
    const MFCellFAB& data1 = a_data1[dit()];
    const MFCellFAB& data2 = a_data2[dit()];
    const Box& box = a_result.disjointBoxLayout().get(dit());

    for (int i = 0; i < result.numPhases(); i++){
      EBCellFAB& result_phase      = result.getPhase(i);
      const EBCellFAB& data1_phase = data1.getPhase(i);
      const EBCellFAB& data2_phase = data2.getPhase(i);

      data_ops::dot_prod(result_phase, data1_phase, data2_phase, box);
    }
  }
}

void data_ops::dot_prod(EBAMRCellData& a_result, const EBAMRCellData& a_data1, const EBAMRCellData& a_data2){
  for (int lvl = 0; lvl < a_result.size(); lvl++){
    data_ops::dot_prod(*a_result[lvl], *a_data1[lvl], *a_data2[lvl]);
  }
}

void data_ops::dot_prod(LevelData<EBCellFAB>& a_result,
			const LevelData<EBCellFAB>& a_data1,
			const LevelData<EBCellFAB>& a_data2){
  const int nc = a_data1.nComp();

  CH_assert(a_data2.nComp() == nc);
  CH_assert(a_result.nComp() == 1);

  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& result      = a_result[dit()];
    const EBCellFAB& data1 = a_data1[dit()];
    const EBCellFAB& data2 = a_data2[dit()];
    const Box& box         = a_result.disjointBoxLayout().get(dit());
    data_ops::dot_prod(result, data1, data2, box);

#if 0
    // Regular cells
    BaseFab<Real>& result_reg      = result.getSingleValuedFAB();
    const BaseFab<Real>& data1_reg = data1.getSingleValuedFAB();
    const BaseFab<Real>& data2_reg = data2.getSingleValuedFAB();
    FORT_DOT_PRODUCT(CHF_FRA1(result_reg, 0),
		     CHF_CONST_FRA(data1_reg),
		     CHF_CONST_FRA(data2_reg),
		     CHF_CONST_INT(nc),
		     CHF_BOX(box));
    

    // Irregular cells
    const EBISBox& ebisbox = result.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);
    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      result(vof, 0) = 0.0;
      for (int comp = 0; comp < nc; comp++){
	result(vof, 0) += data1(vof, comp)*data2(vof,comp);
      }
    }
#endif
  }
}

void data_ops::dot_prod(EBCellFAB& a_result, const EBCellFAB& a_data1, const EBCellFAB& a_data2, const Box& a_box){
  const int nc = a_data1.nComp();
  BaseFab<Real>& result_reg      = a_result.getSingleValuedFAB();
  const BaseFab<Real>& data1_reg = a_data1.getSingleValuedFAB();
  const BaseFab<Real>& data2_reg = a_data2.getSingleValuedFAB();
  FORT_DOT_PRODUCT(CHF_FRA1(result_reg, 0),
		   CHF_CONST_FRA(data1_reg),
		   CHF_CONST_FRA(data2_reg),
		   CHF_CONST_INT(nc),
		   CHF_BOX(a_box));
    

  // Irregular cells
  const EBISBox& ebisbox = a_result.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const IntVectSet& ivs  = ebisbox.getIrregIVS(a_box);
  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    a_result(vof, 0) = 0.0;
    for (int comp = 0; comp < nc; comp++){
      a_result(vof, 0) += a_data1(vof, comp)*a_data2(vof,comp);
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

void data_ops::incr(EBAMRIFData& a_lhs, const EBAMRIFData& a_rhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void data_ops::incr(LevelData<DomainFluxIFFAB>& a_lhs, const LevelData<DomainFluxIFFAB>& a_rhs, const Real& a_scale){
  CH_assert(a_lhs.nComp() == a_rhs.nComp());

  const int ncomp = a_lhs.nComp();
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    DomainFluxIFFAB& lhs       = a_lhs[dit()];
    const DomainFluxIFFAB& rhs = a_rhs[dit()];

    for (int dir = 0; dir < SpaceDim; dir++){
      for (SideIterator sit; sit.ok(); ++sit){
	BaseIFFAB<Real>& cur_lhs       = lhs(dir, sit());
	const BaseIFFAB<Real>& cur_rhs = rhs(dir, sit());

	const IntVectSet& ivs  = cur_lhs.getIVS();
	const EBGraph& ebgraph = cur_lhs.getEBGraph();
	const FaceStop::WhichFaces stop_crit = FaceStop::SurroundingWithBoundary;

	for (FaceIterator faceit(ivs, ebgraph, dir, stop_crit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();

	  for (int comp = 0; comp < ncomp; comp++){
	    cur_lhs(face, comp) += cur_rhs(face, comp)*a_scale;
	  }
	}
	
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

void data_ops::copy(MFAMRCellData& a_dst, const MFAMRCellData& a_src){
  for (int lvl = 0; lvl < a_dst.size(); lvl++){
    if(a_src[lvl] != NULL && a_dst[lvl] != NULL){
      a_src[lvl]->localCopyTo(*a_dst[lvl]);
    }
  }
}

void data_ops::copy(EBAMRCellData& a_dst, const EBAMRCellData& a_src){
  for (int lvl = 0; lvl < a_dst.size(); lvl++){
    if(a_src[lvl] != NULL && a_dst[lvl] != NULL){
      a_src[lvl]->localCopyTo(*a_dst[lvl]);
    }
  }
}

void data_ops::copy(EBAMRIVData& a_dst, const EBAMRIVData& a_src){
  for (int lvl = 0; lvl < a_dst.size(); lvl++){
    if(a_src[lvl] != NULL && a_dst[lvl] != NULL){
      a_src[lvl]->localCopyTo(*a_dst[lvl]);
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
  a_max = -1.E99;
  a_min =  1.E99;
  for (int lvl = 0; lvl < a_E.size(); lvl++){
    Real lvl_max = -1.E99;
    Real lvl_min =  1.E99;
    data_ops::get_max_min(lvl_max, lvl_min, *a_E[lvl], a_comp);

    a_max = Max(a_max, lvl_max);
    a_min = Min(a_min, lvl_min);
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

#if 1 // Optimized code
    EBCellFAB covered_mask(ebisbox, box, 1);
    covered_mask.setVal(1.0);
    covered_mask.setCoveredCellVal(-1.0, 0);

    const BaseFab<Real>& mask = covered_mask.getSingleValuedFAB();
    
    // Maybe this breaks because we should pass a covered flag into the routine
    const BaseFab<Real>& data_reg = data.getSingleValuedFAB();
    FORT_MAX_MIN_NORM(CHF_REAL(a_max),
		      CHF_REAL(a_min),
		      CHF_CONST_FRA(data_reg),
      		      CHF_CONST_FRA1(mask, 0),
		      CHF_CONST_INT(ncomp),
		      CHF_BOX(box));

    // Irregular and multivalued cells
    for (VoFIterator vofit(ebisbox.getIrregIVS(box), ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      Real cur = 0.0;
      for (int comp = 0; comp < ncomp; comp++){
	cur += data(vof, comp)*data(vof, comp);
      }
      cur = sqrt(cur);

      a_max = Max(a_max, cur);
      a_min = Min(a_min, cur);
    }
    
#else // Original code
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
#endif
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

void data_ops::get_max_min_norm(Real& a_max, Real& a_min, EBAMRIVData& a_data){
  a_max = -1.234567E89;
  a_min =  1.234567E89;
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    Real max, min;
    data_ops::get_max_min_norm(max, min, *a_data[lvl]);

    a_max = Max(a_max, max);
    a_min = Min(a_min, min);
  }
}

void data_ops::get_max_min_norm(Real& a_max, Real& a_min, LevelData<BaseIVFAB<Real> >& a_data){
  a_max = -1.234567E89;
  a_min =  1.234567E89;

  const int ncomp = a_data.nComp();
  
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit){
    const Box& box              = a_data.disjointBoxLayout().get(dit());
    const BaseIVFAB<Real>& data = a_data[dit()];
    const IntVectSet ivs(box);

    // Irregular and multivalued cells
    for (VoFIterator vofit(data.getIVS(), data.getEBGraph()); vofit.ok(); ++vofit){
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
    const Box box = lhs.getRegion();
    const EBISBox& ebisbox = lhs.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs(lhs.getRegion());

    const int ncomp = a_lhs.nComp();

#if 1 // Optimized code. 
    // Regular cells. This also does ghost cells
    BaseFab<Real>& lhs_reg = lhs.getSingleValuedFAB();
    FORT_FLOOR(CHF_FRA(lhs_reg),
	       CHF_CONST_INT(ncomp),
	       CHF_CONST_REAL(a_value),
	       CHF_BOX(box));

    // Irregular and multivalued cells
    for (VoFIterator vofit(ebisbox.getIrregIVS(box), ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      for (int comp = 0; comp < a_lhs.nComp(); comp++){
	const Real value = lhs(vof, comp);
	lhs(vof, comp) = Max(value, a_value);
      }
    }
#else // Other code
    // Irregular and multivalued cells
    for (VoFIterator vofit(IntVectSet(box), ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      for (int comp = 0; comp < a_lhs.nComp(); comp++){
	const Real value = lhs(vof, comp);
	lhs(vof, comp) = Max(value, a_value);
      }
    }
#endif
  }
}

void data_ops::floor(EBAMRIVData& a_lhs, const Real a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::floor(*a_lhs[lvl], a_value);
  }
}

void data_ops::floor(LevelData<BaseIVFAB<Real> >& a_lhs, const Real a_value){
  const int ncomp = a_lhs.nComp();
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    BaseIVFAB<Real>& lhs   = a_lhs[dit()];
    for (VoFIterator vofit(lhs.getIVS(), lhs.getEBGraph()); vofit.ok(); ++vofit){
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

void data_ops::laplacian(EBAMRCellData& a_lapl, const EBAMRCellData& a_data){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    data_ops::laplacian(*a_lapl[lvl], *a_data[lvl]);
  }
}

void data_ops::laplacian(LevelData<EBCellFAB>& a_lapl, const LevelData<EBCellFAB>& a_data){
  const int ncomp = a_data.nComp();
  
  for (DataIterator dit = a_lapl.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lapl        = a_lapl[dit()];
    const EBCellFAB& data  = a_data[dit()];
    const EBISBox& ebisbox = data.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const Box box          = a_lapl.disjointBoxLayout().get(dit());
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);

    BaseFab<Real>& lapl_fab       = lapl.getSingleValuedFAB();
    const BaseFab<Real>& data_fab = data.getSingleValuedFAB();

    // Regular stuff
    for (int comp = 0; comp < ncomp; comp++){
      FORT_LAPLACIAN(CHF_FRA1(lapl_fab, comp),
		     CHF_CONST_FRA1(data_fab, comp),
		     CHF_BOX(box));
    }
    

    // Irregular and multi-valued
    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();


      for (int comp = 0; comp < ncomp; comp++){
	lapl(vof, comp) = 0.0;
	for (int dir = 0; dir < SpaceDim; dir++){
	  Real value = 0.0;
	  VoFStencil sten;

	  // Get stencil
	  EBArith::getSecondDerivStencil(sten, vof, ebisbox, dir, 1.0);

	  // Apply it
	  for (int i = 0; i < sten.size(); i++){
	    const VolIndex& ivof = sten.vof(i);
	    const Real iweight   = sten.weight(i);
	    lapl(vof, comp) += iweight*data(ivof,comp);
	  }
	}
      }
    }
  }
}

void data_ops::gen_laplacian(EBAMRCellData& a_lapl, const EBAMRCellData& a_data){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    data_ops::gen_laplacian(*a_lapl[lvl], *a_data[lvl]);
  }
}

void data_ops::gen_laplacian(LevelData<EBCellFAB>& a_lapl, const LevelData<EBCellFAB>& a_data){
  const int ncomp = a_data.nComp();
  
  for (DataIterator dit = a_lapl.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lapl        = a_lapl[dit()];
    const EBCellFAB& data  = a_data[dit()];
    const EBISBox& ebisbox = data.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const Box box          = a_lapl.disjointBoxLayout().get(dit());
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);

    BaseFab<Real>& lapl_fab       = lapl.getSingleValuedFAB();
    const BaseFab<Real>& data_fab = data.getSingleValuedFAB();

    // Regular stuff
    for (int comp = 0; comp < ncomp; comp++){
      FORT_GEN_LAPLACIAN(CHF_FRA1(lapl_fab, comp),
		     CHF_CONST_FRA1(data_fab, comp),
		     CHF_BOX(box));
    }
    

    // Irregular and multi-valued
    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();


      for (int comp = 0; comp < ncomp; comp++){
	lapl(vof, comp) = 0.0;
	for (int dir = 0; dir < SpaceDim; dir++){
	  Real value = 0.0;
	  VoFStencil sten;

	  // Get stencil
	  EBArith::getSecondDerivStencil(sten, vof, ebisbox, dir, 1.0);

	  // Apply it
	  for (int i = 0; i < sten.size(); i++){
	    const VolIndex& ivof = sten.vof(i);
	    const Real iweight   = sten.weight(i);
	    lapl(vof, comp) += iweight*data(ivof,comp);
	  }
	}

#if 1 // Debug
	lapl(vof,comp) = 0.0;
#endif
      }
    }
  }
}

void data_ops::flash_error(EBAMRCellData& a_lapl, const EBAMRCellData& a_data, const Real a_eps){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    data_ops::flash_error(*a_lapl[lvl], *a_data[lvl], a_eps);
  }
}

void data_ops::flash_error(LevelData<EBCellFAB>& a_lapl, const LevelData<EBCellFAB>& a_data, const Real a_eps){
  const int ncomp = a_data.nComp();
  const DisjointBoxLayout& dbl = a_lapl.disjointBoxLayout();
  const ProblemDomain domain   = dbl.physDomain();

  for (DataIterator dit = a_lapl.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lapl        = a_lapl[dit()];
    const EBCellFAB& data  = a_data[dit()];
    const EBISBox& ebisbox = data.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const Box box          = dbl.get(dit());
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);

    BaseFab<Real>& lapl_fab       = lapl.getSingleValuedFAB();
    const BaseFab<Real>& data_fab = data.getSingleValuedFAB();

    // Regular stuff
    for (int comp = 0; comp < ncomp; comp++){
      FORT_FLASH_ERROR(CHF_FRA1(lapl_fab, comp),
		       CHF_CONST_FRA1(data_fab, comp),
		       CHF_CONST_REAL(a_eps),
		       CHF_BOX(box));
    }

    // Can't trust stuff on the sides
    IntVectSet bndry_ivs = ebisbox.getIrregIVS(dbl.get(dit()));
    for (int dir = 0; dir < SpaceDim; dir++){
      Box lo_box, hi_box;
      int has_lo, has_hi;

      EBArith::loHi(lo_box, has_lo, hi_box, has_hi, domain, box, dir);

      if(has_lo) bndry_ivs |= IntVectSet(lo_box);
      if(has_hi) bndry_ivs |= IntVectSet(hi_box);
    }
      
    // Compute stencils for boundary cells
    for (VoFIterator vofit(bndry_ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      for (int comp = 0; comp < ncomp; comp++){
	lapl(vof, comp) = 0.;
      }
    }

    // Set covered to zero
    for (int comp = 0; comp < ncomp; comp++){
      lapl.setCoveredCellVal(0.0, comp);
    }
    

    // Irregular and multi-valued. Set these to zero for now
    IntVectSet irreg = ivs;
    //    irreg.grow(2);
    irreg.grow(2);
    irreg &= box;
    for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      
      for (int comp = 0; comp < ncomp; comp++){
	lapl(vof, comp) = 0.0;
      }
    }
  }
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

void data_ops::multiply(EBAMRFluxData& a_lhs, const EBAMRFluxData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::multiply(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void data_ops::multiply(LevelData<EBFluxFAB>& a_lhs, const LevelData<EBFluxFAB>& a_rhs){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()] *= a_rhs[dit()];
  }
}

void data_ops::multiply_scalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::multiply_scalar(*a_lhs[lvl], *a_rhs[lvl]);
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

void data_ops::multiply_scalar(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::multiply_scalar(*a_lhs[lvl], *a_rhs[lvl]);
  }
}
  
void data_ops::multiply_scalar(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs){
  CH_assert(a_rhs.nComp() == 1);

  const int ncomp = a_lhs.nComp();
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    BaseIVFAB<Real>& lhs       = a_lhs[dit()];
    const BaseIVFAB<Real>& rhs = a_rhs[dit()];
    const EBGraph& ebgraph     = lhs.getEBGraph();
    const IntVectSet& ivs      = lhs.getIVS();

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      for (int comp = 0; comp < ncomp; comp++){
	lhs(vof, comp) *= rhs(vof, 0);
      }
    }
  }
}

void data_ops::norm(Real& a_norm, const LevelData<EBCellFAB>& a_data, const ProblemDomain& a_domain, const int a_p){
  Real volume;

  a_norm = EBLevelDataOps::kappaNorm(volume, a_data, EBLEVELDATAOPS_ALLVOFS, a_domain, a_p);
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

void data_ops::set_value(EBAMRIFData& a_lhs, const Real a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::set_value(*a_lhs[lvl], a_value);
  }
}

void data_ops::set_value(LevelData<DomainFluxIFFAB>& a_lhs, const Real a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    DomainFluxIFFAB& lhs = a_lhs[dit()];
    for (SideIterator sit; sit.ok(); ++sit){
      for (int dir = 0; dir < SpaceDim; dir++){
	BaseIFFAB<Real>& fab = lhs(dir, sit());
	fab.setVal(a_value);
      }
    }
  }
}

void data_ops::sum(Real& a_value){

#ifdef CH_MPI
  Real cur = a_value;
  Real tmp = a_value;
  int result = MPI_Allreduce(&cur, &tmp, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Abort("data_ops::sum - communication error on Allreduce");
  }

  a_value = tmp;
#endif
}

void data_ops::square_root(EBAMRFluxData& a_lhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::square_root(*a_lhs[lvl]);
  }
}

void data_ops::square_root(LevelData<EBFluxFAB>& a_lhs){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){

    const Box& box         = a_lhs.disjointBoxLayout().get(dit());

    
    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& lhs         = a_lhs[dit()][dir];
      BaseFab<Real>& lhs_reg = lhs.getSingleValuedFAB();

      // Face centered box
      Box facebox = box;
      facebox.surroundingNodes(dir);

      // All comps
      for (int comp = 0; comp < lhs.nComp(); comp++){
	FORT_SQUARE_ROOT(CHF_FRA1(lhs_reg, comp),
			 CHF_CONST_INT(dir),
			 CHF_BOX(facebox));

	const FaceStop::WhichFaces stop_crit = FaceStop::SurroundingWithBoundary;
	const EBISBox& ebisbox = lhs.getEBISBox();
	const EBGraph& ebgraph = ebisbox.getEBGraph();
	const IntVectSet& ivs  = ebisbox.getIrregIVS(box);
	
	for (FaceIterator faceit(ivs, ebgraph, dir, stop_crit); faceit.ok(); ++faceit){
	  const FaceIndex& face  = faceit();
	  lhs(face, comp) = sqrt(lhs(face, comp));
	}
      }
    }
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
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs             = a_lhs[dit()];
    const Box& box             = a_lhs.disjointBoxLayout().get(dit());
    const EBCellFAB& rhs       = a_rhs[dit()];

    data_ops::vector_length(lhs, rhs, box);
  }
}

void data_ops::vector_length(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& a_box){
  CH_assert(a_lhs.nComp() == 1);
  CH_assert(a_rhs.nComp() == SpaceDim);

  const int comp = 0;
  const int ncomp = SpaceDim;

  const EBISBox& ebisbox = a_lhs.getEBISBox();

#if 1 // Optimized code
  // Mask for skipping computation on covered cells
  EBCellFAB covered_mask(ebisbox, a_box, 1);
  covered_mask.setVal(1.0);
  covered_mask.setCoveredCellVal(-1.0, 0);
  
  // Regular cells
  BaseFab<Real>& lhs_reg       = a_lhs.getSingleValuedFAB();
  const BaseFab<Real>& rhs_reg = a_rhs.getSingleValuedFAB();
  const BaseFab<Real>& mask    = covered_mask.getSingleValuedFAB();

  FORT_VECTOR_LENGTH(CHF_FRA1(lhs_reg, comp),
		     CHF_CONST_FRA(rhs_reg),
		     CHF_CONST_FRA1(mask,0),
		     CHF_BOX(a_box));


  // Irregular cells and multivalued cells
  for (VoFIterator vofit(ebisbox.getIrregIVS(a_box), ebisbox.getEBGraph()); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();

    a_lhs(vof, comp) = 0.;

    for (int i = 0; i < ncomp; i++){
      a_lhs(vof, comp) += a_rhs(vof, i)*a_rhs(vof, i);
    }

    a_lhs(vof, comp) = sqrt(a_lhs(vof, comp));
  }
#else // Other code
  // Irregular cells and multivalued cells
  for (VoFIterator vofit(IntVectSet(a_box), ebisbox.getEBGraph()); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();

    a_lhs(vof, comp) = 0.;

    for (int i = 0; i < ncomp; i++){
      a_lhs(vof, comp) += a_rhs(vof, i)*a_rhs(vof, i);
    }

    a_lhs(vof, comp) = sqrt(a_lhs(vof, comp));
  }
#endif
}

void data_ops::vector_length2(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::vector_length2(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void data_ops::vector_length2(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs){
  CH_assert(a_lhs.nComp() == 1);
  CH_assert(a_rhs.nComp() == SpaceDim);
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs             = a_lhs[dit()];
    const Box& box             = a_lhs.disjointBoxLayout().get(dit());
    const EBCellFAB& rhs       = a_rhs[dit()];

    data_ops::vector_length2(lhs, rhs, box);
  }
}

void data_ops::vector_length2(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& a_box){
  CH_assert(a_lhs.nComp() == 1);
  CH_assert(a_rhs.nComp() == SpaceDim);

  const int comp = 0;
  const int ncomp = SpaceDim;

  const EBISBox& ebisbox = a_lhs.getEBISBox();

  // Mask for skipping computation on covered cells
  EBCellFAB covered_mask(ebisbox, a_box, 1);
  covered_mask.setVal(1.0);
  covered_mask.setCoveredCellVal(-1.0, 0);
  
  // Regular cells
  BaseFab<Real>& lhs_reg       = a_lhs.getSingleValuedFAB();
  const BaseFab<Real>& rhs_reg = a_rhs.getSingleValuedFAB();
  const BaseFab<Real>& mask    = covered_mask.getSingleValuedFAB();

  FORT_VECTOR_LENGTH2(CHF_FRA1(lhs_reg, comp),
		      CHF_CONST_FRA(rhs_reg),
		      CHF_CONST_FRA1(mask,0),
		      CHF_BOX(a_box));


  // Irregular cells and multivalued cells
  for (VoFIterator vofit(ebisbox.getIrregIVS(a_box), ebisbox.getEBGraph()); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();

    a_lhs(vof, comp) = 0.;

    for (int i = 0; i < ncomp; i++){
      a_lhs(vof, comp) += a_rhs(vof, i)*a_rhs(vof, i);
    }
  }
}
