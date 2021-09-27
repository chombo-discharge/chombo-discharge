/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DataOps.cpp
  @brief  Implementation of CD_DataOps.H
  @author Robert Marskar
*/

#include <EBArith.H>
#include <EBLevelDataOps.H>
#include <MFLevelDataOps.H>
#include <CCProjectorF_F.H>
#include <EBLevelDataOpsF_F.H>
#include <PolyGeom.H>

// Our includes
#include <CD_DataOps.H>
#include <CD_DataOpsF_F.H>
#include <CD_NamespaceHeader.H>

void DataOps::averageCellVectorToFaceScalar(EBAMRFluxData&               a_facedata,
					    const EBAMRCellData&         a_celldata,
					    const Vector<ProblemDomain>& a_domains){

  for (int lvl = 0; lvl < a_facedata.size(); lvl++){

    CH_assert(a_facedata[lvl]->nComp() == 1       );
    CH_assert(a_celldata[lvl]->nComp() == SpaceDim);
    
    DataOps::averageCellVectorToFaceScalar(*a_facedata[lvl], *a_celldata[lvl], a_domains[lvl]);
  }
}

void DataOps::averageCellVectorToFaceScalar(LevelData<EBFluxFAB>&       a_facedata,
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

void DataOps::averageCellToFace(EBAMRFluxData& a_face_data,
				const EBAMRCellData& a_cell_data,
				const Vector<ProblemDomain>& a_domains){
  for (int lvl = 0; lvl < a_face_data.size(); lvl++){
    DataOps::averageCellToFace(*a_face_data[lvl], *a_cell_data[lvl], a_domains[lvl]);
  }
}

void DataOps::averageCellToFace(LevelData<EBFluxFAB>&       a_facedata,
				const LevelData<EBCellFAB>& a_celldata,
				const ProblemDomain&        a_domain){

  CH_assert(a_facedata.nComp() == a_celldata.nComp());

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
					  0,
					  dir,
					  a_domain,
					  icomp,
					  icomp);
      }
    }
  }
}

void DataOps::averageFaceToCell(EBAMRCellData&               a_celldata,
				const EBAMRFluxData&         a_facedata,
				const Vector<ProblemDomain>& a_domains){
  for (int lvl = 0; lvl < a_facedata.size(); lvl++){
    DataOps::averageFaceToCell(*a_celldata[lvl], *a_facedata[lvl], a_domains[lvl]);
  }
}

void DataOps::averageFaceToCell(LevelData<EBCellFAB>&       a_celldata,
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

      int nfaces = 0;
      for (int dir = 0; dir < SpaceDim; dir++){
	for (SideIterator sit; sit.ok(); ++sit){
	  Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit());

	  nfaces += faces.size();

	  for (int iface = 0; iface < faces.size(); iface++){
	    for (int ic = 0; ic < nc; ic++){
	      celldata(vof, ic) += fluxdata[dir](faces[iface],ic);
	    }
	  }
	}
      }
      for (int ic = 0; ic < nc; ic++){
	celldata(vof,ic) *= 1./(nfaces);
      }
    }
  }
}

void DataOps::dotProduct(MFAMRCellData& a_result, const MFAMRCellData& a_data1, const MFAMRCellData& a_data2){
  for (int lvl = 0; lvl < a_result.size(); lvl++){
    DataOps::dotProduct(*a_result[lvl], *a_data1[lvl], *a_data2[lvl]);
  }
}

void DataOps::dotProduct(LevelData<MFCellFAB>& a_result, const LevelData<MFCellFAB>& a_data1, const LevelData<MFCellFAB>& a_data2){
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

      DataOps::dotProduct(result_phase, data1_phase, data2_phase, box);
    }
  }
}

void DataOps::dotProduct(EBAMRCellData& a_result, const EBAMRCellData& a_data1, const EBAMRCellData& a_data2){
  for (int lvl = 0; lvl < a_result.size(); lvl++){
    DataOps::dotProduct(*a_result[lvl], *a_data1[lvl], *a_data2[lvl]);
  }
}

void DataOps::dotProduct(LevelData<EBCellFAB>& a_result,
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
    DataOps::dotProduct(result, data1, data2, box);

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

void DataOps::dotProduct(EBCellFAB& a_result, const EBCellFAB& a_data1, const EBCellFAB& a_data2, const Box& a_box){
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

void DataOps::incr(MFAMRCellData& a_lhs, const MFAMRCellData& a_rhs, const Real a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}
  
void DataOps::incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const Real a_scale){
  MFLevelDataOps::incr(a_lhs, a_rhs, a_scale);
}

void DataOps::incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void DataOps::incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real& a_scale){
  EBLevelDataOps::incr(a_lhs, a_rhs, a_scale);
}

void DataOps::plus(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_srcComp, const int a_dstComp, const int a_numComp){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::plus(*a_lhs[lvl], *a_rhs[lvl], a_srcComp, a_dstComp, a_numComp);
  }
}

void DataOps::plus(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int a_srcComp, const int a_dstComp, const int a_numComp){
  for (DataIterator dit = a_lhs.disjointBoxLayout().dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs       = a_lhs[dit()];
    const EBCellFAB& rhs = a_rhs[dit()];

    lhs.plus(rhs, a_srcComp, a_dstComp, a_numComp);
  }
}

void DataOps::incr(EBAMRFluxData& a_lhs, const EBAMRFluxData& a_rhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void DataOps::incr(LevelData<EBFluxFAB>& a_lhs, const LevelData<EBFluxFAB>& a_rhs, const Real& a_scale){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    DataOps::incr(a_lhs[dit()], a_rhs[dit()], a_scale);
  }
}

void DataOps::incr(EBFluxFAB& a_lhs, const EBFluxFAB& a_rhs, const Real& a_scale){

  EBFluxFAB rhsClone;
  rhsClone.clone(a_rhs);
  
  for (int dir = 0; dir < SpaceDim; dir++){
    EBFaceFAB& lhs       = a_lhs[dir];
    
    rhsClone[dir] *= a_scale;
    lhs += rhsClone[dir];
  }
}

void DataOps::incr(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs, const Real& a_scale){

  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void DataOps::incr(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs, const Real& a_scale){
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

void DataOps::incr(EBAMRIFData& a_lhs, const EBAMRIFData& a_rhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void DataOps::incr(LevelData<DomainFluxIFFAB>& a_lhs, const LevelData<DomainFluxIFFAB>& a_rhs, const Real& a_scale){
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

void DataOps::incr(EBAMRCellData& a_lhs, const EBAMRIVData& a_rhs, const Real a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void DataOps::incr(LevelData<EBCellFAB>& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs, const Real a_scale){
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

void DataOps::incr(EBAMRIVData& a_lhs, const EBAMRCellData& a_rhs, const Real a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void DataOps::incr(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale){
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

void DataOps::copy(MFAMRCellData& a_dst, const MFAMRCellData& a_src){
  for (int lvl = 0; lvl < a_dst.size(); lvl++){
    if(a_src[lvl] != NULL && a_dst[lvl] != NULL){
      a_src[lvl]->localCopyTo(*a_dst[lvl]);
    }
  }
}

void DataOps::copy(EBAMRCellData& a_dst, const EBAMRCellData& a_src){
  for (int lvl = 0; lvl < a_dst.size(); lvl++){
    if(a_src[lvl] != NULL && a_dst[lvl] != NULL){
      a_src[lvl]->localCopyTo(*a_dst[lvl]);
    }
  }
}

void DataOps::copy(EBAMRIVData& a_dst, const EBAMRIVData& a_src){
  for (int lvl = 0; lvl < a_dst.size(); lvl++){
    if(a_src[lvl] != NULL && a_dst[lvl] != NULL){
      a_src[lvl]->localCopyTo(*a_dst[lvl]);
    }
  }
}

void DataOps::exponentiate(EBAMRCellData& a_lhs, const Real a_factor){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::exponentiate(*a_lhs[lvl], a_factor);
  }
}

void DataOps::exponentiate(LevelData<EBCellFAB>& a_lhs, const Real a_factor){
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

void DataOps::getMaxMin(Real& a_max, Real& a_min, EBAMRCellData& a_E, const int a_comp){
  a_max = -1.E99;
  a_min =  1.E99;
  for (int lvl = 0; lvl < a_E.size(); lvl++){
    Real lvl_max = -1.E99;
    Real lvl_min =  1.E99;
    DataOps::getMaxMin(lvl_max, lvl_min, *a_E[lvl], a_comp);

    a_max = Max(a_max, lvl_max);
    a_min = Min(a_min, lvl_min);
  }
}

void DataOps::getMaxMin(Real& a_max, Real& a_min, LevelData<EBCellFAB>& a_E, const int a_comp){
  EBLevelDataOps::getMaxMin(a_max, a_min, a_E, a_comp);
}

void DataOps::getMaxMin(Vector<Real>& a_max, Vector<Real>& a_min, Vector<EBAMRCellData>& a_data){
  a_max.resize(a_data.size(), -1.234567E89);
  a_min.resize(a_data.size(),  1.234567E89);

  const int comp  = 0;
  const int ncomp = 1;
  
  for (int i = 0; i < a_data.size(); i++){
    CH_assert(a_data[i][0]->nComp() == ncomp);
    DataOps::getMaxMin(a_max[i], a_min[i], a_data[i], comp);
  }
}

void DataOps::getMaxMinNorm(Real& a_max, Real& a_min, EBAMRCellData& a_data){
  a_max = -1.234567E89;
  a_min =  1.234567E89;
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    Real max, min;
    DataOps::getMaxMinNorm(max, min, *a_data[lvl]);

    a_max = Max(a_max, max);
    a_min = Min(a_min, min);
  }
}

void DataOps::getMaxMinNorm(Real& a_max, Real& a_min, LevelData<EBCellFAB>& a_data){
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
    MayDay::Error("DataOps::getMaxMinNorm - communication error on norm");
  }
  a_max = tmp;
  
  result = MPI_Allreduce(&a_min, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("DataOps::getMaxMinNorm - communication error on norm");
  }
  a_min = tmp;
#endif
}

void DataOps::getMaxMinNorm(Real& a_max, Real& a_min, EBAMRIVData& a_data){
  a_max = -1.234567E89;
  a_min =  1.234567E89;
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    Real max, min;
    DataOps::getMaxMinNorm(max, min, *a_data[lvl]);

    a_max = Max(a_max, max);
    a_min = Min(a_min, min);
  }
}

void DataOps::getMaxMinNorm(Real& a_max, Real& a_min, LevelData<BaseIVFAB<Real> >& a_data){
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
    MayDay::Error("DataOps::getMaxMinNorm - communication error on norm");
  }
  a_max = tmp;
  
  result = MPI_Allreduce(&a_min, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("DataOps::getMaxMinNorm - communication error on norm");
  }
  a_min = tmp;
#endif
}

void DataOps::invert(EBAMRFluxData& a_data){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    DataOps::invert(*a_data[lvl]);
  }
}

void DataOps::invert(LevelData<EBFluxFAB>& a_data){
  
  const int ncomp = a_data.nComp();

  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit){
    EBFluxFAB& fluxfab = a_data[dit()];
    const EBISBox& ebisbox = fluxfab.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const Box box          = a_data.disjointBoxLayout().get(dit());
    const IntVectSet irreg = ebisbox.getIrregIVS(box);

    // Do each dircetion
    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& data = a_data[dit()][dir];

      // Need a copy because regular Fortran loop inverts irregular face cells
      EBFaceFAB cpy(ebisbox, box, dir, ncomp);
      cpy.setVal(0.0);
      cpy += data;
      
      // Regular cells
      Box facebox = box;
      facebox.surroundingNodes(dir);
      BaseFab<Real>& data_fab = data.getSingleValuedFAB();
      FORT_INVERT(CHF_FRA(data_fab),
		  CHF_CONST_INT(ncomp),
		  CHF_BOX(facebox));
      

      // Irregular cells
      FaceStop::WhichFaces stop_crit = FaceStop::SurroundingWithBoundary;
      for (FaceIterator faceit(irreg, ebgraph, dir, stop_crit); faceit.ok(); ++faceit){
	const FaceIndex& face  = faceit();
	for (int comp = 0; comp < ncomp; comp++){
	  data(face, comp) = 1./cpy(face, comp);
	}
      }
    }
  }
}

void DataOps::scale(MFAMRCellData& a_lhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::scale(*a_lhs[lvl], a_scale);
  }
}

void DataOps::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale){
  MFLevelDataOps::scale(a_lhs, a_scale);
}

void DataOps::scale(EBAMRIVData& a_lhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::scale(*a_lhs[lvl], a_scale);
  }
}

void DataOps::scale(EBAMRCellData& a_lhs, const Real a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::scale(*a_lhs[lvl], a_scale);
  }
}

void DataOps::scale(LevelData<EBCellFAB>& a_lhs, const Real a_scale){
  EBLevelDataOps::scale(a_lhs, a_scale);
}

void DataOps::scale(EBAMRFluxData& a_lhs, const Real a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    EBLevelDataOps::scale(*a_lhs[lvl], a_scale);
  }
}

void DataOps::scale(LevelData<BaseIVFAB<Real> >& a_lhs, const Real& a_scale){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    BaseIVFAB<Real>& lhs = a_lhs[dit()];

    for (VoFIterator vofit(lhs.getIVS(), lhs.getEBGraph()); vofit.ok(); ++vofit){
      for (int comp = 0; comp < a_lhs.nComp(); comp++){
	lhs(vofit(), comp) *= a_scale;
      }
    }
  }
}

void DataOps::divide(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_lcomp, const int a_rcomp){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::divide(*a_lhs[lvl], *a_rhs[lvl], a_lcomp, a_rcomp);
  }
}

void DataOps::divide(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int a_lcomp, const int a_rcomp){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs       = a_lhs[dit()];
    const EBCellFAB& rhs = a_rhs[dit()];

    lhs.divide(rhs, a_rcomp, a_lcomp, 1);
  }
}

void DataOps::divideByScalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::divideByScalar(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void DataOps::divideByScalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs){
  const int lcomps = a_lhs.nComp();
  const int rcomps = a_rhs.nComp();
  
  CH_assert(a_rhs.nComp() == 1);
  CH_assert(a_lhs.nComp() >= 1);

  for (int comp = 0; comp < lcomps; comp++){
    DataOps::divide(a_lhs, a_rhs, comp, 0);
  }
}

void DataOps::divideFallback(EBAMRCellData& a_numerator, const EBAMRCellData& a_denominator, const EBAMRCellData& a_fallback) {
  for (int lvl = 0; lvl < a_numerator.size(); lvl++){
    DataOps::divideFallback(*a_numerator[lvl], *a_denominator[lvl], *a_fallback[lvl]);
  }
}

void DataOps::divideFallback(LevelData<EBCellFAB>& a_numerator, const LevelData<EBCellFAB>& a_denominator, const LevelData<EBCellFAB>& a_fallback) {
  CH_assert(a_numerator.nComp() == a_denominator.nComp());
  CH_assert(a_numerator.nComp() == a_fallback.   nComp());
  
  for (DataIterator dit = a_numerator.dataIterator(); dit.ok(); ++dit){
    EBCellFAB&       numerator   = a_numerator  [dit()];
    const EBCellFAB& denominator = a_denominator[dit()];
    const EBCellFAB& fallback    = a_fallback   [dit()];

    BaseFab<Real>&       regNumerator   = numerator.  getSingleValuedFAB();
    const BaseFab<Real>& regDenominator = denominator.getSingleValuedFAB();
    const BaseFab<Real>& regFallback    = fallback.   getSingleValuedFAB();

    const Box cellBox = a_numerator.disjointBoxLayout()[dit()];

    // I need to clone the input data because the regular kernel will screw with it. 
    EBCellFAB cloneNumerator;
    cloneNumerator.clone(numerator);

    // Regular cells    
    for (int comp = 0; comp < a_numerator.nComp(); comp++){
      FORT_DIVIDE_FALLBACK(CHF_FRA1      (regNumerator,   comp),
			   CHF_CONST_FRA1(regDenominator, comp),
			   CHF_CONST_FRA1(regFallback,    comp),
			   CHF_BOX(cellBox));

    }
    
    // Irregular cells
    const EBISBox&   ebisBox  = numerator.getEBISBox();
    const EBGraph&   ebGraph  = ebisBox.getEBGraph();
    const IntVectSet irregIVS = ebisBox.getIrregIVS(cellBox);
    for (VoFIterator vofit(irregIVS, ebGraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      for (int comp = 0; comp < a_numerator.nComp(); comp++){
	const Real denom = denominator(vof, comp);
	if(std::abs(denom) > 0.0){
	  numerator(vof, comp) = cloneNumerator(vof, comp)/denom;
	}
	else{
	  numerator(vof, comp) = fallback(vof, comp);
	}
      }
    }
  }
}

void DataOps::floor(EBAMRCellData& a_lhs, const Real a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::floor(*a_lhs[lvl], a_value);
  }
}

void DataOps::floor(LevelData<EBCellFAB>& a_lhs, const Real a_value){
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

void DataOps::floor(EBAMRIVData& a_lhs, const Real a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::floor(*a_lhs[lvl], a_value);
  }
}

void DataOps::floor(LevelData<BaseIVFAB<Real> >& a_lhs, const Real a_value){
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

void DataOps::kappaSum(Real& a_mass, const LevelData<EBCellFAB>& a_lhs){


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

void DataOps::kappaScale(EBAMRCellData& a_data){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    DataOps::kappaScale(*a_data[lvl]);
  }
}

void DataOps::kappaScale(LevelData<EBCellFAB>& a_data){
  EBLevelDataOps::kappaWeight(a_data);
}

void DataOps::kappaScale(MFAMRCellData& a_data){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    DataOps::kappaScale(*a_data[lvl]);
  }
}

void DataOps::kappaScale(LevelData<MFCellFAB>& a_data){
  MFLevelDataOps::kappaWeight(a_data);
}

void DataOps::laplacian(EBAMRCellData& a_lapl, const EBAMRCellData& a_data){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    DataOps::laplacian(*a_lapl[lvl], *a_data[lvl]);
  }
}

void DataOps::laplacian(LevelData<EBCellFAB>& a_lapl, const LevelData<EBCellFAB>& a_data){
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

void DataOps::genLaplacian(EBAMRCellData& a_lapl, const EBAMRCellData& a_data){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    DataOps::genLaplacian(*a_lapl[lvl], *a_data[lvl]);
  }
}

void DataOps::genLaplacian(LevelData<EBCellFAB>& a_lapl, const LevelData<EBCellFAB>& a_data){
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

void DataOps::flashError(EBAMRCellData& a_lapl, const EBAMRCellData& a_data, const Real a_eps){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    DataOps::flashError(*a_lapl[lvl], *a_data[lvl], a_eps);
  }
}

void DataOps::flashError(LevelData<EBCellFAB>& a_lapl, const LevelData<EBCellFAB>& a_data, const Real a_eps){
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

void DataOps::ln(EBAMRCellData& a_lhs){
  for (int lvl = 0; lvl <= a_lhs.size(); lvl++){
    DataOps::ln(*a_lhs[lvl]);
  }
}

void DataOps::ln(LevelData<EBCellFAB>& a_lhs){
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

void DataOps::multiply(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::multiply(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void DataOps::multiply(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()] *= a_rhs[dit()];
  }
}

void DataOps::multiply(EBAMRFluxData& a_lhs, const EBAMRFluxData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::multiply(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void DataOps::multiply(LevelData<EBFluxFAB>& a_lhs, const LevelData<EBFluxFAB>& a_rhs){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()] *= a_rhs[dit()];
  }
}

void DataOps::multiply(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::multiply(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void DataOps::multiply(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs){
  const int ncomp = a_lhs.nComp();

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    BaseIVFAB<Real>& lhs       = a_lhs[dit()];
    const BaseIVFAB<Real>& rhs = a_rhs[dit()];
    const EBGraph& ebgraph     = lhs.getEBGraph();
    const IntVectSet& ivs      = lhs.getIVS() & rhs.getIVS();

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      for (int comp = 0; comp < ncomp; comp++){
	lhs(vof, comp) *= rhs(vof, comp);
      }
    }
  }
}

void DataOps::multiplyScalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::multiplyScalar(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void DataOps::multiplyScalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs){
  CH_assert(a_rhs.nComp() == 1);

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs       = a_lhs[dit()];
    const EBCellFAB& rhs = a_rhs[dit()];

    for (int comp = 0; comp < lhs.nComp(); comp++){
      lhs.mult(rhs, 0, comp, 1);
    }
  }
}

void DataOps::multiplyScalar(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::multiplyScalar(*a_lhs[lvl], *a_rhs[lvl]);
  }
}
  
void DataOps::multiplyScalar(LevelData<BaseIVFAB<Real> >& a_lhs, const LevelData<BaseIVFAB<Real> >& a_rhs){
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

void DataOps::norm(Real& a_norm, const LevelData<EBCellFAB>& a_data, const ProblemDomain& a_domain, const int a_p){
  Real volume;

  a_norm = EBLevelDataOps::kappaNorm(volume, a_data, EBLEVELDATAOPS_ALLVOFS, a_domain, a_p);
}

void DataOps::setCoveredValue(EBAMRCellData& a_lhs, const int a_comp, const Real a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::setCoveredValue(*a_lhs[lvl], a_comp, a_value);
  }
}

void DataOps::setCoveredValue(LevelData<EBCellFAB>& a_lhs, const int a_comp, const Real a_value){
  EBLevelDataOps::setCoveredVal(a_lhs, a_comp, a_value);
}

void DataOps::setValue(MFAMRCellData& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Vector<Real>& a_dx, const int a_comp){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::setValue(*a_lhs[lvl], a_function, a_probLo, a_dx[lvl], a_comp);
  }
}

void DataOps::setValue(LevelData<MFCellFAB>& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Real a_dx, const int a_comp){
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    
    MFCellFAB& lhs      = a_lhs[dit()];

    for (int i = 0; i < lhs.numPhases(); i++){
      EBCellFAB& phaseData        = lhs.getPhase(i);
      BaseFab<Real>& phaseDataFAB = phaseData.getSingleValuedFAB();

      const Box box           = dbl[dit()];
      const EBISBox& ebisbox  = phaseData.getEBISBox();
      const EBGraph& ebgraph  = ebisbox.getEBGraph();
      const IntVectSet& irreg = ebisbox.getIrregIVS(box);

      // Regular cells
      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv = bit();

	const RealVect pos = a_probLo + (0.5*RealVect::Unit + RealVect(iv))*a_dx;

	phaseDataFAB(iv, a_comp) = a_function(pos);
      }

      // Irregular cells
      for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const IntVect& iv   = vof.gridIndex();
      
	const RealVect pos = a_probLo + (0.5*RealVect::Unit + RealVect(iv))*a_dx + ebisbox.centroid(vof)*a_dx;

	phaseData(vof, a_comp) = a_function(pos);
      }
      
    }
  }
}

void DataOps::setValue(EBAMRCellData& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Vector<Real>& a_dx, const int a_comp){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::setValue(*a_lhs[lvl], a_function, a_probLo, a_dx[lvl], a_comp);
  }
}

void DataOps::setValue(LevelData<EBCellFAB>& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Real a_dx, const int a_comp){
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    EBCellFAB& lhs        = a_lhs[dit()];
    BaseFab<Real>& lhsFAB = lhs.getSingleValuedFAB();
    
    const Box box           = dbl[dit()];
    const EBISBox& ebisbox  = lhs.getEBISBox();
    const EBGraph& ebgraph  = ebisbox.getEBGraph();
    const IntVectSet& irreg = ebisbox.getIrregIVS(box);

    // Regular cells
    for (BoxIterator bit(box); bit.ok(); ++bit){
      const IntVect iv = bit();

      const RealVect pos = a_probLo + (0.5*RealVect::Unit + RealVect(iv))*a_dx;

      lhsFAB(iv, a_comp) = a_function(pos);
    }

    // Irregular cells
    for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const IntVect& iv   = vof.gridIndex();
      
      const RealVect pos = a_probLo + (0.5*RealVect::Unit + RealVect(iv))*a_dx + ebisbox.centroid(vof)*a_dx;

      lhs(vof, a_comp) = a_function(pos);
    }
  }
}

void DataOps::setValue(EBAMRFluxData& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Vector<Real>& a_dx, const int a_comp){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::setValue(*a_lhs[lvl], a_function, a_probLo, a_dx[lvl], a_comp);
  }
}

void DataOps::setValue(LevelData<EBFluxFAB>& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Real a_dx, const int a_comp){
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& lhs        = a_lhs[dit()][dir];
      BaseFab<Real>& lhsFAB = lhs.getSingleValuedFAB();
    
      const Box box           = dbl[dit()];
      const EBISBox& ebisbox  = lhs.getEBISBox();
      const EBGraph& ebgraph  = ebisbox.getEBGraph();
      const IntVectSet& irreg = ebisbox.getIrregIVS(box);

      Box facebox = box;
      facebox.surroundingNodes(dir);

      // Regular cells
      for (BoxIterator bit(facebox); bit.ok(); ++bit){
	const IntVect iv = bit();

	const RealVect pos = a_probLo + RealVect(iv)*a_dx;

	lhsFAB(iv, a_comp) = a_function(pos);
      }

      // Irregular cells. 
      const FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;
      for (FaceIterator faceIt(irreg, ebgraph, dir, stopCrit); faceIt.ok(); ++faceIt){
	const FaceIndex& face = faceIt();
	
	const RealVect pos = EBArith::getFaceLocation(face, a_dx, a_probLo);

	lhs(face, a_comp) = a_function(pos);
      }
    }
  }
}

void DataOps::setValue(EBAMRIVData& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Vector<Real>& a_dx, const int a_comp){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::setValue(*a_lhs[lvl], a_function, a_probLo, a_dx[lvl], a_comp);
  }
}

void DataOps::setValue(LevelData<BaseIVFAB<Real> >& a_lhs, const std::function<Real(const RealVect)>& a_function, const RealVect a_probLo, const Real a_dx, const int a_comp){
  // As we don't specify where the function should be evaluated, this routine sets a_lhs to be evaluated at the cell center.
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    BaseIVFAB<Real>& lhs    = a_lhs[dit()];

    const EBGraph& ebgraph  = lhs.getEBGraph();
    const IntVectSet& irreg = lhs.getIVS();

    for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
      const VolIndex vof = vofit();
      const IntVect  iv  = vof.gridIndex();
      const RealVect pos = a_probLo + (0.5*RealVect::Unit + RealVect(iv))*a_dx;

      lhs(vof, a_comp) = a_function(pos);
    }
  }
}

void DataOps::setValue(EBAMRCellData& a_lhs, const std::function<RealVect(const RealVect)>& a_function, const RealVect a_probLo, const Vector<Real>& a_dx){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::setValue(*a_lhs[lvl], a_function, a_probLo, a_dx[lvl]);
  }
}

void DataOps::setValue(LevelData<EBCellFAB>& a_lhs, const std::function<RealVect(const RealVect)>& a_function, const RealVect a_probLo, const Real a_dx){
  CH_assert(a_lhs.nComp() == SpaceDim);
  
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    EBCellFAB& lhs        = a_lhs[dit()];
    BaseFab<Real>& lhsFAB = lhs.getSingleValuedFAB();
    
    const Box box           = dbl[dit()];
    const EBISBox& ebisbox  = lhs.getEBISBox();
    const EBGraph& ebgraph  = ebisbox.getEBGraph();
    const IntVectSet& irreg = ebisbox.getIrregIVS(box);

    // Regular cells
    for (BoxIterator bit(box); bit.ok(); ++bit){
      const IntVect iv = bit();

      const RealVect pos = a_probLo + (0.5*RealVect::Unit + RealVect(iv))*a_dx;
      const RealVect val = a_function(pos);
      
      for (int comp = 0; comp < SpaceDim; comp++){
	lhsFAB(iv, comp) = val[comp];
      }
    }

    // Irregular cells
    for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const IntVect& iv   = vof.gridIndex();
      
      const RealVect pos = a_probLo + (0.5*RealVect::Unit + RealVect(iv))*a_dx + ebisbox.centroid(vof)*a_dx;
      const RealVect val = a_function(pos);

      for (int comp = 0; comp < SpaceDim; comp++){
	lhs(vof, comp) = val[comp];
      }
    }
  }
}

void DataOps::setValue(EBAMRCellData& a_data, const Real& a_value){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    EBLevelDataOps::setVal(*a_data[lvl], a_value);
  }
}

void DataOps::setValue(EBAMRCellData& a_lhs, const Real a_value, const int a_comp){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::setValue(*a_lhs[lvl], a_value, a_comp);
  }
}

void DataOps::setValue(LevelData<EBCellFAB>& a_lhs, const Real a_value, const int a_comp){
  EBLevelDataOps::setVal(a_lhs, a_value, a_comp);
}

void DataOps::setValue(LevelData<EBCellFAB>& a_lhs, const Real a_value){
  EBLevelDataOps::setVal(a_lhs, a_value);
}

void DataOps::setValue(LevelData<EBFluxFAB>& a_lhs, const Real a_value){
  EBLevelDataOps::setVal(a_lhs, a_value);
}
  
void DataOps::setValue(LevelData<BaseIVFAB<Real> >& a_lhs, const Real a_value){
  EBLevelDataOps::setVal(a_lhs, a_value);
}

void DataOps::setValue(EBAMRFluxData& a_data, const Real& a_value){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    EBLevelDataOps::setVal(*a_data[lvl], a_value);
  }
}

void DataOps::setValue(EBAMRIVData& a_data, const Real& a_value){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    EBLevelDataOps::setVal(*a_data[lvl], a_value);
  }
}

void DataOps::setValue(MFAMRCellData& a_lhs, const Real& a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::setValue(*a_lhs[lvl], a_value);
  }
}

void DataOps::setValue(LevelData<MFCellFAB>& a_lhs, const Real& a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    MFCellFAB& lhs = a_lhs[dit()];
    lhs.setVal(a_value);
  }
}

void DataOps::setValue(MFAMRFluxData& a_lhs, const Real& a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::setValue(*a_lhs[lvl] , a_value);
  }
}

void DataOps::setValue(LevelData<MFFluxFAB>& a_lhs, const Real& a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    MFFluxFAB& lhs = a_lhs[dit()];
    lhs.setVal(a_value);
  }
}

void DataOps::setValue(MFAMRIVData& a_lhs, const Real& a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::setValue(*a_lhs[lvl] , a_value);
  }
}

void DataOps::setValue(LevelData<MFBaseIVFAB>& a_lhs, const Real& a_value){
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()].setVal(a_value);
  }
}

void DataOps::setValue(EBAMRIFData& a_lhs, const Real a_value){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::setValue(*a_lhs[lvl], a_value);
  }
}

void DataOps::setValue(LevelData<DomainFluxIFFAB>& a_lhs, const Real a_value){
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

void DataOps::sum(Real& a_value){

#ifdef CH_MPI
  Real cur = a_value;
  Real tmp = a_value;
  int result = MPI_Allreduce(&cur, &tmp, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Abort("DataOps::sum - communication error on Allreduce");
  }

  a_value = tmp;
#endif
}

void DataOps::squareRoot(EBAMRFluxData& a_lhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::squareRoot(*a_lhs[lvl]);
  }
}

void DataOps::squareRoot(LevelData<EBFluxFAB>& a_lhs){
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

void DataOps::vectorLength(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::vectorLength(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void DataOps::vectorLength(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs){
  CH_assert(a_lhs.nComp() == 1);
  CH_assert(a_rhs.nComp() == SpaceDim);
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs             = a_lhs[dit()];
    const Box& box             = a_lhs.disjointBoxLayout().get(dit());
    const EBCellFAB& rhs       = a_rhs[dit()];

    DataOps::vectorLength(lhs, rhs, box);
  }
}

void DataOps::vectorLength(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& a_box){
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

void DataOps::vectorLength2(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::vectorLength2(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void DataOps::vectorLength2(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs){
  CH_assert(a_lhs.nComp() == 1);
  CH_assert(a_rhs.nComp() == SpaceDim);
  
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs             = a_lhs[dit()];
    const Box& box             = a_lhs.disjointBoxLayout().get(dit());
    const EBCellFAB& rhs       = a_rhs[dit()];

    DataOps::vectorLength2(lhs, rhs, box);
  }
}

void DataOps::vectorLength2(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& a_box){
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

void DataOps::computeMinValidBox(RealVect& a_lo, RealVect& a_hi, const RealVect a_normal, const RealVect a_centroid){
  const int num_segments = 10;

  // Default values
  a_lo = -0.5*RealVect::Unit;
  a_hi =  0.5*RealVect::Unit;

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const RealVect plane_normal = (sit() == Side::Lo) ? BASISREALV(dir) : -BASISREALV(dir);
      const RealVect plane_point  = RealVect::Zero - 0.5*plane_normal;      // Center point on plane
      const RealVect base_shift   = plane_normal/(1.0*num_segments);

#if CH_SPACEDIM == 2
      Vector<RealVect> corners(2);
      const int otherDir = (dir + 1) % SpaceDim;
      corners[0] = plane_point - 0.5*RealVect(BASISV(otherDir));
      corners[1] = plane_point + 0.5*RealVect(BASISV(otherDir));
#elif CH_SPACEDIM == 3
      Vector<RealVect> corners(4);
      const int otherDir1 = (dir + 1) % SpaceDim;
      const int otherDir2 = (dir + 2) % SpaceDim;
      corners[0] = plane_point - 0.5*RealVect(BASISV(otherDir1)) - 0.5*RealVect(BASISV(otherDir2));
      corners[1] = plane_point - 0.5*RealVect(BASISV(otherDir1)) + 0.5*RealVect(BASISV(otherDir2));
      corners[2] = plane_point + 0.5*RealVect(BASISV(otherDir1)) - 0.5*RealVect(BASISV(otherDir2));
      corners[3] = plane_point + 0.5*RealVect(BASISV(otherDir1)) + 0.5*RealVect(BASISV(otherDir2));
#endif

      // Shift corners in direction plane_normal with length base_shift. Keep track of the total
      // displacement of the plane. 
      RealVect shift_vector = RealVect::Zero;
      bool allInside = DataOps::allCornersInsideEb(corners, a_normal, a_centroid);

      while(allInside){

	// Shift the corners
	DataOps::shiftCorners(corners, base_shift);
	shift_vector += base_shift;

	// Check if shifted corners are inside EB
	allInside = DataOps::allCornersInsideEb(corners, a_normal, a_centroid);

	// If they are, we can change some components of a_lo
	if(allInside) {
	  if(sit() == Side::Lo){
	    a_lo[dir] = -0.5 + shift_vector[dir]; 
	  }
	  else if(sit() == Side::Hi){
	    a_hi[dir] = 0.5 + shift_vector[dir];
	  }
	}
      }
    }
  }
}

bool DataOps::allCornersInsideEb(const Vector<RealVect>& a_corners, const RealVect a_normal, const RealVect a_centroid){
  bool ret = true;

  // If any point it outside the EB, i.e. inside the domain boundary, return false. 
  for (int i = 0; i < a_corners.size(); i++){
    if(PolyGeom::dot((a_corners[i]-a_centroid), a_normal) > 0.0){
      ret = false;
    }
  }

  return ret;
}

void DataOps::shiftCorners(Vector<RealVect>& a_corners, const RealVect& a_distance){
  for(int i = 0; i < a_corners.size(); i++){
    a_corners[i] += a_distance;
  }
}

void DataOps::filterSmooth(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_stride, const Real a_alpha){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    DataOps::filterSmooth(*a_lhs[lvl], *a_rhs[lvl], a_stride, a_alpha);
  }
}

void DataOps::filterSmooth(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int a_stride, const Real a_alpha){

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  const IntVect ghostVec       = a_lhs.ghostVect();
  const int ncomp              = a_lhs.nComp();

  // Dummy check. 
  for (int dir = 0; dir < SpaceDim; dir++){
    if(ghostVec[dir] < a_stride) MayDay::Abort("DataOps::filterSmooth - stride is too large!");
  }

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& lhs       = a_lhs[dit()];
    const EBCellFAB& rhs = a_rhs[dit()];
    
    const Box box  = dbl[dit()];


    BaseFab<Real>&       lhs_fab = lhs.getSingleValuedFAB();
    const BaseFab<Real>& rhs_fab = rhs.getSingleValuedFAB();

    for (int icomp = 0; icomp < ncomp; icomp++){
      FORT_FILTER_SMOOTH(CHF_FRA1(lhs_fab, icomp),
			 CHF_CONST_FRA1(rhs_fab, icomp),
			 CHF_CONST_INT(a_stride),
			 CHF_CONST_REAL(a_alpha),
			 CHF_BOX(box));
    }

    // No filtering of irregular cells. 
    const EBISBox& ebisbox = lhs.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      for (int icomp = 0; icomp < ncomp; icomp++){
	lhs(vof, icomp) = rhs(vof, icomp);
      }
    }
  }
}

void DataOps::computeParticleWeights(unsigned long long&      a_weight,
				     unsigned long long&      a_num,
				     unsigned long long&      a_remainder,
				     const unsigned long long a_numPhysicalParticles,
				     const int                a_ppc) {

  if(a_numPhysicalParticles <= a_ppc){  
    a_weight    = 1;
    a_remainder = 0;
    a_num       = a_numPhysicalParticles; 
  }
  else{ // Add superparticles
    a_weight    = a_numPhysicalParticles/a_ppc;
    a_remainder = a_numPhysicalParticles%a_ppc;
    a_num       = (a_remainder == 0) ? a_ppc : a_ppc - 1;
  }
}

#include <CD_NamespaceFooter.H>
