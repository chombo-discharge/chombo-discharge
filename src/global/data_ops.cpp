/*!
  @file data_ops.cpp
  @brief Implementation of data_ops.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "data_ops.H"
#include "EBLevelDataOps.H"
#include "MFLevelDataOps.H"

void data_ops::incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    EBLevelDataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void data_ops::copy(EBAMRCellData& a_dst, const EBAMRCellData& a_src){
  for (int lvl = 0; lvl < a_dst.size(); lvl++){
    if(a_src[lvl] != NULL && a_dst[lvl] != NULL){
      a_src[lvl]->copyTo(*a_dst[lvl]);
    }
  }
}

void data_ops::scale(MFAMRCellData& a_lhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    MFLevelDataOps::scale(*a_lhs[lvl], a_scale);
  }
}

void data_ops::scale(EBAMRIVData& a_lhs, const Real& a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    data_ops::scale(*a_lhs[lvl], a_scale);
  }
}

void data_ops::scale(EBAMRCellData& a_lhs, const Real a_scale){
  for (int lvl = 0; lvl < a_lhs.size(); lvl++){
    EBLevelDataOps::scale(*a_lhs[lvl], a_scale);
  }
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

void data_ops::set_value(EBAMRCellData& a_data, const Real& a_value){
  for (int lvl = 0; lvl < a_data.size(); lvl++){
    EBLevelDataOps::setVal(*a_data[lvl], a_value);
  }
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
