/*!
  @file irreg_amr_stencilI.hpp
  @brief Implementation of irreg_amr_stencil
  @author Robert Marskar
  @date Aug. 207
*/

#include "irreg_amr_stencil.H"
#include "centroid_interp.H"
#include "eb_centroid_interp.H"


template <class T> irreg_amr_stencil<T>::irreg_amr_stencil(){
  CH_TIME("irreg_amr_stencil::irreg_amr_stencil");
  m_defined = false;
}

//
template<class T> irreg_amr_stencil<T>::irreg_amr_stencil(const Vector<DisjointBoxLayout>& a_grids,
							  const Vector<EBISLayout>&        a_ebisl,
							  const Vector<ProblemDomain>&     a_domains,
							  const Vector<Real>&              a_dx,
							  const int                        a_finestLevel,
							  const int                        a_order,
							  const int                        a_radius){
  CH_TIME("irreg_amr_stencil::irreg_amr_stencil");
  
  this->define(a_grids, a_ebisl, a_domains, a_dx, a_finestLevel, a_order, a_radius);
}

template <class T> irreg_amr_stencil<T>::~irreg_amr_stencil(){
  CH_TIME("irreg_amr_stencil::~irreg_amr_stencil");
}

//
template<class T> void irreg_amr_stencil<T>::define(const Vector<DisjointBoxLayout>& a_grids,
						    const Vector<EBISLayout>&        a_ebisl,
						    const Vector<ProblemDomain>&     a_domains,
						    const Vector<Real>&              a_dx,
						    const int                        a_finestLevel,
						    const int                        a_order,
						    const int                        a_radius){
  CH_TIME("irreg_amr_stencil::define");

  m_finest_level = a_finestLevel;
  m_radius      = a_radius;
  m_order       = a_order;

  m_grids.resize(m_finest_level + 1);
  m_ebisl.resize(m_finest_level + 1);
  m_domains.resize(m_finest_level + 1);
  m_dx.resize(m_finest_level + 1);
  m_stencils.resize(m_finest_level + 1);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_grids[lvl]   = a_grids[lvl];
    m_ebisl[lvl]   = a_ebisl[lvl];
    m_domains[lvl] = a_domains[lvl];
    m_dx[lvl]      = a_dx[lvl];
    m_stencils[lvl] = RefCountedPtr<irreg_stencil>(new T(m_grids[lvl],
							 m_ebisl[lvl],
							 m_domains[lvl],
							 m_dx[lvl],
							 m_order,
							 m_radius));
  }
  
  m_defined = true;
}

//
template<class T> void irreg_amr_stencil<T>::apply(LevelData<EBCellFAB>&       a_dst,
							    const LevelData<EBCellFAB>& a_src,
							    const int                   a_lvl,
							    const bool                  a_limit) const {
  CH_TIME("irreg_amr_stencil::apply");

  const DisjointBoxLayout& dbl  = m_grids[a_lvl];
  const EBISLayout& ebisl       = m_ebisl[a_lvl];
  const irreg_stencil& irr_sten = *m_stencils[a_lvl];
  const int nComp               = a_dst.nComp();

  DataIterator dit = dbl.dataIterator();
  const int nbox = dit.size();
  
  for (int ibox = 0; ibox < nbox; ibox++){
    const BaseIVFAB<VoFStencil>& stencils = irr_sten[dit[ibox]];
    
    const Box& box          = dbl.get(dit[ibox]);
    const EBISBox& ebisbox  = ebisl[dit[ibox]];
    const EBGraph& ebgraph  = ebisbox.getEBGraph();
    const IntVectSet& ivs   = stencils.getIVS();

    EBCellFAB& dst_fab       = a_dst[dit[ibox]];
    const EBCellFAB& src_fab = a_src[dit[ibox]];

    dst_fab.copy(src_fab);

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof       = vofit();
      const VoFStencil& vofSten = stencils(vof, 0);

      //
      for (int comp = 0; comp < nComp; comp++){
    	dst_fab(vof, comp) = 0.;

    	for (int i = 0; i < vofSten.size(); i++){
    	  const VolIndex& ivof = vofSten.vof(i);
    	  const Real& iweight  = vofSten.weight(i);
    	  dst_fab(vof,comp) += src_fab(ivof, comp)*iweight;
    	}


	if(a_limit){ // Limit 
	  const bool signChange = dst_fab(vof, comp)*src_fab(vof, comp) < 0.;
	  dst_fab(vof, comp) = signChange ? 0. : dst_fab(vof, comp);
	}
      }
    }
  }
}

//
template<class T> void irreg_amr_stencil<T>::apply(LevelData<EBCellFAB>& a_data,
							    const int             a_lvl,
							    const bool            a_limit) const {
  CH_TIME("irreg_amr_stencil::apply");

  const DisjointBoxLayout& dbl = m_grids[a_lvl];
  const EBISLayout& ebisl      = m_ebisl[a_lvl];
  const irreg_stencil& irrSten  = *m_stencils[a_lvl];
  const int nComp              = a_data.nComp();

  DataIterator dit = dbl.dataIterator();
  const int nbox = dit.size();
  
  for (int ibox = 0; ibox < nbox; ibox++){
    const Box& box         = dbl.get(dit[ibox]);
    const EBISBox& ebisbox = ebisl[dit[ibox]];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const BaseIVFAB<VoFStencil>& stencils = irrSten[dit[ibox]];
    const IntVectSet& ivs  = stencils.getIVS();

    EBCellFAB& data = a_data[dit[ibox]];
    EBCellFAB  cpy;

    cpy.clone(data); // Clone data
    cpy.setVal(0.0); // Copy
    cpy += data;     // data

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof       = vofit();
      const VoFStencil& vofSten = stencils(vof, 0);

      //
      for (int comp = 0; comp < nComp; comp++){
    	data(vof, comp) = 0.;

    	for (int i = 0; i < vofSten.size(); i++){
    	  const VolIndex& ivof = vofSten.vof(i);
    	  const Real& iweight  = vofSten.weight(i);
    	  data(vof,comp) += cpy(ivof, comp)*iweight;
    	}


	if(a_limit){ // Limit extrapolation
	  const bool signChange = data(vof, comp)*cpy(vof, comp) < 0.;
	  data(vof, comp) = signChange ? 0. : data(vof, comp);
	}
      }
    }
  }
}

//
template<class T> void irreg_amr_stencil<T>::apply(LevelData<BaseIVFAB<Real> >& a_dst,
							    const LevelData<EBCellFAB>&  a_src,
							    const int                    a_lvl,
							    const bool                   a_limit) const {
  CH_TIME("irreg_amr_stencil::apply");

  //
  const DisjointBoxLayout& dbl = m_grids[a_lvl];
  const EBISLayout& ebisl      = m_ebisl[a_lvl];
  const irreg_stencil& irrSten  = *m_stencils[a_lvl];
  const int nComp              = a_dst.nComp();

  DataIterator dit = dbl.dataIterator();
  const int nbox = dit.size();
  
  for (int ibox = 0; ibox < nbox; ibox++){
    const Box& box          = dbl.get(dit[ibox]);
    const EBISBox& ebisbox  = ebisl[dit[ibox]];
    const EBGraph& ebgraph  = ebisbox.getEBGraph();
    const BaseIVFAB<VoFStencil>& stencils = irrSten[dit[ibox]];
    const IntVectSet& ivs   = stencils.getIVS();
    const EBCellFAB& src    = a_src[dit[ibox]];

    BaseIVFAB<Real>& dst    = a_dst[dit[ibox]];

      
    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof       = vofit();
      const VoFStencil& vofSten = stencils(vof, 0);

      //
      for (int comp = 0; comp < nComp; comp++){
	dst(vof, comp) = 0.0;

	for (int i = 0; i < vofSten.size(); i++){
	  const VolIndex& ivof = vofSten.vof(i);
	  const Real& iweight  = vofSten.weight(i);
	  dst(vof,comp) += src(ivof, comp)*iweight;
	}

	if(a_limit){ // Limit extrapolation
	  const bool signChange = dst(vof, comp)*src(vof, comp) < 0.;
	  dst(vof, comp) = signChange ? 0. : dst(vof, comp);
	}
      }
    }
  }
}

template<class T> void irreg_amr_stencil<T>::apply(EBAMRCellData& a_data,
							    const bool   a_limit) const {
  CH_TIME("irreg_amr_stencil::apply");

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    irreg_amr_stencil::apply(*a_data[lvl], lvl, a_limit);
  }
}

//
template<class T> void irreg_amr_stencil<T>::apply(EBAMRCellData&       a_dst,
							    const EBAMRCellData& a_src,
							    const bool         a_limit) const {
  CH_TIME("irreg_amr_stencil::apply");

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    irreg_amr_stencil::apply(*a_dst[lvl], *a_src[lvl], lvl, a_limit);
  }
}

//
template<class T> void irreg_amr_stencil<T>::apply(EBAMRIVData&         a_dst,
							    const EBAMRCellData& a_src,
							    const bool         a_limit) const {
  CH_TIME("irreg_amr_stencil::apply");

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    irreg_amr_stencil::apply(*a_dst[lvl], *a_src[lvl], lvl, a_limit);
  }
}

//
template<class T> void irreg_amr_stencil<T>::apply(Vector<EBAMRIVData*>&         a_dst,
							    const Vector<EBAMRCellData*>& a_src,
							    const bool                  a_limit) const {
  CH_TIME("irreg_amr_stencil::apply");

  const int dstSize = a_dst.size();
  const int srcSize = a_src.size();
  CH_assert(dstSize == srcSize);

  for (int i = 0; i < dstSize; i++){
    EBAMRIVData& dst = *a_dst[i];
    const EBAMRCellData& src = *(const_cast<const EBAMRCellData*> (a_src[i]));

    irreg_amr_stencil::apply(dst, src, a_limit);
  }
}

template<class T> void irreg_amr_stencil<T>::apply(Vector<EBAMRCellData*>& a_data,
							    const bool            a_limit) const {
  CH_TIME("irreg_amr_stencil::apply");

  const int dstSize = a_data.size();

  for (int i = 0; i < dstSize; i++){
    irreg_amr_stencil::apply(*a_data[i], a_limit);
  }
}

template<class T> void irreg_amr_stencil<T>::apply(Vector<EBAMRCellData*>&       a_dst,
							    const Vector<EBAMRCellData*>& a_src,
							    const bool                  a_limit) const {
  CH_TIME("irreg_amr_stencil::apply");

  const int dstSize = a_dst.size();
  const int srcSize = a_src.size();
  CH_assert(dstSize == srcSize);

  for (int i = 0; i < dstSize; i++){
    EBAMRCellData& dst = *a_dst[i];
    const EBAMRCellData& src = *(const_cast<const EBAMRCellData*> (a_src[i]));

    apply(dst, src, a_limit);
  }
}

template class irreg_amr_stencil<centroid_interp>;
template class irreg_amr_stencil<eb_centroid_interp>;
// template class irreg_amr_stencil<IrregCenterInterpolationStencil>;
