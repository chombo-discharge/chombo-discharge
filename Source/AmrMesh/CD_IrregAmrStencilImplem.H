/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_IrregAmrStencilImplem.H
  @brief  Implementation of CD_IrregAmrStencil.H
  @author Robert Marskar
*/

#ifndef CD_IrregAmrStencilImplem_H
#define CD_IrregAmrStencilImplem_H

// Our includes
#include <CD_IrregAmrStencil.H>
#include <CD_NamespaceHeader.H>

template <class IrregSten>
constexpr int IrregAmrStencil<IrregSten>::stenComp;

template <class IrregSten>
IrregAmrStencil<IrregSten>::IrregAmrStencil(){
  CH_TIME("IrregAmrStencil::IrregAmrStencil");
  m_isDefined = false;
}

template<class IrregSten>
IrregAmrStencil<IrregSten>::IrregAmrStencil(const Vector<DisjointBoxLayout>& a_grids,
					    const Vector<EBISLayout>&        a_ebisl,
					    const Vector<ProblemDomain>&     a_domains,
					    const Vector<Real>&              a_dx,
					    const int                        a_finestLevel,
					    const int                        a_order,
					    const int                        a_radius,
					    const IrregStencil::StencilType  a_type){
  CH_TIME("IrregAmrStencil::IrregAmrStencil");
  
  this->define(a_grids, a_ebisl, a_domains, a_dx, a_finestLevel, a_order, a_radius, a_type);
}

template <class IrregSten>
IrregAmrStencil<IrregSten>::~IrregAmrStencil(){
  CH_TIME("IrregAmrStencil::~IrregAmrStencil");
}

template <class IrregSten>
const IrregStencil& IrregAmrStencil<IrregSten>::operator[](const int& a_lvl) const {
  return *m_stencils[a_lvl];
}

template <class IrregSten>
IrregStencil& IrregAmrStencil<IrregSten>::operator[](const int& a_lvl) {
  return *m_stencils[a_lvl];
}

template<class IrregSten>
void IrregAmrStencil<IrregSten>::define(const Vector<DisjointBoxLayout>& a_grids,
					const Vector<EBISLayout>&        a_ebisl,
					const Vector<ProblemDomain>&     a_domains,
					const Vector<Real>&              a_dx,
					const int                        a_finestLevel,
					const int                        a_order,
					const int                        a_radius,
					const IrregStencil::StencilType  a_type){
  CH_TIME("IrregAmrStencil::define");

  // Store input arguments. 
  m_grids       = a_grids;
  m_ebisl       = a_ebisl;  
  m_domains     = a_domains;
  m_dx          = a_dx;
  m_finestLevel = a_finestLevel;
  m_order       = a_order;
  m_radius      = a_radius;  
  m_stencilType = a_type;  

  // Define stencils on each level
  m_stencils.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    m_stencils[lvl] = RefCountedPtr<IrregStencil>(new IrregSten(m_grids  [lvl],
								m_ebisl  [lvl],
								m_domains[lvl],
								m_dx     [lvl],
								m_order,
								m_radius,
								m_stencilType));
  }
  
  m_isDefined = true;
}

template<class IrregSten>
void IrregAmrStencil<IrregSten>::apply(      LevelData<EBCellFAB>& a_dst,
				       const LevelData<EBCellFAB>& a_src,
				       const int                   a_lvl) const {
  CH_TIME("IrregAmrStencil::apply");

  CH_assert(m_isDefined);

  const DisjointBoxLayout& dbl     =  m_grids   [a_lvl];
  const EBISLayout&        ebisl   =  m_ebisl   [a_lvl];
  const IrregStencil&      irrSten = *m_stencils[a_lvl];
  const int                nComp   =  a_dst.nComp();
  
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const BaseIVFAB<VoFStencil>& stencils = irrSten[dit()];    

    EBCellFAB& dst       = a_dst[dit()];
    const EBCellFAB& src = a_src[dit()];

    // Apply stencil in each cut-cell. 
    VoFIterator& vofit = irrSten.getVofIterator(dit());
    for (vofit.reset(); vofit.ok(); ++vofit){    
      const VolIndex& vof       = vofit();
      const VoFStencil& vofSten = stencils(vof, stenComp);

      for (int comp = 0; comp < nComp; comp++){
	dst(vof, comp) = 0.;

	for (int i = 0; i < vofSten.size(); i++){
	  const VolIndex& ivof    = vofSten.vof(i);
	  const Real&     iweight = vofSten.weight(i);
	  
	  dst(vof,comp) += iweight * src(ivof, comp);
	}
      }
    }
  }
}

template<class IrregSten>
void IrregAmrStencil<IrregSten>::apply(LevelData<EBCellFAB>& a_data,
				       const int             a_lvl) const {
  CH_TIME("IrregAmrStencil::apply");

  CH_assert(m_isDefined);  

  const DisjointBoxLayout& dbl     =  m_grids   [a_lvl];
  const EBISLayout&        ebisl   =  m_ebisl   [a_lvl];
  const IrregStencil&      irrSten = *m_stencils[a_lvl];
  const int                nComp   =  a_data.nComp();
  
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const BaseIVFAB<VoFStencil>& stencils = irrSten[dit()];

    EBCellFAB& data = a_data[dit()];
    EBCellFAB  cpy;

    // Need a clone to avoid race condition. 
    cpy.clone(data);
    cpy.setVal(0.0);
    cpy += data;    

    // Apply stencil in each cut-cell. 
    VoFIterator& vofit = irrSten.getVofIterator(dit());
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof       = vofit();
      const VoFStencil& vofSten = stencils(vof, stenComp);

      for (int comp = 0; comp < nComp; comp++){
	data(vof, comp) = 0.;

	for (int i = 0; i < vofSten.size(); i++){
	  const VolIndex& ivof    = vofSten.vof(i);
	  const Real&     iweight = vofSten.weight(i);
	  
	  data(vof,comp) += iweight * cpy(ivof, comp);
	}
      }
    }
  }
}

template<class IrregSten>
void IrregAmrStencil<IrregSten>::apply(LevelData<BaseIVFAB<Real> >& a_dst,
				       const LevelData<EBCellFAB>&  a_src,
				       const int                    a_lvl) const {
  CH_TIME("IrregAmrStencil::apply");

  CH_assert(m_isDefined);

  const DisjointBoxLayout& dbl     =  m_grids   [a_lvl];
  const EBISLayout&        ebisl   =  m_ebisl   [a_lvl];
  const IrregStencil&      irrSten = *m_stencils[a_lvl];
  const int                nComp   =  a_dst.nComp();  

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const BaseIVFAB<VoFStencil>& stencils = irrSten[dit()];
    
    const EBCellFAB& src = a_src[dit()];
    BaseIVFAB<Real>& dst = a_dst[dit()];

    // Apply stencil in each cut-cell. 
    VoFIterator& vofit = irrSten.getVofIterator(dit());
    for (vofit.reset(); vofit.ok(); ++vofit){          
      const VolIndex& vof       = vofit();
      const VoFStencil& vofSten = stencils(vof, stenComp);

      for (int comp = 0; comp < nComp; comp++){
	dst(vof, comp) = 0.0;

	for (int i = 0; i < vofSten.size(); i++){
	  const VolIndex& ivof    = vofSten.vof(i);
	  const Real&     iweight = vofSten.weight(i);
	  
	  dst(vof,comp) += iweight * src(ivof, comp);
	}
      }
    }
  }
}

template<class IrregSten>
void IrregAmrStencil<IrregSten>::apply(EBAMRCellData& a_data) const {
  CH_TIME("IrregAmrStencil::apply");

  CH_assert(m_isDefined);    

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    this->apply(*a_data[lvl], lvl);
  }
}

template<class IrregSten>
void IrregAmrStencil<IrregSten>::apply(      EBAMRCellData& a_dst,
				       const EBAMRCellData& a_src) const {
  CH_TIME("IrregAmrStencil::apply");

  CH_assert(m_isDefined);    

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    this->apply(*a_dst[lvl], *a_src[lvl], lvl);
  }
}

template<class IrregSten>
void IrregAmrStencil<IrregSten>::apply(      EBAMRIVData&   a_dst,
				       const EBAMRCellData& a_src) const {
  CH_TIME("IrregAmrStencil::apply");

  CH_assert(m_isDefined);    

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    this->apply(*a_dst[lvl], *a_src[lvl], lvl);
  }
}

template<class IrregSten>
void IrregAmrStencil<IrregSten>::apply(      Vector<EBAMRIVData*>&   a_dst,
				       const Vector<EBAMRCellData*>& a_src) const {
  CH_TIME("IrregAmrStencil::apply");
  
  const int dstSize = a_dst.size();
  const int srcSize = a_src.size();

  CH_assert(m_isDefined);      
  CH_assert(dstSize == srcSize);

  for (int i = 0; i < dstSize; i++){
    EBAMRIVData& dst = *a_dst[i];
    const EBAMRCellData& src = *(const_cast<const EBAMRCellData*> (a_src[i]));

    this->apply(dst, src);
  }
}

template<class IrregSten>
void IrregAmrStencil<IrregSten>::apply(Vector<EBAMRCellData*>& a_data) const {
  CH_TIME("IrregAmrStencil::apply");

  CH_assert(m_isDefined);        

  const int dstSize = a_data.size();

  for (int i = 0; i < dstSize; i++){
    this->apply(*a_data[i]);
  }
}

template<class IrregSten>
void IrregAmrStencil<IrregSten>::apply(      Vector<EBAMRCellData*>& a_dst,
				       const Vector<EBAMRCellData*>& a_src) const {
  CH_TIME("IrregAmrStencil::apply");

  const int dstSize = a_dst.size();
  const int srcSize = a_src.size();

  CH_assert(m_isDefined);          
  CH_assert(dstSize == srcSize);

  for (int i = 0; i < dstSize; i++){
    EBAMRCellData& dst = *a_dst[i];
    const EBAMRCellData& src = *(const_cast<const EBAMRCellData*> (a_src[i]));

    this->apply(dst, src);
  }
}

#include <CD_NamespaceFooter.H>

#endif
