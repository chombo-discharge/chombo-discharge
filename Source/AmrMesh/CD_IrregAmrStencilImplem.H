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
IrregAmrStencil<IrregSten>::IrregAmrStencil()
{
  CH_TIME("IrregAmrStencil::IrregAmrStencil");
  m_isDefined = false;
}

template <class IrregSten>
IrregAmrStencil<IrregSten>::IrregAmrStencil(const Vector<DisjointBoxLayout>& a_grids,
                                            const Vector<EBISLayout>&        a_ebisl,
                                            const Vector<ProblemDomain>&     a_domains,
                                            const Vector<Real>&              a_dx,
                                            const int                        a_finestLevel,
                                            const int                        a_order,
                                            const int                        a_radius,
                                            const IrregStencil::StencilType  a_type)
{
  CH_TIME("IrregAmrStencil::IrregAmrStencil");

  this->define(a_grids, a_ebisl, a_domains, a_dx, a_finestLevel, a_order, a_radius, a_type);
}

template <class IrregSten>
IrregAmrStencil<IrregSten>::~IrregAmrStencil()
{
  CH_TIME("IrregAmrStencil::~IrregAmrStencil");
}

template <class IrregSten>
const IrregStencil&
IrregAmrStencil<IrregSten>::operator[](const int& a_lvl) const
{
  return *m_stencils[a_lvl];
}

template <class IrregSten>
IrregStencil&
IrregAmrStencil<IrregSten>::operator[](const int& a_lvl)
{
  return *m_stencils[a_lvl];
}

template <class IrregSten>
void
IrregAmrStencil<IrregSten>::define(const Vector<DisjointBoxLayout>& a_grids,
                                   const Vector<EBISLayout>&        a_ebisl,
                                   const Vector<ProblemDomain>&     a_domains,
                                   const Vector<Real>&              a_dx,
                                   const int                        a_finestLevel,
                                   const int                        a_order,
                                   const int                        a_radius,
                                   const IrregStencil::StencilType  a_type)
{
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
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    m_stencils[lvl] = RefCountedPtr<IrregStencil>(
      new IrregSten(m_grids[lvl], m_ebisl[lvl], m_domains[lvl], m_dx[lvl], m_order, m_radius, m_stencilType));
  }

  m_isDefined = true;
}

template <class IrregSten>
void
IrregAmrStencil<IrregSten>::apply(LevelData<EBCellFAB>& a_dst, const LevelData<EBCellFAB>& a_src, const int a_lvl) const
{
  CH_TIME("IrregAmrStencil::apply(LD<EBCellFABx2, int)");

  CH_assert(m_isDefined);

  const DataIterator& dit = m_grids[a_lvl].dataIterator();

  const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    m_stencils[a_lvl]->apply(a_dst[din], a_src[din], din);
  }
}

template <class IrregSten>
void
IrregAmrStencil<IrregSten>::apply(LevelData<EBCellFAB>& a_data, const int a_lvl) const
{
  CH_TIME("IrregAmrStencil::apply(LD<EBCellFAB, int)");

  CH_assert(m_isDefined);

  const DataIterator& dit = m_grids[a_lvl].dataIterator();

  const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    EBCellFAB& data = a_data[din];
    EBCellFAB  cpy;

    // Need a clone to avoid race condition.
    cpy.clone(data);
    cpy.setVal(0.0);
    cpy += data;

    m_stencils[a_lvl]->apply(data, cpy, din);
  }
}

template <class IrregSten>
void
IrregAmrStencil<IrregSten>::apply(LevelData<BaseIVFAB<Real>>& a_dst,
                                  const LevelData<EBCellFAB>& a_src,
                                  const int                   a_lvl) const
{
  CH_TIME("IrregAmrStencil::apply(LD<BaseIVFAB>, LD<EBCellFAB>, int)");

  CH_assert(m_isDefined);

  const DataIterator& dit = m_grids[a_lvl].dataIterator();

  const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    m_stencils[a_lvl]->apply(a_dst[din], a_src[din], din);
  }
}

template <class IrregSten>
void
IrregAmrStencil<IrregSten>::apply(EBAMRCellData& a_data) const
{
  CH_TIME("IrregAmrStencil::apply");

  CH_assert(m_isDefined);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    this->apply(*a_data[lvl], lvl);
  }
}

template <class IrregSten>
void
IrregAmrStencil<IrregSten>::apply(EBAMRCellData& a_dst, const EBAMRCellData& a_src) const
{
  CH_TIME("IrregAmrStencil::apply");

  CH_assert(m_isDefined);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    this->apply(*a_dst[lvl], *a_src[lvl], lvl);
  }
}

template <class IrregSten>
void
IrregAmrStencil<IrregSten>::apply(EBAMRIVData& a_dst, const EBAMRCellData& a_src) const
{
  CH_TIME("IrregAmrStencil::apply");

  CH_assert(m_isDefined);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    this->apply(*a_dst[lvl], *a_src[lvl], lvl);
  }
}

template <class IrregSten>
void
IrregAmrStencil<IrregSten>::apply(Vector<EBAMRIVData*>& a_dst, const Vector<EBAMRCellData*>& a_src) const
{
  CH_TIME("IrregAmrStencil::apply");

  const int dstSize = a_dst.size();
  const int srcSize = a_src.size();

  CH_assert(m_isDefined);
  CH_assert(dstSize == srcSize);

  for (int i = 0; i < dstSize; i++) {
    this->apply(*a_dst[i], (const EBAMRCellData&)*a_src[i]);
  }
}

template <class IrregSten>
void
IrregAmrStencil<IrregSten>::apply(Vector<EBAMRCellData*>& a_data) const
{
  CH_TIME("IrregAmrStencil::apply");

  CH_assert(m_isDefined);

  const int dstSize = a_data.size();

  for (int i = 0; i < dstSize; i++) {
    this->apply(*a_data[i]);
  }
}

template <class IrregSten>
void
IrregAmrStencil<IrregSten>::apply(Vector<EBAMRCellData*>& a_dst, const Vector<EBAMRCellData*>& a_src) const
{
  CH_TIME("IrregAmrStencil::apply");

  const int dstSize = a_dst.size();
  const int srcSize = a_src.size();

  CH_assert(m_isDefined);
  CH_assert(dstSize == srcSize);

  for (int i = 0; i < dstSize; i++) {
    this->apply(*a_dst[i], (const EBAMRCellData&)*a_src[i]);
  }
}

#include <CD_NamespaceFooter.H>

#endif
