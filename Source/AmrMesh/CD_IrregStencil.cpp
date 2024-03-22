/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   IrregStencil.cpp
  @brief  Implementation of IrregStencil.H
  @author Robert Marskar
*/

// Chombo includes
#include <NeighborIterator.H>

// Our includes
#include <CD_IrregStencil.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int IrregStencil::m_defaultStenComp;
constexpr int IrregStencil::m_defaultNumSten;

IrregStencil::IrregStencil()
{}

IrregStencil::~IrregStencil()
{}

IrregStencil::IrregStencil(const DisjointBoxLayout&        a_dbl,
                           const EBISLayout&               a_ebisl,
                           const ProblemDomain&            a_domain,
                           const Real&                     a_dx,
                           const int                       a_order,
                           const int                       a_radius,
                           const IrregStencil::StencilType a_type)
{

  this->define(a_dbl, a_ebisl, a_domain, a_dx, a_order, a_radius, a_type);
}

const BaseIVFAB<VoFStencil>&
IrregStencil::operator[](const DataIndex& a_dit) const
{
  return *m_stencils[a_dit];
}

BaseIVFAB<VoFStencil>&
IrregStencil::operator[](const DataIndex& a_dit)
{
  return *m_stencils[a_dit];
}

void
IrregStencil::define(const DisjointBoxLayout&        a_dbl,
                     const EBISLayout&               a_ebisl,
                     const ProblemDomain&            a_domain,
                     const Real&                     a_dx,
                     const int                       a_order,
                     const int                       a_radius,
                     const IrregStencil::StencilType a_type)
{
  CH_TIME("IrregStencil::define");

  m_dbl         = a_dbl;
  m_ebisl       = a_ebisl;
  m_dx          = a_dx;
  m_radius      = a_radius;
  m_order       = a_order;
  m_stencilType = a_type;

  m_stencils.define(m_dbl);
  m_vofIter.define(m_dbl);

  const DataIterator& dit  = m_dbl.dataIterator();
  const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box&        box     = m_dbl[din];
    const EBISBox&    ebisbox = m_ebisl[din];
    const EBGraph&    ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs     = ebisbox.getIrregIVS(box);

    // Build the coarse-fine interface around this box
    IntVectSet       cfivs = IntVectSet(grow(box, 1) & m_domain);
    NeighborIterator nit(m_dbl);
    for (nit.begin(din); nit.ok(); ++nit) {
      cfivs -= m_dbl[nit()];
    }

    m_stencils[din] = RefCountedPtr<BaseIVFAB<VoFStencil>>(new BaseIVFAB<VoFStencil>(ivs, ebgraph, m_defaultNumSten));

    VoFIterator& vofit = m_vofIter[din];
    vofit.define(ivs, ebgraph);

    auto kernel = [&](const VolIndex& vof) -> void {
      VoFStencil& stencil = (*m_stencils[din])(vof, 0);
      this->buildStencil(stencil, vof, m_dbl, m_domain, ebisbox, box, m_dx, cfivs);

#if 0 // Safety test
      Real sum = 0.0;
      for (int i = 0; i < stencil.size(); i++){
	sum += stencil.weight(i);
      }

      if(Abs(sum - 1.0) > 1.E-6){
	MayDay::Warning("IrregStencil::define - weights do not sum to 1. Something has gone wrong with one of your stencils");
	stencil *= 1./sum;
      }
#endif
    };

    BoxLoops::loop(vofit, kernel);
  }
}

void
IrregStencil::apply(EBCellFAB& a_dst, const EBCellFAB& a_src, const DataIndex& a_dit) const
{
  CH_TIME("IrregStencil::apply");

  const BaseIVFAB<VoFStencil>& stencils = *m_stencils[a_dit];

  VoFIterator& vofit = m_vofIter[a_dit];

  auto kernel = [&](const VolIndex& vof) -> void {
    const VoFStencil& stencil = stencils(vof, m_defaultStenComp);

    for (int comp = 0; comp < a_dst.nComp(); comp++) {

      a_dst(vof, comp) = 0.0;

      for (int i = 0; i < stencil.size(); i++) {
        const VolIndex& ivof    = stencil.vof(i);
        const Real&     iweight = stencil.weight(i);

        a_dst(vof, comp) += iweight * a_src(ivof, comp);
      }
    }
  };

  BoxLoops::loop(vofit, kernel);
}

void
IrregStencil::apply(BaseIVFAB<Real>& a_dst, const EBCellFAB& a_src, const DataIndex& a_dit) const
{
  CH_TIME("IrregStencil::apply");

  const BaseIVFAB<VoFStencil>& stencils = *m_stencils[a_dit];

  VoFIterator& vofit = m_vofIter[a_dit];

  auto kernel = [&](const VolIndex& vof) -> void {
    const VoFStencil& stencil = stencils(vof, m_defaultStenComp);

    for (int comp = 0; comp < a_dst.nComp(); comp++) {

      a_dst(vof, comp) = 0.0;

      for (int i = 0; i < stencil.size(); i++) {
        const VolIndex& ivof    = stencil.vof(i);
        const Real&     iweight = stencil.weight(i);

        a_dst(vof, comp) += iweight * a_src(ivof, comp);
      }
    }
  };

  BoxLoops::loop(vofit, kernel);
}

#include <CD_NamespaceFooter.H>
