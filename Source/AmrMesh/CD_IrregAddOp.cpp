/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_IrregAddOp.cpp
  @brief  Implementation of CD_IrregAddOp.H
  @author Robert Marskar
*/

// Chombo includes
#include <VoFIterator.H>
#include <CH_Timer.H>

// Our includes
#include <CD_IrregAddOp.H>
#include <CD_NamespaceHeader.H>

IrregAddOp::IrregAddOp() noexcept { CH_TIME("IrregAddOp::IrregAddOp"); }

IrregAddOp::~IrregAddOp() noexcept { CH_TIME("IrregAddOp::~IrregAddOp"); }

void
IrregAddOp::linearIn(BaseIVFAB<Real>& a_data, void* a_buffer, const Box& a_region, const Interval& a_comps) const
{
  CH_TIME("IrregAddOp::linearIn");

  // Clone the data holder and linearize the buffer onto it.
  BaseIVFAB<Real> clone(a_data.getIVS(), a_data.getEBGraph(), a_comps.size());
  clone.linearIn(a_buffer, a_region, a_comps);

  // Increment a_data with the buffer.
  VoFIterator vofit(a_data.getIVS() & a_region, a_data.getEBGraph());

  const int srcBegin = a_comps.begin();
  const int dstBegin = a_comps.begin();
  const int numComp  = a_comps.size();

  for (int icomp = 0; icomp < numComp; icomp++) {
    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex& vof = vofit();

      a_data(vof, dstBegin + icomp) += clone(vof, srcBegin + icomp);
    }
  }
}

void
IrregAddOp::op(BaseIVFAB<Real>&       a_dst,
               const Box&             a_regionFrom,
               const Interval&        a_dstVars,
               const Box&             a_regionTo,
               const BaseIVFAB<Real>& a_src,
               const Interval&        a_srcVars) const
{
  CH_TIME("IrregAddOp::op");

  CH_assert(a_dstVars.size() == a_srcVars.size());

  const IntVectSet ivs     = a_dst.getIVS() & a_regionTo;
  const EBGraph&   ebgraph = a_dst.getEBGraph();

  const int srcBegin = a_srcVars.begin();
  const int dstBegin = a_dstVars.begin();
  const int numComp  = a_dstVars.size();

  VoFIterator vofit(ivs, ebgraph);

  for (int icomp = 0; icomp < numComp; icomp++) {
    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex& vof = vofit();

      a_dst(vof, dstBegin + icomp) += a_src(vof, srcBegin + icomp);
    }
  }
}

#include <CD_NamespaceFooter.H>
