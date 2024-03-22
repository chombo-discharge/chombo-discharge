/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAddOp.cpp
  @brief  Implementation of CD_EBAddOp.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBAddOp.H>
#include <CD_NamespaceHeader.H>

EBAddOp::EBAddOp()
{
  CH_TIME("EBAddOp::EBAddOp");
}

EBAddOp::~EBAddOp()
{
  CH_TIME("EBAddOp::~EBAddOp");
}

void
EBAddOp::linearIn(EBCellFAB& a_data, void* a_buffer, const Box& a_region, const Interval& a_comps) const
{
  CH_TIME("EBAddOp::linearIn");

  // Linearize the input buffer onto an EBCellFAB.
  EBCellFAB incr;
  incr.clone(a_data);
  incr.setVal(0.);
  incr.linearIn(a_buffer, a_region, a_comps);

  // Add the input buffer to a_data.
  const int isrc = a_comps.begin();
  const int idst = a_comps.begin();
  const int inco = a_comps.size();

  a_data.plus(incr, isrc, idst, inco);
}

void
EBAddOp::op(EBCellFAB&       a_dst,
            const Box&       a_regionFrom,
            const Interval&  a_dstVars,
            const Box&       a_regionTo,
            const EBCellFAB& a_src,
            const Interval&  a_srcVars) const
{
  CH_TIME("EBAddOp::op");

  CH_assert(a_dstVars.size() == a_srcVars.size());

  // I think restricting the addition to a_regionTo is correct -- we are only trying to add into valid region
  // and excluding adding into ghost cells.

  const int isrc = a_srcVars.begin();
  const int idst = a_dstVars.begin();
  const int inco = a_dstVars.size();

  //  a_dst.plus(a_src, isrc, idst, inco);
  a_dst.plus(a_src, a_regionTo, isrc, idst, inco);
}

#include <CD_NamespaceFooter.H>
