/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CentroidInterpolationStencil.cpp
  @brief  Implementation of CD_CentroidInterpolationStencil.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_CentroidInterpolationStencil.H>
#include <CD_LinearStencil.H>
#include <CD_LeastSquares.H>
#include <CD_NamespaceHeader.H>

#define DEBUG_CENTROID_INTERP 0

CentroidInterpolationStencil::CentroidInterpolationStencil(const DisjointBoxLayout&        a_dbl,
                                                           const EBISLayout&               a_ebisl,
                                                           const ProblemDomain&            a_domain,
                                                           const Real&                     a_dx,
                                                           const int                       a_order,
                                                           const int                       a_radius,
                                                           const IrregStencil::StencilType a_type)
  : IrregStencil()
{

  CH_TIME("CentroidInterpolationStencil::CentroidInterpolationStencil");

  this->define(a_dbl, a_ebisl, a_domain, a_dx, a_order, a_radius, a_type);
}

CentroidInterpolationStencil::~CentroidInterpolationStencil()
{
  CH_TIME("CentroidInterpolationStencil::~CentroidInterpolationStencil");
}

void
CentroidInterpolationStencil::buildStencil(VoFStencil&              a_sten,
                                           const VolIndex&          a_vof,
                                           const DisjointBoxLayout& a_dbl,
                                           const ProblemDomain&     a_domain,
                                           const EBISBox&           a_ebisbox,
                                           const Box&               a_box,
                                           const Real&              a_dx,
                                           const IntVectSet&        a_cfivs)
{
  CH_TIME("CentroidInterpolationStencil::buildStencil");

  bool found_stencil = false;

  IntVectSet noCFIVS = IntVectSet();

  switch (m_stencilType) {
  case IrregStencil::StencilType::Linear: {
    found_stencil = LinearStencil::getLinearInterpStencil(a_sten,
                                                          a_ebisbox.centroid(a_vof),
                                                          a_vof,
                                                          a_domain,
                                                          a_ebisbox);
    break;
  }
  case IrregStencil::StencilType::TaylorExtrapolation: {
    found_stencil = this
                      ->getTaylorExtrapolationStencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, noCFIVS);
    break;
  }
  case IrregStencil::StencilType::LeastSquares: {
    found_stencil = this->getLeastSquaresStencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, noCFIVS);
    break;
  }
  case IrregStencil::StencilType::PiecewiseLinear: {
    found_stencil = this->getPiecewiseLinearStencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, noCFIVS);
    break;
  }
  default: {
    MayDay::Abort("CentroidInterpolationStencil::buildStencil - Unsupported stencil type");
    break;
  }
  }

  // Drop to zeroth order if we couldn't find a stencil.
  if (!found_stencil) {
    a_sten.clear();
    a_sten.add(a_vof, 1.0);
  }

#if DEBUG_CENTROID_INTERP
  for (int i = 0; i < a_sten.size(); i++) {
    const Real w = a_sten.weight(i);
    if (w < 0.0)
      MayDay::Warning("CentroidInterpolationStencil::buildStencils - I got negative weights!!!");
  }
#endif
}

bool
CentroidInterpolationStencil::getTaylorExtrapolationStencil(VoFStencil&              a_sten,
                                                            const VolIndex&          a_vof,
                                                            const DisjointBoxLayout& a_dbl,
                                                            const ProblemDomain&     a_domain,
                                                            const EBISBox&           a_ebisbox,
                                                            const Box&               a_box,
                                                            const Real&              a_dx,
                                                            const IntVectSet&        a_cfivs)
{
  CH_TIME("CentroidInterpolationStencil::getTaylorExtrapolationStencil");

  const int       comp     = 0;
  const RealVect& centroid = a_ebisbox.centroid(a_vof);
  IntVectSet*     cfivs    = const_cast<IntVectSet*>(&a_cfivs);

  int order;
  if (m_order == 1) {
    order = EBArith::getFirstOrderExtrapolationStencil(a_sten,
                                                       centroid * a_dx,
                                                       a_dx * RealVect::Unit,
                                                       a_vof,
                                                       a_ebisbox,
                                                       -1,
                                                       cfivs,
                                                       comp);
  }
  else if (m_order == 2) {
    order = EBArith::getExtrapolationStencil(a_sten,
                                             centroid * a_dx,
                                             a_dx * RealVect::Unit,
                                             a_vof,
                                             a_ebisbox,
                                             -1,
                                             cfivs,
                                             comp);
  }
  else {
    MayDay::Abort(
      "CentroidInterpolationStencil::getTaylorExtrapolationStencil - bad order requested. Only first and second order is supported");
  }

  return (order > 0);
}

bool
CentroidInterpolationStencil::getLeastSquaresStencil(VoFStencil&              a_sten,
                                                     const VolIndex&          a_vof,
                                                     const DisjointBoxLayout& a_dbl,
                                                     const ProblemDomain&     a_domain,
                                                     const EBISBox&           a_ebisbox,
                                                     const Box&               a_box,
                                                     const Real&              a_dx,
                                                     const IntVectSet&        a_cfivs)
{

  const int  weightingPower = 0;
  const bool useStartVof    = true;

  a_sten = LeastSquares::getInterpolationStencil(Location::Cell::Centroid,
                                                 Location::Cell::Center,
                                                 LeastSquares::Connectivity::MonotonePath,
                                                 a_vof,
                                                 a_ebisbox,
                                                 a_dx,
                                                 weightingPower,
                                                 m_radius,
                                                 m_order,
                                                 useStartVof);

  return (a_sten.size() > 0);
}

bool
CentroidInterpolationStencil::getPiecewiseLinearStencil(VoFStencil&              a_sten,
                                                        const VolIndex&          a_vof,
                                                        const DisjointBoxLayout& a_dbl,
                                                        const ProblemDomain&     a_domain,
                                                        const EBISBox&           a_ebisbox,
                                                        const Box&               a_box,
                                                        const Real&              a_dx,
                                                        const IntVectSet&        a_cfivs)
{
  a_sten.clear();
  a_sten.add(a_vof, 1.0);

  constexpr Real tolerance = 1.E-8;

  const RealVect centroid = a_ebisbox.centroid(a_vof);
  for (int dir = 0; dir < SpaceDim; dir++) {

    bool     hasLo, hasLower, hasHi, hasHigher;
    VolIndex loVoF, lowerVoF, hiVoF, higherVoF;
    EBArith::getVoFsDir(hasLo, loVoF, hasLower, lowerVoF, a_ebisbox, a_vof, dir, Side::Lo, (IntVectSet*)&a_cfivs);
    EBArith::getVoFsDir(hasHi, hiVoF, hasHigher, higherVoF, a_ebisbox, a_vof, dir, Side::Hi, (IntVectSet*)&a_cfivs);

    if (hasLo || hasHi) {
      if (centroid[dir] > tolerance && hasHi) { // Deriv is to the right
        a_sten.add(hiVoF, centroid[dir]);
        a_sten.add(a_vof, -centroid[dir]);
      }
      else if (centroid[dir] < -tolerance && hasLo) { // Deriv is to the left
        a_sten.add(loVoF, -centroid[dir]);
        a_sten.add(a_vof, centroid[dir]);
      }
      else if (Abs(centroid[dir]) < tolerance) {
        // No deriv in this direction
      }
    }
  }

  return (a_sten.size() > 1);
}

#include <CD_NamespaceFooter.H>
