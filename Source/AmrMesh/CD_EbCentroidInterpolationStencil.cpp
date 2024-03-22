/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbCentroidInterpolationStencil.cpp
  @brief  Implementation of CD_EbCentroidInterpolationStencil.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_EbCentroidInterpolationStencil.H>
#include <CD_LinearStencil.H>
#include <CD_LeastSquares.H>
#include <CD_NamespaceHeader.H>

EbCentroidInterpolationStencil::EbCentroidInterpolationStencil(const DisjointBoxLayout&        a_dbl,
                                                               const EBISLayout&               a_ebisl,
                                                               const ProblemDomain&            a_domain,
                                                               const Real&                     a_dx,
                                                               const int                       a_order,
                                                               const int                       a_radius,
                                                               const IrregStencil::StencilType a_type)
  : IrregStencil()
{
  CH_TIME("EbCentroidInterpolationStencil::EbCentroidInterpolationStencil");

  this->define(a_dbl, a_ebisl, a_domain, a_dx, a_order, a_radius, a_type);
}

EbCentroidInterpolationStencil::~EbCentroidInterpolationStencil()
{
  CH_TIME("EbCentroidInterpolationStencil::~EbCentroidInterpolationStencil");
}

void
EbCentroidInterpolationStencil::buildStencil(VoFStencil&              a_sten,
                                             const VolIndex&          a_vof,
                                             const DisjointBoxLayout& a_dbl,
                                             const ProblemDomain&     a_domain,
                                             const EBISBox&           a_ebisbox,
                                             const Box&               a_box,
                                             const Real&              a_dx,
                                             const IntVectSet&        a_cfivs)
{
  CH_TIME("EbCentroidInterpolationStencil::buildStencil");

  bool foundStencil = false;

  IntVectSet noCFIVS = IntVectSet();

  // Try to build the preferred stencil.
  switch (m_stencilType) {
  case IrregStencil::StencilType::Linear: {
    const RealVect centroid = a_ebisbox.bndryCentroid(a_vof);
    foundStencil            = LinearStencil::getLinearInterpStencil(a_sten, centroid, a_vof, a_domain, a_ebisbox);
    break;
  }
  case IrregStencil::StencilType::TaylorExtrapolation: {
    foundStencil = this->getTaylorExtrapolationStencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, noCFIVS);
    break;
  }
  case IrregStencil::StencilType::LeastSquares: {
    foundStencil = this->getLeastSquaresInterpolationStencil(a_sten,
                                                             a_vof,
                                                             a_dbl,
                                                             a_domain,
                                                             a_ebisbox,
                                                             a_box,
                                                             a_dx,
                                                             noCFIVS);
    break;
  }
  case IrregStencil::StencilType::PiecewiseLinear: {
    foundStencil = this->getPiecewiseLinearStencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, noCFIVS);
    break;
  }
  default: {
    MayDay::Error("EbCentroidInterpolationStencil::buildStencil - Unsupported stencil type");
    break;
  }
  }

  // If we couldn't find a stencil, try other types in this order
  if (!foundStencil) {
    const RealVect centroid = a_ebisbox.bndryCentroid(a_vof);
    foundStencil            = LinearStencil::getLinearInterpStencil(a_sten, centroid, a_vof, a_domain, a_ebisbox);
  }
  if (!foundStencil) {
    foundStencil = this->getTaylorExtrapolationStencil(a_sten, a_vof, a_dbl, a_domain, a_ebisbox, a_box, a_dx, noCFIVS);
  }
  if (!foundStencil) {
    foundStencil = this->getLeastSquaresInterpolationStencil(a_sten,
                                                             a_vof,
                                                             a_dbl,
                                                             a_domain,
                                                             a_ebisbox,
                                                             a_box,
                                                             a_dx,
                                                             noCFIVS);
  }

  if (!foundStencil) { // Drop to zeroth order.
    a_sten.clear();
    a_sten.add(a_vof, 1.0);
  }
}

bool
EbCentroidInterpolationStencil::getTaylorExtrapolationStencil(VoFStencil&              a_sten,
                                                              const VolIndex&          a_vof,
                                                              const DisjointBoxLayout& a_dbl,
                                                              const ProblemDomain&     a_domain,
                                                              const EBISBox&           a_ebisbox,
                                                              const Box&               a_box,
                                                              const Real&              a_dx,
                                                              const IntVectSet&        a_cfivs)
{
  CH_TIME("EbCentroidInterpolationStencil::getTaylorExtrapolationStencil");

  const int       comp     = 0;
  const RealVect& centroid = a_ebisbox.bndryCentroid(a_vof);
  IntVectSet*     cfivs    = const_cast<IntVectSet*>(&a_cfivs);

  if (m_order == 1) {
    EBArith::getFirstOrderExtrapolationStencil(a_sten,
                                               centroid * a_dx,
                                               a_dx * RealVect::Unit,
                                               a_vof,
                                               a_ebisbox,
                                               -1,
                                               cfivs,
                                               comp);
  }
  else if (m_order == 2) {
    EBArith::getExtrapolationStencil(a_sten, centroid * a_dx, a_dx * RealVect::Unit, a_vof, a_ebisbox, -1, cfivs, comp);
  }
  else {
    MayDay::Error(
      "EbCentroidInterpolationStencil::getTaylorExtrapolationStencil - bad order requested. Only first and second order is supported");
  }

  return a_sten.size() > 0;
}

bool
EbCentroidInterpolationStencil::getLeastSquaresInterpolationStencil(VoFStencil&              a_sten,
                                                                    const VolIndex&          a_vof,
                                                                    const DisjointBoxLayout& a_dbl,
                                                                    const ProblemDomain&     a_domain,
                                                                    const EBISBox&           a_ebisbox,
                                                                    const Box&               a_box,
                                                                    const Real&              a_dx,
                                                                    const IntVectSet&        a_cfivs)
{
  CH_TIME("EbCentroidInterpolationStencil::getLeastSquaresInterpolationStencil");

  const int  weightingPower = 0;
  const bool useStartVof    = true;

  // Have LeastSquares compute an interpolation stencil to specified order. This will compute the pseudo-inverse
  // (using the Moore-Penrose). If this routine returns a stencil with size 0, it failed to find a stencil.
  a_sten = LeastSquares::getInterpolationStencil(Location::Cell::Boundary,
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
EbCentroidInterpolationStencil::getPiecewiseLinearStencil(VoFStencil&              a_sten,
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

  constexpr Real thresh = 1.E-8;

  const RealVect centroid = a_ebisbox.bndryCentroid(a_vof);
  for (int dir = 0; dir < SpaceDim; dir++) {

    bool     hasLo, hasLower, hasHi, hasHigher;
    VolIndex loVoF, lowerVoF, hiVoF, higherVoF;
    EBArith::getVoFsDir(hasLo, loVoF, hasLower, lowerVoF, a_ebisbox, a_vof, dir, Side::Lo, (IntVectSet*)&a_cfivs);
    EBArith::getVoFsDir(hasHi, hiVoF, hasHigher, higherVoF, a_ebisbox, a_vof, dir, Side::Hi, (IntVectSet*)&a_cfivs);

    if (hasLo || hasHi) {
      if (centroid[dir] > thresh && hasHi) { // Deriv is to the right
        a_sten.add(hiVoF, centroid[dir]);
        a_sten.add(a_vof, -centroid[dir]);
      }
      else if (centroid[dir] < -thresh && hasLo) { // Deriv is to the left
        a_sten.add(loVoF, -centroid[dir]);
        a_sten.add(a_vof, centroid[dir]);
      }
    }
  }

  // We will always be able to find this stencil.
  return true;
}

#include <CD_NamespaceFooter.H>
