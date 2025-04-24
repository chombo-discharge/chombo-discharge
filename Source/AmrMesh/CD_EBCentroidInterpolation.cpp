/* chombo-discharge
 * Copyright © 2025 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBCentroidInterpolation.cpp
  @brief  Implementation of CD_EBCentroidInterpolation.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <EBArith.H>

// Our includes
#include <CD_EBCentroidInterpolation.H>
#include <CD_LeastSquares.H>
#include <CD_LinearStencil.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBCentroidInterpolation::EBCentroidInterpolation() noexcept
{
  CH_TIME("EBCentroidInterpolation:EBCentroidInterpolation(weak)");

  m_dx        = -1.0;
  m_isDefined = false;
}

EBCentroidInterpolation::EBCentroidInterpolation(const EBLevelGrid& a_eblg,
                                                 const Real&        a_dx,
                                                 const Type&        a_interpolationType) noexcept
  : EBCentroidInterpolation()
{
  CH_TIME("EBCentroidInterpolation::EBCentroidInterpolation(weak)");

  this->define(a_eblg, a_dx, a_interpolationType);
}

EBCentroidInterpolation::~EBCentroidInterpolation() noexcept
{

  CH_TIME("EBCentroidInterpolation::~EBCentroidInterpolation");
}

void
EBCentroidInterpolation::define(const EBLevelGrid& a_eblg, const Real& a_dx, const Type& a_interpolationType) noexcept
{
  CH_TIME("EBCentroidInterpolation::define");

  CH_assert(a_eblg.isDefined());
  CH_assert(a_dx >= 0.0);

  m_eblg              = a_eblg;
  m_dx                = a_dx;
  m_interpolationType = a_interpolationType;

  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const EBISLayout&        ebisl  = m_eblg.getEBISL();
  const ProblemDomain&     domain = m_eblg.getDomain();
  const DataIterator&      dit    = dbl.dataIterator();

  m_vofIterator.define(dbl);
  m_interpStencils.define(dbl);

  const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {

    const DataIndex& din = dit[mybox];

    const EBISBox&   ebisBox  = ebisl[din];
    const Box&       cellBox  = dbl[din];
    const EBGraph&   ebGraph  = ebisBox.getEBGraph();
    const IntVectSet irregIVS = ebisBox.getIrregIVS(cellBox);

    VoFIterator&           vofit    = m_vofIterator[din];
    BaseIVFAB<VoFStencil>& stencils = m_interpStencils[din];

    vofit.define(irregIVS, ebGraph);
    stencils.define(irregIVS, ebGraph, 1);

    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex& vof     = vofit();
      VoFStencil&     stencil = stencils(vof, 0);

      auto defaultStencil = [&]() -> void {
        stencil.clear();
        stencil.add(vof, 1.0);
      };

      switch (m_interpolationType) {
      case Type::Constant: {
        stencil.add(vof, 1.0);

        break;
      }
      case Type::Linear: {
        const bool foundStencil = this->getLinearStencil(stencil, vof, ebisBox, domain);

        if (!foundStencil) {
          defaultStencil();
        }

        break;
      }
      case Type::Taylor: {
        const bool foundStencil = this->getTaylorExtrapolationStencil(stencil, vof, ebisBox, domain);

        if (!foundStencil) {
          defaultStencil();
        }

        break;
      }
      case Type::LeastSquares: {
        const bool foundStencil = this->getLeastSquaresStencil(stencil, vof, ebisBox, domain);

        if (!foundStencil) {
          defaultStencil();
        }

        break;
      }
      case Type::PiecewiseLinear: {
        const bool foundStencil = this->getPiecewiseLinearStencil(stencil, vof, ebisBox, domain);

        if (!foundStencil) {
          defaultStencil();
        }

        break;
      }
      case Type::MinMod: {
        stencil.clear(); // No need for a stencil.

        break;
      }
      case Type::MonotonizedCentral: {
        stencil.clear(); // No need for a stencil.

        break;
      }
      case Type::Superbee: {
        stencil.clear(); // No need for a stencil.

        break;
      }
      default: {
        defaultStencil();

        break;
      }
      }
    }
  }

  m_isDefined = true;
}

bool
EBCentroidInterpolation::getLinearStencil(VoFStencil&          a_stencil,
                                          const VolIndex&      a_vof,
                                          const EBISBox&       a_ebisBox,
                                          const ProblemDomain& a_domain) const noexcept
{
  CH_TIME("EBCentroidInterpolation::getLinearStencil");

  const bool foundStencil = LinearStencil::getLinearInterpStencil(a_stencil,
                                                                  a_ebisBox.bndryCentroid(a_vof),
                                                                  a_vof,
                                                                  a_domain,
                                                                  a_ebisBox);

  return foundStencil;
}

bool
EBCentroidInterpolation::getTaylorExtrapolationStencil(VoFStencil&          a_stencil,
                                                       const VolIndex&      a_vof,
                                                       const EBISBox&       a_ebisBox,
                                                       const ProblemDomain& a_domain) const noexcept
{
  CH_TIME("EBCentroidInterpolation::getTaylorExtrapolationStencil");

  const int order = EBArith::getFirstOrderExtrapolationStencil(a_stencil,
                                                               a_ebisBox.bndryCentroid(a_vof) * m_dx,
                                                               m_dx * RealVect::Unit,
                                                               a_vof,
                                                               a_ebisBox,
                                                               -1,
                                                               nullptr,
                                                               0);

  return (order > 0);
}

bool
EBCentroidInterpolation::getLeastSquaresStencil(VoFStencil&          a_stencil,
                                                const VolIndex&      a_vof,
                                                const EBISBox&       a_ebisBox,
                                                const ProblemDomain& a_domain) const noexcept
{
  CH_TIME("EBCentroidInterpolation::getLeastSquaresStencil");

  const int radius = 1;
  const int order  = 1;
  const int weight = 0;

  const bool useStartVoF = true;

  a_stencil = LeastSquares::getInterpolationStencil(Location::Cell::Boundary,
                                                    Location::Cell::Center,
                                                    LeastSquares::Connectivity::MonotonePath,
                                                    a_vof,
                                                    a_ebisBox,
                                                    m_dx,
                                                    weight,
                                                    radius,
                                                    order,
                                                    useStartVoF);

  return (a_stencil.size() > 0);
}

bool
EBCentroidInterpolation::getPiecewiseLinearStencil(VoFStencil&          a_stencil,
                                                   const VolIndex&      a_vof,
                                                   const EBISBox&       a_ebisBox,
                                                   const ProblemDomain& a_domain) const noexcept
{
  CH_TIME("EBCentroidInterpolation::getPiecewiseLinearStencil");

  a_stencil.clear();
  a_stencil.add(a_vof, 1.0);

  const Real     tolerance = 1.E-8;
  const RealVect centroid  = a_ebisBox.bndryCentroid(a_vof);

  for (int dir = 0; dir < SpaceDim; dir++) {

    bool hasLo     = false;
    bool hasLower  = false;
    bool hasHi     = false;
    bool hasHigher = false;

    VolIndex loVoF;
    VolIndex lowerVoF;
    VolIndex hiVoF;
    VolIndex higherVoF;

    EBArith::getVoFsDir(hasLo, loVoF, hasLower, lowerVoF, a_ebisBox, a_vof, dir, Side::Lo, nullptr);
    EBArith::getVoFsDir(hasHi, hiVoF, hasHigher, higherVoF, a_ebisBox, a_vof, dir, Side::Hi, nullptr);

    if (hasLo || hasHi) {
      if (centroid[dir] > tolerance && hasHi) { // Deriv is to the right
        a_stencil.add(hiVoF, centroid[dir]);
        a_stencil.add(a_vof, -centroid[dir]);
      }
      else if (centroid[dir] < -tolerance && hasLo) { // Deriv is to the left
        a_stencil.add(loVoF, -centroid[dir]);
        a_stencil.add(a_vof, centroid[dir]);
      }
      else if (std::abs(centroid[dir]) < tolerance) {
        // No deriv in this direction
      }
    }
  }

  return (a_stencil.size() > 1);
}

void
EBCentroidInterpolation::interpolate(LevelData<BaseIVFAB<Real>>& a_centroidData,
                                     const LevelData<EBCellFAB>& a_cellData) const noexcept
{
  CH_TIME("EBCentroidInterpolation::interpolate(LD<BaseIVFAB<Real>>, LD<EBCellFAB>)");

  CH_assert(m_isDefined);
  CH_assert(a_centroidData.isDefined());
  CH_assert(a_cellData.isDefined());
  CH_assert(a_centroidData.disjointBoxLayout() == m_eblg.getDBL());
  CH_assert(a_cellData.disjointBoxLayout() == m_eblg.getDBL());
  CH_assert(a_centroidData.nComp() == a_cellData.nComp());

  const DisjointBoxLayout& dbl       = m_eblg.getDBL();
  const ProblemDomain&     domain    = m_eblg.getDomain();
  const EBISLayout&        ebisl     = m_eblg.getEBISL();
  const Box&               domainBox = domain.domainBox();
  const DataIterator&      dit       = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    this->interpolate(a_centroidData[din], a_cellData[din], din);
  }
}

void
EBCentroidInterpolation::interpolate(BaseIVFAB<Real>& a_centroidData,
                                     const EBCellFAB& a_cellData,
                                     const DataIndex& a_din) const noexcept
{
  CH_TIME("EBCentroidInterpolation::interpolate<T>)");

  CH_assert(m_isDefined);
  CH_assert(a_centroidData.isDefined());
  CH_assert(a_cellData.isDefined());
  CH_assert(a_centroidData.nComp() == a_cellData.nComp());

  const DisjointBoxLayout& dbl       = m_eblg.getDBL();
  const ProblemDomain&     domain    = m_eblg.getDomain();
  const EBISLayout&        ebisl     = m_eblg.getEBISL();
  const Box&               domainBox = domain.domainBox();

  const EBISBox&               ebisBox  = ebisl[a_din];
  const Box&                   cellBox  = dbl[a_din];
  const BaseIVFAB<VoFStencil>& stencils = m_interpStencils[a_din];

  const int nComp = a_cellData.nComp();

  for (int comp = 0; comp < nComp; comp++) {

    // This is the kernel that is used when the interpolation is expressable as a stencil.
    auto stencilKernel = [&](const VolIndex& vof) -> void {
      a_centroidData(vof, comp) = 0.0;

      const VoFStencil& stencil = stencils(vof, 0);
      for (int i = 0; i < stencil.size(); i++) {
        const VolIndex& ivof    = stencil.vof(i);
        const Real&     iweight = stencil.weight(i);

        a_centroidData(vof, comp) += iweight * a_cellData(ivof, comp);
      }
    };

    // Kernel used when interpolation is done with a slope limiter.
    auto slopeKernel = [&](const VolIndex& vof) -> void {
      const IntVect iv = vof.gridIndex();

      a_centroidData(vof, comp) = a_cellData(vof, comp);

      for (int dir = 0; dir < SpaceDim; dir++) {
        Real slope = 0.0;

        const bool onLoSide = (iv[dir] == domainBox.smallEnd(dir));
        const bool onHiSide = (iv[dir] == domainBox.bigEnd(dir));

        const bool hasFacesLeft  = (ebisBox.numFaces(vof, dir, Side::Lo) == 1) && !onLoSide;
        const bool hasFacesRight = (ebisBox.numFaces(vof, dir, Side::Hi) == 1) && !onHiSide;

        Vector<FaceIndex> facesLeft;
        Vector<FaceIndex> facesRight;

        VolIndex vofLeft;
        VolIndex vofRight;

        Real dwl = 0.0;
        Real dwr = 0.0;

        // Compute left and right slope
        if (hasFacesLeft) {
          facesLeft = ebisBox.getFaces(vof, dir, Side::Lo);
          vofLeft   = facesLeft[0].getVoF(Side::Lo);
          dwl       = a_cellData(vof, comp) - a_cellData(vofLeft, comp);
        }
        if (hasFacesRight) {
          facesRight = ebisBox.getFaces(vof, dir, Side::Hi);
          vofRight   = facesRight[0].getVoF(Side::Hi);
          dwr        = a_cellData(vofRight, comp) - a_cellData(vof, comp);
        }

        if (!hasFacesLeft && hasFacesRight) {
          dwl = dwr;
        }
        else if (hasFacesLeft && !hasFacesRight) {
          dwr = dwl;
        }

        // Limit the slopes.
        switch (m_interpolationType) {
        case Type::MinMod: {
          slope = this->MinMod(dwl, dwr);

          break;
        }
        case Type::MonotonizedCentral: {
          slope = this->MonotonizedCentral(dwl, dwr);

          break;
        }
        case Type::Superbee: {
          slope = this->Superbee(dwl, dwr);

          break;
        }
        default: {
          MayDay::Abort("CD_EBCentroidInterpolation::interpolate(BaseIVFAB) - logic bust");

          break;
        }
        }

        const Real dx = ebisBox.bndryCentroid(vof)[dir];

        a_centroidData(vof, comp) += slope * dx;
      }
    };

    switch (m_interpolationType) {
    case Type::MinMod: {
      BoxLoops::loop(m_vofIterator[a_din], slopeKernel);

      break;
    }
    case Type::MonotonizedCentral: {
      BoxLoops::loop(m_vofIterator[a_din], slopeKernel);

      break;
    }
    case Type::Superbee: {
      BoxLoops::loop(m_vofIterator[a_din], slopeKernel);

      break;
    }
    default: {
      BoxLoops::loop(m_vofIterator[a_din], stencilKernel);

      break;
    }
    }
  }
}

Real
EBCentroidInterpolation::MinMod(const Real& a_dwl, const Real& a_dwr) const noexcept
{
  Real slope = 0.0;

  if (a_dwl * a_dwr > 0.0) {
    slope = std::abs(a_dwl) < std::abs(a_dwr) ? a_dwl : a_dwr;
  }

  return slope;
}

Real
EBCentroidInterpolation::MonotonizedCentral(const Real& a_dwl, const Real& a_dwr) const noexcept
{
  Real slope = 0.0;

  if (a_dwl * a_dwr > 0.0) {
    const Real dwc = a_dwl + a_dwr;
    const Real sgn = Real((dwc > 0.0) - (dwc < 0.0));

    slope = sgn * std::min(0.5 * std::abs(dwc), 2.0 * std::min(std::abs(a_dwl), std::abs(a_dwr)));
  }

  return slope;
}

Real
EBCentroidInterpolation::Superbee(const Real& a_dwl, const Real& a_dwr) const noexcept
{
  Real slope = 0.0;

  if (a_dwl * a_dwr > 0.0) {
    const Real s1 = this->MinMod(a_dwl, 2 * a_dwr);
    const Real s2 = this->MinMod(a_dwr, 2 * a_dwl);

    if (s1 * s2 > 0.0) {
      slope = std::abs(s1) > std::abs(s2) ? s1 : s2;
    }
  }

  return slope;
}

#include <CD_NamespaceFooter.H>
