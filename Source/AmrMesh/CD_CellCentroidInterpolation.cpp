/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CellCentroidInterpolation.cpp
 * @brief  Implementation of CD_CellCentroidInterpolation.H
 * @author Robert Marskar
 */

// Chombo includes
#include <CH_Timer.H>
#include <EBArith.H>

// Our includes
#include <CD_CellCentroidInterpolation.H>
#include <CD_LeastSquares.H>
#include <CD_LinearStencil.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

CellCentroidInterpolation::CellCentroidInterpolation() noexcept : m_isDefined(false), m_dx(-1.0)
{
  CH_TIME("CellCentroidInterpolation:CellCentroidInterpolation(weak)");
}

CellCentroidInterpolation::CellCentroidInterpolation(const EBLevelGrid& a_eblg,
                                                     const Real&        a_dx,
                                                     const Type&        a_interpolationType) noexcept
  : CellCentroidInterpolation()
{
  CH_TIME("CellCentroidInterpolation::CellCentroidInterpolation(weak)");

  this->define(a_eblg, a_dx, a_interpolationType);
}

CellCentroidInterpolation::~CellCentroidInterpolation() noexcept
{
  CH_TIME("CellCentroidInterpolation::~CellCentroidInterpolation");
}

void
CellCentroidInterpolation::define(const EBLevelGrid& a_eblg, const Real& a_dx, const Type& a_interpolationType) noexcept
{
  CH_TIME("CellCentroidInterpolation::define");

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

  // Stencil-based interpolations use one stencil per cut cell. The slope-limited ones store the cell's neighbors in
  // 2*SpaceDim components instead.
  const bool slopeLimited = (m_interpolationType == Type::MinMod) ||
                            (m_interpolationType == Type::MonotonizedCentral) ||
                            (m_interpolationType == Type::Superbee);

  const int numStencilComps = slopeLimited ? 2 * SpaceDim : 1;

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
    stencils.define(irregIVS, ebGraph, numStencilComps);

    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex& vof     = vofit();
      VoFStencil&     stencil = stencils(vof, 0);

      auto defaultStencil = [&]() -> void {
        stencil.clear();
        stencil.add(vof, 1.0);
      };

      // For the slope-limited interpolations, components 2*dir and 2*dir + 1 hold the neighbor across the low and high
      // face in direction dir, if there is exactly one such face and the cell is not on the domain boundary. The
      // neighbors are found once here, so that interpolate() does not query the graph for faces.
      auto slopeNeighbors = [&]() -> void {
        const IntVect iv = vof.gridIndex();

        for (int dir = 0; dir < SpaceDim; dir++) {
          for (SideIterator sit; sit.ok(); ++sit) {
            const Side::LoHiSide side     = sit();
            const bool           onDomain = (side == Side::Lo) ? (iv[dir] == domain.domainBox().smallEnd(dir))
                                                               : (iv[dir] == domain.domainBox().bigEnd(dir));

            VoFStencil& neighbor = stencils(vof, 2 * dir + ((side == Side::Lo) ? 0 : 1));

            neighbor.clear();

            if (ebisBox.numFaces(vof, dir, side) == 1 && !onDomain) {
              const Vector<FaceIndex> faces = ebisBox.getFaces(vof, dir, side);

              neighbor.add(faces[0].getVoF(side), 1.0);
            }
          }
        }
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
        slopeNeighbors();

        break;
      }
      case Type::MonotonizedCentral: {
        slopeNeighbors();

        break;
      }
      case Type::Superbee: {
        slopeNeighbors();

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
CellCentroidInterpolation::getLinearStencil(VoFStencil&          a_stencil,
                                            const VolIndex&      a_vof,
                                            const EBISBox&       a_ebisBox,
                                            const ProblemDomain& a_domain) const noexcept
{
  const bool foundStencil = LinearStencil::getLinearInterpStencil(a_stencil,
                                                                  a_ebisBox.centroid(a_vof),
                                                                  a_vof,
                                                                  a_domain,
                                                                  a_ebisBox);

  return foundStencil;
}

bool
CellCentroidInterpolation::getTaylorExtrapolationStencil(VoFStencil&     a_stencil,
                                                         const VolIndex& a_vof,
                                                         const EBISBox&  a_ebisBox,
                                                         const ProblemDomain& /*a_domain*/) const noexcept
{
  const int order = EBArith::getFirstOrderExtrapolationStencil(a_stencil,
                                                               a_ebisBox.centroid(a_vof) * m_dx,
                                                               m_dx * RealVect::Unit,
                                                               a_vof,
                                                               a_ebisBox,
                                                               -1,
                                                               nullptr,
                                                               0);

  return (order > 0);
}

bool
CellCentroidInterpolation::getLeastSquaresStencil(VoFStencil&     a_stencil,
                                                  const VolIndex& a_vof,
                                                  const EBISBox&  a_ebisBox,
                                                  const ProblemDomain& /*a_domain*/) const noexcept
{
  const int radius = 1;
  const int order  = 1;
  const int weight = 0;

  const bool useStartVoF = true;

  a_stencil = LeastSquares::getInterpolationStencil(Location::Cell::Centroid,
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
CellCentroidInterpolation::getPiecewiseLinearStencil(VoFStencil&     a_stencil,
                                                     const VolIndex& a_vof,
                                                     const EBISBox&  a_ebisBox,
                                                     const ProblemDomain& /*a_domain*/) const noexcept
{
  a_stencil.clear();
  a_stencil.add(a_vof, 1.0);

  const Real     tolerance = 1.E-8;
  const RealVect centroid  = a_ebisBox.centroid(a_vof);

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
CellCentroidInterpolation::interpolate(LevelData<BaseIVFAB<Real>>& a_centroidData,
                                       const LevelData<EBCellFAB>& a_cellData) const noexcept
{
  CH_TIME("CellCentroidInterpolation::interpolate(LD<BaseIVFAB<Real>>, LD<EBCellFAB>)");

  CH_assert(m_isDefined);
  CH_assert(a_centroidData.isDefined());
  CH_assert(a_cellData.isDefined());
  CH_assert(a_centroidData.disjointBoxLayout() == m_eblg.getDBL());
  CH_assert(a_cellData.disjointBoxLayout() == m_eblg.getDBL());
  CH_assert(a_centroidData.nComp() == a_cellData.nComp());

  const DisjointBoxLayout& dbl   = m_eblg.getDBL();
  const EBISLayout&        ebisl = m_eblg.getEBISL();
  const DataIterator&      dit   = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    this->interpolate<BaseIVFAB<Real>>(a_centroidData[din], a_cellData[din], din);
  }
}

void
CellCentroidInterpolation::interpolate(LevelData<EBCellFAB>&       a_centroidData,
                                       const LevelData<EBCellFAB>& a_cellData) const noexcept
{
  CH_TIME("CellCentroidInterpolation::interpolate(LD<EBCellFAB>, LD<EBCellFAB>)");

  CH_assert(m_isDefined);
  CH_assert(a_centroidData.isDefined());
  CH_assert(a_cellData.isDefined());
  CH_assert(a_centroidData.disjointBoxLayout() == m_eblg.getDBL());
  CH_assert(a_cellData.disjointBoxLayout() == m_eblg.getDBL());
  CH_assert(a_centroidData.nComp() == a_cellData.nComp());

  const DisjointBoxLayout& dbl   = m_eblg.getDBL();
  const EBISLayout&        ebisl = m_eblg.getEBISL();
  const DataIterator&      dit   = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    a_centroidData[din].copy(a_cellData[din]);

    this->interpolate<EBCellFAB>(a_centroidData[din], a_cellData[din], din);
  }
}

void
CellCentroidInterpolation::interpolate(LevelData<EBCellFAB>& a_data) const noexcept
{
  CH_TIME("CellCentroidInterpolation::interpolate(LD<EBCellFAB>)");

  CH_assert(m_isDefined);
  CH_assert(a_data.isDefined());
  CH_assert(a_data.disjointBoxLayout() == m_eblg.getDBL());

  const DisjointBoxLayout& dbl   = m_eblg.getDBL();
  const EBISLayout&        ebisl = m_eblg.getEBISL();
  const DataIterator&      dit   = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    EBCellFAB tmp;
    tmp.clone(a_data[din]);

    this->interpolate<EBCellFAB>(a_data[din], tmp, din);
  }
}

#include <CD_NamespaceFooter.H>
