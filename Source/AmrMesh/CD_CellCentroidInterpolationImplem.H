/* chombo-discharge
 * Copyright © 2025 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CellCentroidInterpolationImplem.H
  @brief  Implementation of CD_CellCentroidInterpolation.H
  @author Robert Marskar
*/

#ifndef CD_CellCentroidInterpolationImplem_H
#define CD_CellCentroidInterpolationImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_CellCentroidInterpolation.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

template <typename T>
void
CellCentroidInterpolation::interpolate(T&               a_centroidData,
                                       const EBCellFAB& a_cellData,
                                       const DataIndex& a_din) const noexcept
{
  CH_TIME("CellCentroidInterpolation::interpolate<T>)");

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
          MayDay::Abort("CD_CellCentroidInterpolation::interpolate(BaseIVFAB) - logic bust");

          break;
        }
        }

        const Real dx = ebisBox.centroid(vof)[dir];

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

inline Real
CellCentroidInterpolation::MinMod(const Real& a_dwl, const Real& a_dwr) const noexcept
{
  Real slope = 0.0;

  if (a_dwl * a_dwr > 0.0) {
    slope = std::abs(a_dwl) < std::abs(a_dwr) ? a_dwl : a_dwr;
  }

  return slope;
}

inline Real
CellCentroidInterpolation::MonotonizedCentral(const Real& a_dwl, const Real& a_dwr) const noexcept
{
  Real slope = 0.0;

  if (a_dwl * a_dwr > 0.0) {
    const Real dwc = a_dwl + a_dwr;
    const Real sgn = Real((dwc > 0.0) - (dwc < 0.0));

    slope = sgn * std::min(0.5 * std::abs(dwc), 2.0 * std::min(std::abs(a_dwl), std::abs(a_dwr)));
  }

  return slope;
}

inline Real
CellCentroidInterpolation::Superbee(const Real& a_dwl, const Real& a_dwr) const noexcept
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

#endif
