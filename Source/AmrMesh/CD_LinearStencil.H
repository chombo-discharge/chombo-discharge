/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_LinearStencil.H
  @brief  Routines for computing a bi/tri-linear interpolation stencil in cut-cells. 
  @author Robert Marskar
*/

#ifndef CD_LINEARSTENCIL_H
#define CD_LINEARSTENCIL_H

// Chombo includes
#include <Stencil.H>
#include <ProblemDomain.H>
#include <EBISBox.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Class which computes a bi/tri-linear stencil in cut-cells. This class makes appropriate modifications if one of the
  cells don't exist, and even if there are multivalued cells present. 
  @note This routine assumes that the data is fundamentally cell-centered. 
*/
class LinearStencil
{
public:
  /**
    @brief Get a bi(tri)linear interpolation stencil in 2D(3D)
    @param[out] a_stencil  Interpolation stencil
    @param[in]  a_pos      Relative position inside cell
    @param[in]  a_vof      Cut-cell
    @param[in]  a_domain   Domain
    @return Linear interp stencil
    @param[in] a_ebisbox Ebisbox
  */
  static bool
  getLinearInterpStencil(VoFStencil&          a_stencil,
                         const RealVect&      a_pos,
                         const VolIndex&      a_vof,
                         const ProblemDomain& a_domain,
                         const EBISBox&       a_ebisbox);

protected:
  /**
    @brief Tolerance for when we switch between 1D/2D/3D stencils
  */
  static constexpr Real tolerance = 1.E-6;

  /**
    @brief Compute a linear interpolation stencil. 
    @details This only does interpolation along the input coordinate. The other
    coordinates are ignored. 
    @param[out] a_stencil   Interpolation stencil
    @param[in]  a_pos       Relative position inside cell
    @param[in]  a_vof       Cut-cell
    @param[in]  a_domain    Domain
    @param[in]  a_interpDir Interpolation coordinate
    @return Computed interp stencil1 d
    @param[in] a_ebisbox Ebisbox
  */
  static bool
  computeInterpStencil1D(VoFStencil&          a_stencil,
                         const RealVect&      a_pos,
                         const VolIndex&      a_vof,
                         const ProblemDomain& a_domain,
                         const EBISBox&       a_ebisbox,
                         const int            a_interpDir);
  /**
    @brief Compute a bilinear interpolation stencil
    @details This ignores interpolation along the input coordinate. 
    @param[out] a_stencil   Interpolation stencil
    @param[in]  a_pos  Relative position inside cell
    @param[in]  a_vof       Cut-cell
    @param[in]  a_domain    Domain
    @return Computed interp stencil2 d
    @param[in] a_ebisbox Ebisbox
    @param[in] a_noInterpDir No interp dir
  */
  static bool
  computeInterpStencil2D(VoFStencil&          a_stencil,
                         const RealVect&      a_pos,
                         const VolIndex&      a_vof,
                         const ProblemDomain& a_domain,
                         const EBISBox&       a_ebisbox,
                         const int            a_noInterpDir = 2);

#if CH_SPACEDIM == 3
  /**
    @brief Compute a trilinear interpolation stencil
    @param[out] a_stencil  Interpolation stencil
    @param[in]  a_pos      Relative position inside cell
    @param[in]  a_vof      Cut-cell
    @param[in]  a_domain   Domain
    @param[in]  a_ebisBox  EBISBox
  */
  static bool
  computeInterpStencil3D(VoFStencil&          a_stencil,
                         const RealVect&      a_pos,
                         const VolIndex&      a_vof,
                         const ProblemDomain& a_domain,
                         const EBISBox&       a_ebisbox);
#endif
};

#include <CD_NamespaceFooter.H>

#endif
