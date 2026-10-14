/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_IrregularDeposition.H
 * @brief  Declaration of cut-cell deposition strategies.
 * @author Robert Marskar
 */

#ifndef CD_IRREGULARDEPOSITION_H
#define CD_IRREGULARDEPOSITION_H

// Std includes
#include <string>

// Chombo includes
#include <MayDay.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
 * @brief How a deposition scheme treats the cut cells.
 * @details Particle deposition divides by dx^D and never by the volume fraction kappa, so a cut cell holds kappa*n
 * where the fluid solvers hold n. This enum selects what, if anything, is done about that. It is deliberately a single
 * selector rather than a set of independent flags: the alternatives are mutually exclusive answers to one question, and
 * a combination such as Mirror together with Redistribute is a silent 1/kappa error which no assertion further down can
 * catch. Making it unrepresentable is the point of the type.
 * @details This is orthogonal to DepositionType, which selects the deposition kernel (NGP/CIC/TSC), and to
 * CoarseFineDeposition, which selects how the refinement boundary is handled.
 */
enum class IrregularDeposition
{
  /**
   * @brief Deposit as-is, with no cut-cell treatment at all.
   * @details phi is then an extended state that runs into the embedded boundary, and a cut cell holds kappa*n.
   */
  Native,

  /**
   * @brief Put the particle's entire cloud in its own cell when that cell is a cut cell.
   * @details Nothing leaks into the solid, but the cut cell still holds kappa*n and the deposition is no longer a
   * smooth function of the particle position.
   */
  NGP,

  /**
   * @brief Even extension of the density about the embedded boundary: deposit each particle's cloud and the cloud of
   * its image reflected across the boundary.
   * @details Under this option A CUT CELL HOLDS n, NOT kappa*n -- unlike every other option here, and unlike what the
   * rest of the library has historically assumed. That is the whole point (the field is then a true density up to the
   * wall, so a quasi-neutral plasma against a wall produces no spurious charge layer), but it is a change of meaning
   * that consumers of phi must opt into knowingly. Readers that assume kappa*n in a cut cell need attention before they
   * are used together with this option; DataOps::filterSmooth (CD_DataOps.cpp:679) and the arithmeticAverage sites are
   * the known ones.
   * @details Incompatible with DepositionType::NGP, which puts the whole cloud in the particle's own cell -- an image
   * in that same cell would then simply double it. ItoSolver::parseDeposition rejects the combination.
   */
  Mirror,

  /**
   * @brief Hybrid divergence: the cut-cell mass excess deltaM = (1-kappa)*dc is redistributed to the neighbours.
   * @details Conservative, but kappa appears in a denominator, so small kappa is a singular limit.
   */
  Redistribute,

  /**
   * @brief As Redistribute, but the conservative divergence is first blended with the non-conservative one.
   */
  RedistributeBlended
};

/**
 * @brief Map a cut-cell deposition selector string to an IrregularDeposition.
 * @details Valid input options are:
 * 'native' -> Native
 * 'ngp' -> NGP
 * 'mirror' -> Mirror
 * 'redistribute' -> Redistribute
 * 'redistribute_blended' -> RedistributeBlended
 * @param[in] a_str String selector. Must take one of the above strings.
 * @return The corresponding strategy; aborts on an unknown selector.
 */
inline IrregularDeposition
irregularDepositionFromString(const std::string& a_str) noexcept;

#include <CD_NamespaceFooter.H>

#include <CD_IrregularDepositionImplem.H>

#endif
