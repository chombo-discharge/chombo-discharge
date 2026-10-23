/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_RefinementBoxTagger.H
 * @brief  Declaration of a cell tagger that refines user-specified boxes.
 * @author Robert Marskar
 */

#ifndef CD_REFINEMENTBOXTAGGER_H
#define CD_REFINEMENTBOXTAGGER_H

// Std includes
#include <string>

// Chombo includes
#include <RefCountedPtr.H>

// Our includes
#include <CD_CellTagger.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Cell tagger that only refines the boxes given as RefinementBoxTagger.ref_box* in the input file.
 * @details This gives static refinement where neither geometric refinement nor a physics tagger applies, e.g. on a
 * geometry without embedded boundaries. A cell is tagged on level lvl if its center lies in a refinement box that asks
 * for a level finer than lvl. The box edges are included (CellTagger::getManualRefinementLevel), but an edge that
 * coincides exactly with a cell center is decided by floating-point rounding, so place the edges between cell centers.
 * The boxes, and the tag buffer, are parsed as in CellTagger.
 */
class RefinementBoxTagger : public CellTagger
{
public:
  /**
   * @brief Disallowed default constructor.
   */
  RefinementBoxTagger() = delete;

  /**
   * @brief Constructor.
   * @param[in] a_amr   AMR mesh
   * @param[in] a_realm Realm whose grids are tagged
   */
  RefinementBoxTagger(const RefCountedPtr<AmrMesh>& a_amr, const std::string& a_realm) noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~RefinementBoxTagger() noexcept;

  /**
   * @brief Parse the verbosity, the refinement boxes and the buffer.
   */
  virtual void
  parseOptions() override;

  /**
   * @brief Nothing to regrid: the tagger holds no mesh data.
   */
  virtual void
  regrid() override;

  /**
   * @brief Tag every cell whose center lies in a refinement box that asks for a finer level.
   * @param[in,out] a_tags Cell tags on each level
   * @return True if any cell was tagged.
   */
  virtual bool
  tagCells(EBAMRTags& a_tags) override;

protected:
  /**
   * @brief AMR mesh.
   */
  RefCountedPtr<AmrMesh> m_amr;

  /**
   * @brief Realm whose grids are tagged.
   */
  std::string m_realm;
};

#include <CD_NamespaceFooter.H>

#endif
