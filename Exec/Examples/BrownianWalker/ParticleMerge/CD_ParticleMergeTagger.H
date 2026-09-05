/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ParticleMergeTagger.H
 * @brief  Declaration of a cell tagger that refines a fixed, manually specified slab
 * @author Robert Marskar
 */

#ifndef CD_PARTICLEMERGETAGGER_H
#define CD_PARTICLEMERGETAGGER_H

// Our includes
#include <CD_AmrMesh.H>
#include <CD_CellTagger.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Cell tagger that refines a fixed slab of the domain, specified in the input script.
 * @details The slab spans the domain in every direction but one, and is bounded in that one direction
 * by two physical coordinates. Because the refined region is fixed and known, the coarse-fine boundary
 * it produces is a flat plane at a position the user chose, which is what makes it useful for measuring
 * what a merge does on either side of a refinement boundary: no feedback from the solution, no drifting
 * boundary, and the same geometry on every step.
 */
class ParticleMergeTagger : public CellTagger
{
public:
  /**
   * @brief Default constructor. Leaves the object undefined; call define().
   */
  ParticleMergeTagger() noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~ParticleMergeTagger() noexcept;

  /**
   * @brief Define the tagger.
   * @param[in] a_amr AMR mesh, used for the grids the tags are laid out on.
   */
  virtual void
  define(const RefCountedPtr<AmrMesh>& a_amr) noexcept;

  /**
   * @brief Parse the slab bounds and the direction it is bounded in.
   * @details Reads ParticleMergeTagger.refine_dir, .refine_lo and .refine_hi.
   */
  void
  parseOptions() override;

  /**
   * @brief Regrid the tagger. Nothing to do -- the tagged region does not depend on the grids.
   */
  void
  regrid() override;

  /**
   * @brief Tag every cell whose centre falls inside the slab.
   * @param[in,out] a_tags Per-level, per-patch tags to fill.
   * @return True, always: the tags are recomputed from the (fixed) slab every call.
   */
  bool
  tagCells(EBAMRTags& a_tags) override;

  /**
   * @brief Turn the per-patch tags into the per-level sets that AmrMesh::regridAmr consumes.
   * @details AmrMesh takes one IntVectSet per level while CellTagger fills a per-patch
   * LayoutData<DenseIntVectSet>, so something has to bridge the two. Driver does this internally; an
   * application that drives AmrMesh directly, as this example does, has to do it itself.
   * @param[out] a_levelTags Per-level tags.
   * @param[in]  a_tags      Per-level, per-patch tags as filled by tagCells().
   */
  void
  gatherTags(Vector<IntVectSet>& a_levelTags, const EBAMRTags& a_tags) const noexcept;

protected:
  /**
   * @brief AMR mesh.
   */
  RefCountedPtr<AmrMesh> m_amr;

  /**
   * @brief Direction the refined slab is bounded in. It spans the domain in every other direction.
   */
  int m_refineDir;

  /**
   * @brief Lower bound of the refined slab, as a physical coordinate along m_refineDir.
   */
  Real m_refineLo;

  /**
   * @brief Upper bound of the refined slab, as a physical coordinate along m_refineDir.
   */
  Real m_refineHi;
};

#include <CD_NamespaceFooter.H>

#endif
