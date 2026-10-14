/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ParticleMergeTagger.H
 * @brief  Declaration of a cell tagger that refines only where the user says to
 * @author Robert Marskar
 */

#ifndef CD_PARTICLEMERGETAGGER_H
#define CD_PARTICLEMERGETAGGER_H

// Our includes
#include <CD_AmrMesh.H>
#include <CD_CellTagger.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Cell tagger that refines only the manual refinement boxes given in the input script.
 * @details CellTagger already knows how to read manual refinement boxes -- parseRefinementBoxes() reads
 * ParticleMergeTagger.num_ref_boxes and the ref_boxN_lo/hi/lvl triplets, and
 * getManualRefinementLevel() answers, for a position, which level the user asked it to be refined to.
 * This class does nothing but apply that answer, so the refined region is exactly what the input script
 * names and nothing else.
 *
 * That is the point for a merge study: a refinement box does not depend on the solution, so the
 * coarse-fine boundary is a surface at a position the user chose and it does not move as the run
 * proceeds. A merge's sub-cell artifacts are locked to the local cell width, so the two sides of that
 * boundary alias on lattices whose periods differ by the refinement ratio; holding the boundary still is
 * what makes the two sides comparable from one step to the next.
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
   * @param[in] a_amr AMR mesh, which supplies the grids the tags are laid out on.
   */
  virtual void
  define(const RefCountedPtr<AmrMesh>& a_amr) noexcept;

  /**
   * @brief Parse the manual refinement boxes and the tag buffer.
   */
  void
  parseOptions() override;

  /**
   * @brief Regrid the tagger. Nothing to do -- the tagged region does not depend on the grids or the
   * solution.
   */
  void
  regrid() override;

  /**
   * @brief Tag every cell the manual refinement boxes ask to be refined beyond its own level.
   * @param[in,out] a_tags Per-level, per-patch tags to fill.
   * @return True, always: the tags are recomputed from the (fixed) boxes on every call.
   */
  bool
  tagCells(EBAMRTags& a_tags) override;

  /**
   * @brief Turn the per-patch tags into the per-level sets that AmrMesh::regridAmr consumes.
   * @details CellTagger fills a per-patch LayoutData<DenseIntVectSet> while AmrMesh takes one IntVectSet
   * per level, so something has to bridge the two. Driver does this internally; an application that
   * drives AmrMesh directly has to do it itself.
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
};

#include <CD_NamespaceFooter.H>

#endif
