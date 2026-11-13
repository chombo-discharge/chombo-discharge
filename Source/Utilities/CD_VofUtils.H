/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_VofUtils.H
 * @brief  Various functions for getting Vofs near cut-cells
 * @author Robert Marskar
 */

#ifndef CD_VOFUTILS_H
#define CD_VOFUTILS_H

// Chombo includes
#include <Stencils.H>
#include <EBISBox.H>
#include <EBISLayout.H>
#include <FaceIndex.H>
#include <IntVectSet.H>
#include <RealVect.H>
#include <ProblemDomain.H>
#include <VolIndex.H>
#include <LoHiSide.H>
#include <DenseIntVectSet.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Static class which contains some routines for fetching VoFs using various algorithms. Very useful in
 * combination with least squares.
 */
class VofUtils
{
public:
  /**
   * @brief Class for specifying connectivity of Vofs.
   * @details MontonePath:     Vofs must be reachable with monotone path.
   * @details SimplyConnected: All vofs must be connected
   * @details All:             All vofs allowed
   */
  enum class Connectivity
  {
    MonotonePath,
    SimplyConnected,
    All
  };

  /**
   * @brief Neighborhood type
   * @param[in] Quadrant.   Only fetch vofs in a specified quadrant.
   * @param[in] Radius.     Fetch all vofs in a specified radius.
   * @param[in] SemiCircle. Fetch all vofs in a half-circle with specified radius.
   */
  enum class Neighborhood
  {
    Quadrant,
    Radius,
    SemiCircle
  };

  /**
   * @brief Disallowed
   */
  VofUtils() = delete;

  /**
   * @brief Fetch vofs in a specified radius and with a specified connectivitiy between them.
   * @param[in] a_startVof     Starting Vof.
   * @param[in] a_ebisbox      EBISBox.
   * @param[in] a_radius       Radius.
   * @param[in] a_connectivity Connectivity specification (e.g. monotone path or simplyconnected)
   * @param[in] a_addStartVof  If false, a_startVof will not be included in the list.
   * @return Vofs in radius
   */
  static Vector<VolIndex>
  getVofsInRadius(const VolIndex&    a_startVof,
                  const EBISBox&     a_ebisbox,
                  const int          a_radius,
                  const Connectivity a_connectivity,
                  const bool         a_addStartVof);

  /**
   * @brief Fetch vofs in a specified radius and with a specified connectivitiy between them.
   * @param[in] a_startVof     Starting Vof.
   * @param[in] a_ebisbox      EBISBox.
   * @param[in] a_radius       Radius.
   * @param[in] a_validCells   Valid cells to include
   * @param[in] a_connectivity Connectivity specification (e.g. monotone path or simplyconnected)
   * @param[in] a_addStartVof  If false, a_startVof will not be included in the list.
   * @return Vofs in radius
   */
  static Vector<VolIndex>
  getVofsInRadius(const VolIndex&        a_startVof,
                  const EBISBox&         a_ebisbox,
                  const int              a_radius,
                  const DenseIntVectSet& a_validCells,
                  const Connectivity     a_connectivity,
                  const bool             a_addStartVof);

  /**
   * @brief Fetch vofs in a quadrant, using a specified radius and with a specified connectivitiy between them.
   * @param[in] a_startVof     Starting Vof.
   * @param[in] a_ebisbox      EBISBox.
   * @param[in] a_normal       Normal vector
   * @param[in] a_radius       Radius.
   * @param[in] a_connectivity Connectivity specification (e.g. monotone path or simplyconnected)
   * @param[in] a_addStartVof  If false, a_startVof will not be included in the list.
   * @note This routine will fetch in a "symmetric" quadrant if the normal vector aligns with a coordinate direction
   * @return Vofs in quadrant
   */
  static Vector<VolIndex>
  getVofsInQuadrant(const VolIndex&    a_startVof,
                    const EBISBox&     a_ebisbox,
                    const RealVect&    a_normal,
                    const int          a_radius,
                    const Connectivity a_connectivity,
                    const bool         a_addStartVof);

  /**
   * @brief Extract the subset of vofs that lie within the quadrant/octant defined by the normal vector
   * @param[in] a_vofs Trial Vofs.
   * @param[in] a_startVof Starting vof.
   * @param[in] a_ebisbox EBISBox holding the geometric information.
   * @param[in] a_normal Normal vector.
   * @param[in] a_radius Radius of Vof neighborhood.
   * @return Returns the subset of vofs that lie within the quadrant/octant.
   */
  static Vector<VolIndex>
  restrictVofsToQuadrant(const Vector<VolIndex>& a_vofs,
                         const VolIndex&         a_startVof,
                         const EBISBox&          a_ebisbox,
                         const RealVect&         a_normal,
                         const int               a_radius);

  /**
   * @brief Get all vofs lying in a semi-circle (sphere in 3D) from the indicated position in the starting vof
   * @param[in] a_startVof     Starting Vof.
   * @param[in] a_ebisbox      EBISBox.
   * @param[in] a_normal       Normal vector
   * @param[in] a_radius       Radius.
   * @param[in] a_deltaThresh  Threshold for the normal distance betw
   * @param[in] a_connectivity Connectivity specification (e.g. monotone path or simplyconnected)
   * @param[in] a_vofLocation  Point inside the current vof
   * @param[in] a_cellLocation Other cell locations to test with.
   * @param[in] a_addStartVof  If false, a_startVof will not be included in the list.
   * @note This routine will fetch vofs within a_radius and discard the ones that fall on the ''wrong'' side of the
   * semi-circle (sphere in 3D). Here, the wrong side is defined as all points (x-x0)*n < a_deltaThresh
   * @return Vofs in semi circle
   */
  static Vector<VolIndex>
  getVofsInSemiCircle(const VolIndex&      a_startVof,
                      const EBISBox&       a_ebisbox,
                      const RealVect&      a_normal,
                      const int            a_radius,
                      const Real           a_deltaThresh,
                      const Connectivity   a_connectivity,
                      const Location::Cell a_vofLocation,
                      const Location::Cell a_cellLocation,
                      const bool           a_addStartVof);

  /**
   * @brief Extract the subset of VoFs that lie in a semi-circle ahead of the input vof (a_startVof)
   * @note This routine will discard vofs that fall on the ''wrong'' side of the semi-circle (sphere in 3D). Here, the
   * wrong side is defined as all points (x-x0)*n < a_deltaThresh
   * @param[in] a_vofs         Trial Vofs.
   * @param[in] a_startVof     Starting Vof.
   * @param[in] a_ebisbox      EBISBox.
   * @param[in] a_normal       Normal vector
   * @param[in] a_deltaThresh  Threshold for the normal distance betw
   * @param[in] a_vofLocation  Point inside the current vof
   * @param[in] a_cellLocation Other cell locations to test with.
   * @return Returns the subset of vofs that lie within the semi-circle/sphere ahead of the starting vof.
   */
  static Vector<VolIndex>
  restrictVofsToSemiCircle(const Vector<VolIndex>& a_vofs,
                           const VolIndex&         a_startVof,
                           const EBISBox&          a_ebisbox,
                           const RealVect&         a_normal,
                           const Real              a_deltaThresh,
                           const Location::Cell    a_vofLocation,
                           const Location::Cell    a_cellLocation);

  /**
   * @brief Modifies a_vofs so it doesn't include cells in a_excludeIVS
   * @param[in] a_vofs List of vofs to be modified
   * @param[in] a_excludeBox Box to be excluded
   */
  static void
  excludeCells(Vector<VolIndex>& a_vofs, const Box& a_excludeBox);

  /**
   * @brief Modifies a_vofs so it doesn't include cells in a_excludeIVS
   * @param[in] a_vofs List of vofs to be modified
   * @param[in] a_excludeIVS List of cells to be excluded.
   */
  static void
  excludeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_excludeIVS);

  /**
   * @brief Modifies a_vofs so it only includes cells in a box
   * @param[in] a_vofs       List of vofs to be modified
   * @param[in] a_includeBox List of cells to be excluded.
   */
  static void
  includeCells(Vector<VolIndex>& a_vofs, const Box& a_includeBox);

  /**
   * @brief Modifies a_vofs so it only includes cells in a box
   * @param[in] a_vofs       List of vofs to be modified
   * @param[in] a_includeIVS List of cells to be included
   */
  static void
  includeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_includeIVS);

  /**
   * @brief Modifies a_vofs so it only includes cells in a box
   * @param[in] a_vofs       List of vofs to be modified
   * @param[in] a_includeIVS List of cells to be included
   */
  static void
  includeCells(Vector<VolIndex>& a_vofs, const DenseIntVectSet& a_includeIVS);

  /**
   * @brief Discard duplicates
   * @param[in,out] a_vofs List of Vofs
   */
  static void
  onlyUnique(Vector<VolIndex>& a_vofs);

  /**
   * @brief Get the VoFs on a finer layout that coarsen to a coarse VoF.
   * @details Finds the VoFs among the cells of refine(coarse cell, a_refRat) whose coarsening by a_refRat is a_coarVoF.
   * They are listed in the order EBISLayout::refine lists them for a graph that Chombo coarsened: at each factor of
   * two, the VoFs one level finer are gathered in BoxIterator order of their cells, and each connected set is listed in
   * the order it grows from its first VoF, as in EBGraph::getVoFSets. VoFs on intermediate levels are identified
   * through coarsening, and two of them are connected when some pair of the finer VoFs under them shares a face. Only
   * coarsening links are read, so the coarse graph needs no links to the finer level.
   * @param[in] a_coarVoF   Coarse VoF.
   * @param[in] a_ebislFine Finer layout. Shares DataIndex with the layout of a_coarVoF, holds every cell of
   * refine(coarse cell, a_refRat), and has a maximum coarsening ratio of at least a_refRat.
   * @param[in] a_din       Grid index on both layouts.
   * @param[in] a_refRat    Refinement ratio between the two layouts. A power of two.
   * @return The finer VoFs that coarsen to a_coarVoF.
   */
  static Vector<VolIndex>
  getFinerVoFs(const VolIndex& a_coarVoF, const EBISLayout& a_ebislFine, const DataIndex& a_din, const int a_refRat);

  /**
   * @brief Get the faces on a finer layout that coarsen to a coarse face.
   * @details For an interior face, lists the faces between the finer VoFs of its low and high VoFs, pairing them in the
   * order of getFinerVoFs at each factor of two as EBISLayout::refine does. For a face on the domain boundary, lists
   * the boundary faces of the finer VoFs of its interior VoF that lie against the boundary.
   * @param[in] a_coarFace  Coarse face.
   * @param[in] a_ebislFine Finer layout. Shares DataIndex with the layout of a_coarFace, holds every cell of the
   * refined cells on both sides of the face, and has a maximum coarsening ratio of at least a_refRat.
   * @param[in] a_din       Grid index on both layouts.
   * @param[in] a_refRat    Refinement ratio between the two layouts. A power of two.
   * @return The finer faces that coarsen to a_coarFace.
   */
  static Vector<FaceIndex>
  getFinerFaces(const FaceIndex& a_coarFace, const EBISLayout& a_ebislFine, const DataIndex& a_din, const int a_refRat);

  /**
   * @brief Get the VoFs on a finer layout that coarsen to a coarse VoF and lie in the finer EBISBox, in no particular
   * order.
   * @details For callers that need only the set of finer VoFs, and may ask about coarse cells whose finer cells reach
   * past the finer EBISBox. Finer cells outside the EBISBox are left out. The VoFs are listed in BoxIterator order of
   * their cells, then by cell index. This costs less than getFinerVoFs, which also orders them.
   * @param[in] a_coarVoF   Coarse VoF.
   * @param[in] a_ebislFine Finer layout. Shares DataIndex with the layout of a_coarVoF, and has a maximum coarsening
   * ratio of at least a_refRat.
   * @param[in] a_din       Grid index on both layouts.
   * @param[in] a_refRat    Refinement ratio between the two layouts. A power of two.
   * @return The finer VoFs in the finer EBISBox that coarsen to a_coarVoF.
   */
  static Vector<VolIndex>
  getFinerVoFsInRegion(const VolIndex&   a_coarVoF,
                       const EBISLayout& a_ebislFine,
                       const DataIndex&  a_din,
                       const int         a_refRat);

protected:
  /**
   * @brief Get the VoFs on a finer layout that coarsen to a coarse VoF, in gathering order.
   * @details Gathering order is BoxIterator order of the cells of refine(coarse cell, a_refRat), then cell index.
   * @param[in] a_coarVoF   Coarse VoF.
   * @param[in] a_ebislFine Finer layout.
   * @param[in] a_din       Grid index on both layouts.
   * @param[in] a_refRat    Refinement ratio between the two layouts.
   * @return The finer VoFs that coarsen to a_coarVoF.
   */
  static Vector<VolIndex>
  gatherFinerVoFs(const VolIndex& a_coarVoF, const EBISLayout& a_ebislFine, const DataIndex& a_din, const int a_refRat);

  /**
   * @brief Group the finest VoFs under one VoF by their ancestor one factor of two finer than it, and order the groups.
   * @details The groups are ordered as EBGraph::getVoFSets lists the VoFs of one connected set: gathered in BoxIterator
   * order of the ancestors' cells, then cell index, and listed in the order the set grows from its first member. Two
   * groups are connected when some pair of their finest VoFs shares a face.
   * @param[in] a_finest    The finest VoFs under one VoF.
   * @param[in] a_ebislFine Finest layout.
   * @param[in] a_din       Grid index.
   * @param[in] a_ratio     Ratio between the level of the VoF the finest VoFs lie under and the finest level. At least
   * two.
   * @return The finest VoFs of each group, groups in order.
   */
  static Vector<Vector<VolIndex>>
  orderFinerGroups(const Vector<VolIndex>& a_finest,
                   const EBISLayout&       a_ebislFine,
                   const DataIndex&        a_din,
                   const int               a_ratio);

  /**
   * @brief Order the finest VoFs under one VoF as EBISLayout::refine lists them.
   * @param[in] a_finest    The finest VoFs under one VoF.
   * @param[in] a_ebislFine Finest layout.
   * @param[in] a_din       Grid index.
   * @param[in] a_ratio     Ratio between the level of the VoF the finest VoFs lie under and the finest level.
   * @return a_finest in order.
   */
  static Vector<VolIndex>
  orderFinerVoFs(const Vector<VolIndex>& a_finest,
                 const EBISLayout&       a_ebislFine,
                 const DataIndex&        a_din,
                 const int               a_ratio);

  /**
   * @brief Append the finest faces between the finest VoFs under two VoFs that share an interior face.
   * @param[in,out] a_faces     Faces, appended to.
   * @param[in]     a_finestLo  The finest VoFs under the low VoF.
   * @param[in]     a_finestHi  The finest VoFs under the high VoF.
   * @param[in]     a_ebislFine Finest layout.
   * @param[in]     a_din       Grid index.
   * @param[in]     a_ratio     Ratio between the level of the two VoFs and the finest level.
   * @param[in]     a_dir       Face direction.
   */
  static void
  appendFinerFaces(Vector<FaceIndex>&      a_faces,
                   const Vector<VolIndex>& a_finestLo,
                   const Vector<VolIndex>& a_finestHi,
                   const EBISLayout&       a_ebislFine,
                   const DataIndex&        a_din,
                   const int               a_ratio,
                   const int               a_dir);

  /**
   * @brief Append the finest boundary faces under a VoF's face on the domain boundary.
   * @param[in,out] a_faces     Faces, appended to.
   * @param[in]     a_finest    The finest VoFs under the VoF.
   * @param[in]     a_ebislFine Finest layout.
   * @param[in]     a_din       Grid index.
   * @param[in]     a_ratio     Ratio between the level of the VoF and the finest level.
   * @param[in]     a_dir       Face direction.
   * @param[in]     a_side      Side of the VoF on which the domain boundary lies.
   */
  static void
  appendFinerBoundaryFaces(Vector<FaceIndex>&      a_faces,
                           const Vector<VolIndex>& a_finest,
                           const EBISLayout&       a_ebislFine,
                           const DataIndex&        a_din,
                           const int               a_ratio,
                           const int               a_dir,
                           const Side::LoHiSide    a_side);

  /**
   * @brief Check if a quadrant is well-defined for a direction in space.
   * @param[in] a_normal The direction in space which defines the quadrant.
   * @return True if the quadrant is well defined and false otherwise.
   */
  static bool
  isQuadrantWellDefined(const RealVect& a_normal);

  /**
   * @brief Returns the coordinate direction
   * @return Cardinal direction
   * @param[in] a_normal Normal
   */
  static std::pair<int, Side::LoHiSide>
  getCardinalDirection(const RealVect& a_normal);

  /**
   * @brief Get all VoFs that can be found with a monotone path with specified radius
   * @param[in] a_startVof    Starting Vof.
   * @param[in] a_ebisbox     EBISBox.
   * @param[in] a_radius      Radius.
   * @param[in] a_addStartVof If false, a_startVof will not be included in the list.
   * @return Returns list of VoFs that can be found with a monotone path of radius <= a_radius.
   */
  static Vector<VolIndex>
  getVofsInMonotonePath(const VolIndex& a_startVof,
                        const EBISBox&  a_ebisbox,
                        const int       a_radius,
                        const bool      a_addStartVof);

  /**
   * @brief Get all connected Vofs in radius which can be reached without crossing a covered face or reaching into a
   * domain boundary.
   * @param[in] a_startVof    Starting Vof
   * @param[in] a_ebisbox     EBISBox
   * @param[in] a_radius      Search radius
   * @param[in] a_addStartVof If false, a_startVof will not be in the returned list.
   * @return Returns all Vofs within radius a_radius from a_startVof
   */
  static Vector<VolIndex>
  getConnectedVofsInRadius(const VolIndex& a_startVof,
                           const EBISBox&  a_ebisbox,
                           const int       a_radius,
                           const bool      a_addStartVof);

  /**
   * @brief Get all Vofs in a radius from a starting Vof
   * @param[in] a_startVof   Starting Vof
   * @param[in] a_ebisbox    EBISBox
   * @param[in] a_radius     Search radius
   * @param[in] a_addStartVof If false, a_startVof will not be included in the list.
   * @return All vofs in radius
   */
  static Vector<VolIndex>
  getAllVofsInRadius(const VolIndex& a_startVof,
                     const EBISBox&  a_ebisbox,
                     const int       a_radius,
                     const bool      a_addStartVof);

  /**
   * @brief Get all Vofs that are connected to a specific starting Vof. The connection does not have to be direct.
   * @param[in] a_startVof Starting Vof
   * @param[in] a_allVofs  All Vofs. This can, but does not have to, include a_startVof.
   * @param[in] a_ebisbox  EBISBox.
   * @return List of Vofs that are connected to a_startVof. This list does *not* include a_startVof
   */
  static Vector<VolIndex>
  connectedVofsOnly(const VolIndex& a_startVof, const Vector<VolIndex>& a_allVofs, const EBISBox& a_ebisbox);

  /**
   * @brief Get quadrant defined by direction in space.
   * @param[in] a_vof      Starting vof
   * @param[in] a_ebisbox  EBISBox
   * @param[in] a_radius   Quadrant radius
   * @return Quadrant as a box, with one of the corners on a_vof.
   * @param[in] a_normal Normal
   */
  static Box
  getQuadrant(const RealVect& a_normal, const VolIndex& a_vof, const EBISBox& a_ebisbox, const Real a_radius);

  /**
   * @brief Get "symmetric" quadrant defined by direction in space.
   * @param[in] a_cardinal The direction in space which defines the quadrant.
   * @param[in] a_vof      Starting vof
   * @param[in] a_radius   Quadrant radius
   * @return Symmetric quadrant with a_vof on one of the sides.
   * @param[in] a_ebisbox Description
   */
  static Box
  getSymmetricQuadrant(const std::pair<int, Side::LoHiSide>& a_cardinal,
                       const VolIndex&                       a_vof,
                       const EBISBox&                        a_ebisbox,
                       const Real                            a_radius);

  /**
   * @brief Get all Vofs in a monotone path. Does not move more than a_radius in any direction.
   * @param[in,out] a_vofList    List of vofs that were within the monotone path
   * @param[in]    a_startVof   Starting vof.
   * @param[in]    a_ebisbox    EBISBox.
   * @param[in]    a_radius     Maximum radius.
   * @param[in]    a_timesMoved Number of times moved in any direction
   * @param[in]    a_pathSign   Path sign.
   */
  static void
  getVofsInMonotonePath(Vector<VolIndex>& a_vofList,
                        const VolIndex&   a_startVof,
                        const EBISBox&    a_ebisbox,
                        const int         a_radius,
                        const IntVect&    a_timesMoved,
                        const IntVect&    a_pathSign);

  /**
   * @brief Get all Vofs in a monotone path. Does not move more than a_radius in any direction.
   * @param[in,out] a_vofSet     List of vofs that were within the monotone path
   * @param[in]    a_startVof   Starting vof.
   * @param[in]    a_ebisbox    EBISBox.
   * @param[in]    a_radius     Maximum radius.
   * @param[in]    a_timesMoved Number of times moved in any direction
   * @param[in]    a_pathSign   Path sign.
   * @param[in] a_validBox Description
   */
  static void
  getVofsInMonotonePath(std::set<VolIndex>& a_vofSet,
                        const VolIndex&     a_startVof,
                        const EBISBox&      a_ebisbox,
                        const Box&          a_validBox,
                        const int&          a_radius,
                        const IntVect&      a_timesMoved,
                        const IntVect&      a_pathSign);
};

#include <CD_NamespaceFooter.H>

#endif
