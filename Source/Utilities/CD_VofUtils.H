/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_VofUtils.H
  @brief  Various functions for getting Vofs near cut-cells
  @author Robert Marskar
*/

#ifndef CD_VofUtils_H
#define CD_VofUtils_H

// Chombo includes
#include <Stencils.H>
#include <EBISBox.H>
#include <IntVectSet.H>
#include <RealVect.H>
#include <ProblemDomain.H>
#include <VolIndex.H>
#include <LoHiSide.H>
#include <DenseIntVectSet.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Static class which contains some routines for fetching VoFs using various algorithms. Very useful in combination with least squares. 
*/
class VofUtils
{
public:
  /*!
    @brief Class for specifying connectivity of Vofs. 
    @details MontonePath:     Vofs must be reachable with monotone path. 
    @details SimplyConnected: All vofs must be connected
    @details All:             All vofs allowed
  */
  enum class Connectivity
  {
    MonotonePath,
    SimplyConnected,
    All
  };

  /*!
    @brief Neighborhood type
    @param[in] Quadrant.   Only fetch vofs in a specified quadrant. 
    @param[in] Radius.     Fetch all vofs in a specified radius. 
    @param[in] SemiCircle. Fetch all vofs in a half-circle with specified radius. 
  */
  enum class Neighborhood
  {
    Quadrant,
    Radius,
    SemiCircle
  };

  /*!
    @brief Disallowed
  */
  VofUtils() = delete;

  /*!
    @brief Fetch vofs in a specified radius and with a specified connectivitiy between them.
    @param[in] a_startVof     Starting Vof. 
    @param[in] a_ebisbox      EBISBox.
    @param[in] a_radius       Radius.
    @param[in] a_connectivity Connectivity specification (e.g. monotone path or simplyconnected)
    @param[in] a_addStartVof  If false, a_startVof will not be included in the list. 
  */
  static Vector<VolIndex>
  getVofsInRadius(const VolIndex&    a_startVof,
                  const EBISBox&     a_ebisbox,
                  const int          a_radius,
                  const Connectivity a_connectivity,
                  const bool         a_addStartVof);

  /*!
    @brief Fetch vofs in a specified radius and with a specified connectivitiy between them.
    @param[in] a_startVof     Starting Vof. 
    @param[in] a_ebisbox      EBISBox.
    @param[in] a_radius       Radius.
    @param[in] a_validCells   Valid cells to include
    @param[in] a_connectivity Connectivity specification (e.g. monotone path or simplyconnected)
    @param[in] a_addStartVof  If false, a_startVof will not be included in the list. 
  */
  static Vector<VolIndex>
  getVofsInRadius(const VolIndex&        a_startVof,
                  const EBISBox&         a_ebisbox,
                  const int              a_radius,
                  const DenseIntVectSet& a_validCells,
                  const Connectivity     a_connectivity,
                  const bool             a_addStartVof);

  /*!
    @brief Fetch vofs in a quadrant, using a specified radius and with a specified connectivitiy between them.
    @param[in] a_startVof     Starting Vof. 
    @param[in] a_ebisbox      EBISBox.
    @param[in] a_normal       Normal vector
    @param[in] a_radius       Radius.
    @param[in] a_connectivity Connectivity specification (e.g. monotone path or simplyconnected)
    @param[in] a_addStartVof  If false, a_startVof will not be included in the list. 
    @note This routine will fetch in a "symmetric" quadrant if the normal vector aligns with a coordinate direction
  */
  static Vector<VolIndex>
  getVofsInQuadrant(const VolIndex&    a_startVof,
                    const EBISBox&     a_ebisbox,
                    const RealVect&    a_normal,
                    const int          a_radius,
                    const Connectivity a_connectivity,
                    const bool         a_addStartVof);

  /*!
    @brief Extract the subset of vofs that lie within the quadrant/octant defined by the normal vector
    @param[in] a_vofs Trial Vofs.
    @param[in] a_startVof Starting vof.
    @param[in] a_ebisbox EBISBox holding the geometric information.
    @param[in] a_normal Normal vector.
    @param[in] a_radius Radius of Vof neighborhood.
    @return Returns the subset of vofs that lie within the quadrant/octant.
  */
  static Vector<VolIndex>
  restrictVofsToQuadrant(const Vector<VolIndex>& a_vofs,
                         const VolIndex&         a_startVof,
                         const EBISBox&          a_ebisbox,
                         const RealVect&         a_normal,
                         const int               a_radius);

  /*!
    @brief Get all vofs lying in a semi-circle (sphere in 3D) from the indicated position in the starting vof
    @param[in] a_startVof     Starting Vof. 
    @param[in] a_ebisbox      EBISBox.
    @param[in] a_normal       Normal vector
    @param[in] a_radius       Radius.
    @param[in] a_deltaThresh  Threshold for the normal distance betw
    @param[in] a_connectivity Connectivity specification (e.g. monotone path or simplyconnected)
    @param[in] a_vofLocation  Point inside the current vof
    @param[in] a_cellLocation Other cell locations to test with. 
    @param[in] a_addStartVof  If false, a_startVof will not be included in the list. 
    @note This routine will fetch vofs within a_radius and discard the ones that fall on the ''wrong'' side of the semi-circle (sphere in 3D). Here, the
    wrong side is defined as all points (x-x0)*n < a_deltaThresh
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

  /*!
    @brief Extract the subset of VoFs that lie in a semi-circle ahead of the input vof (a_startVof)
    @note This routine will discard vofs that fall on the ''wrong'' side of the semi-circle (sphere in 3D). Here, the
    wrong side is defined as all points (x-x0)*n < a_deltaThresh    
    @param[in] a_vofs         Trial Vofs.
    @param[in] a_startVof     Starting Vof. 
    @param[in] a_ebisbox      EBISBox.
    @param[in] a_normal       Normal vector
    @param[in] a_deltaThresh  Threshold for the normal distance betw
    @param[in] a_vofLocation  Point inside the current vof
    @param[in] a_cellLocation Other cell locations to test with.
    @return Returns the subset of vofs that lie within the semi-circle/sphere ahead of the starting vof.
  */
  static Vector<VolIndex>
  restrictVofsToSemiCircle(const Vector<VolIndex>& a_vofs,
                           const VolIndex&         a_startVof,
                           const EBISBox&          a_ebisbox,
                           const RealVect&         a_normal,
                           const Real              a_deltaThresh,
                           const Location::Cell    a_vofLocation,
                           const Location::Cell    a_cellLocation);

  /*!
    @brief Modifies a_vofs so it doesn't include cells in a_excludeIVS
    @param[in] a_vofs List of vofs to be modified
    @param[in] a_excludeBox Box to be excluded
    @return Modifed list of vofs which do not include cell in a box
  */
  static void
  excludeCells(Vector<VolIndex>& a_vofs, const Box& a_excludeBox);

  /*!
    @brief Modifies a_vofs so it doesn't include cells in a_excludeIVS
    @param[in] a_vofs List of vofs to be modified
    @param[in] a_excludeIVS List of cells to be excluded. 
    @return Modifed list of vofs which do not include a_excludeIVS cells. 
  */
  static void
  excludeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_excludeIVS);

  /*!
    @brief Modifies a_vofs so it only includes cells in a box
    @param[in] a_vofs       List of vofs to be modified
    @param[in] a_includeBox List of cells to be excluded. 
    @return Modifed list of vofs which are guaranteed to lie in a_box
  */
  static void
  includeCells(Vector<VolIndex>& a_vofs, const Box& a_includeBox);

  /*!
    @brief Modifies a_vofs so it only includes cells in a box
    @param[in] a_vofs       List of vofs to be modified
    @param[in] a_includeIVS List of cells to be included
    @return Modifed list of vofs which are guaranteed to lie in specified grid cells
  */
  static void
  includeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_includeIVS);

  /*!
    @brief Modifies a_vofs so it only includes cells in a box
    @param[in] a_vofs       List of vofs to be modified
    @param[in] a_includeIVS List of cells to be included
    @return Modifed list of vofs which are guaranteed to lie in specified grid cells
  */
  static void
  includeCells(Vector<VolIndex>& a_vofs, const DenseIntVectSet& a_includeIVS);

  /*!
    @brief Discard duplicates
    @param[inout] a_vofs List of Vofs
    @return Returns the same vofs with discarded duplicates. 
  */
  static void
  onlyUnique(Vector<VolIndex>& a_vofs);

protected:
  /*!
    @brief Check if a quadrant is well-defined for a direction in space.
    @param[in] a_normal The direction in space which defines the quadrant. 
    @return True if the quadrant is well defined and false otherwise. 
  */
  static bool
  isQuadrantWellDefined(const RealVect a_normal);

  /*!
    @brief Returns the coordinate direction
    @param[in] Non-zero normal vector
    @param[in] Returns dirction where ||a_normal[dir]|| == 1. If a_normal[dir] == -1 then dir => -dir
  */
  static std::pair<int, Side::LoHiSide>
  getCardinalDirection(const RealVect a_normal);

  /*!
    @brief Get all VoFs that can be found with a monotone path with specified radius
    @param[in] a_startVof    Starting Vof. 
    @param[in] a_ebisbox     EBISBox.
    @param[in] a_radius      Radius.
    @param[in] a_addStartVof If false, a_startVof will not be included in the list. 
    @return Returns list of VoFs that can be found with a monotone path of radius <= a_radius. 
  */
  static Vector<VolIndex>
  getVofsInMonotonePath(const VolIndex& a_startVof,
                        const EBISBox&  a_ebisbox,
                        const int       a_radius,
                        const bool      a_addStartVof);

  /*!
    @brief Get all connected Vofs in radius which can be reached without crossing a covered face or reaching into a domain boundary. 
    @param[in] a_startVof    Starting Vof
    @param[in] a_ebisbox     EBISBox
    @param[in] a_radius      Search radius
    @param[in] a_addStartVof If false, a_startVof will not be in the returned list. 
    @return Returns all Vofs within radius a_radius from a_startVof
  */
  static Vector<VolIndex>
  getConnectedVofsInRadius(const VolIndex& a_startVof,
                           const EBISBox&  a_ebisbox,
                           const int       a_radius,
                           const bool      a_addStartVof);

  /*!
    @brief Get all Vofs in a radius from a starting Vof
    @param[in] a_startVof   Starting Vof
    @param[in] a_ebisbox    EBISBox
    @param[in] a_radius     Search radius
    @param[in] a_addStartVof If false, a_startVof will not be included in the list. 
  */
  static Vector<VolIndex>
  getAllVofsInRadius(const VolIndex& a_startVof,
                     const EBISBox&  a_ebisbox,
                     const int       a_radius,
                     const bool      a_addStartVof);

  /*!
    @brief Get all Vofs that are connected to a specific starting Vof. The connection does not have to be direct. 
    @param[in] a_startVof Starting Vof
    @param[in] a_allVofs  All Vofs. This can, but does not have to, include a_startVof.
    @param[in] a_ebisbox  EBISBox. 
    @return List of Vofs that are connected to a_startVof. This list does *not* include a_startVof
  */
  static Vector<VolIndex>
  connectedVofsOnly(const VolIndex& a_startVof, const Vector<VolIndex>& a_allVofs, const EBISBox& a_ebisbox);

  /*!
    @brief Get quadrant defined by direction in space. 
    @param[in] a_cardinal The direction in space which defines the quadrant. 
    @param[in] a_vof      Starting vof
    @param[in] a_ebisbox  EBISBox
    @param[in] a_radius   Quadrant radius
    @return Quadrant as a box, with one of the corners on a_vof. 
  */
  static Box
  getQuadrant(const RealVect& a_normal, const VolIndex& a_vof, const EBISBox& a_ebisbox, const Real a_radius);

  /*!
    @brief Get "symmetric" quadrant defined by direction in space. 
    @param[in] a_cardinal The direction in space which defines the quadrant. 
    @param[in] a_vof      Starting vof
    @param[in] a_radius   Quadrant radius
    @return Symmetric quadrant with a_vof on one of the sides. 
  */
  static Box
  getSymmetricQuadrant(const std::pair<int, Side::LoHiSide>& a_cardinal,
                       const VolIndex&                       a_vof,
                       const EBISBox&                        a_ebisbox,
                       const Real                            a_radius);

  /*!
    @brief Get all Vofs in a monotone path. Does not move more than a_radius in any direction. 
    @param[inout] a_vofList    List of vofs that were within the monotone path
    @param[in]    a_startVof   Starting vof. 
    @param[in]    a_ebisbox    EBISBox.
    @param[in]    a_radius     Maximum radius. 
    @param[in]    a_timesMoved Number of times moved in any direction
    @param[in]    a_pathSign   Path sign. 
  */
  static void
  getVofsInMonotonePath(Vector<VolIndex>& a_vofList,
                        const VolIndex&   a_startVof,
                        const EBISBox&    a_ebisbox,
                        const int         a_radius,
                        const IntVect&    a_timesMoved,
                        const IntVect&    a_pathSign);

  /*!
    @brief Get all Vofs in a monotone path. Does not move more than a_radius in any direction. 
    @param[inout] a_vofList    List of vofs that were within the monotone path
    @param[in]    a_startVof   Starting vof. 
    @param[in]    a_ebisbox    EBISBox.
    @param[in]    a_radius     Maximum radius. 
    @param[in]    a_timesMoved Number of times moved in any direction
    @param[in]    a_pathSign   Path sign. 
  */
  static void
  getVofsInMonotonePath(std::set<VolIndex>& a_vofList,
                        const VolIndex&     a_startVof,
                        const EBISBox&      a_ebisbox,
                        const Box&          a_validBox,
                        const int&          a_radius,
                        const IntVect&      a_timesMoved,
                        const IntVect&      a_pathSign);
};

#include <CD_NamespaceFooter.H>

#endif
