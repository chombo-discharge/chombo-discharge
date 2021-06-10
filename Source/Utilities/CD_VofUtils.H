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

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Static class which contains some routines for fetching VoFs using various algorithms. Very useful in combination with least squares. 
*/
class VofUtils {
public:

  /*!
    @brief Disallowed
  */
  VofUtils() = delete;

  /*!
    @brief Check if a quadrant is well-defined for a direction in space.
    @param[in] a_normal The direction in space which defines the quadrant. 
    @return True if the quadrant is well defined and false otherwise. 
  */
  static bool isQuadrantWellDefined(const RealVect a_normal);

  /*!
    @brief Returns the coordinate direction
    @param[in] Non-zero normal vector
    @param[in] Returns dirction where ||a_normal[dir]|| == 1. If a_normal[dir] == -1 then dir => -dir
  */
  static std::pair<int, Side::LoHiSide> getCardinalDirection(const RealVect a_normal);

  /*!
    @brief Get quadrant defined by direction in space. 
    @param[in] a_normal The direction in space which defines the quadrant. 
    @return Quadrant as IntVect. The IntVect is ill-defined if the quadrant is not well-defined. 
  */
  static IntVect getQuadrant(const RealVect a_normal);

  /*!
    @brief Get all Vofs in a radius from a starting Vof
    @param[in] a_startVof   Starting Vof
    @param[in] a_ebisbox    EBISBox
    @param[in] a_radius     Search radius
  */
  static Vector<VolIndex> getAllVofsInRadius(const VolIndex& a_startVof, const EBISBox& a_ebisbox, const int a_radius, const bool a_addStartVof);

  /*!
    @brief Get all Vofs within a box.
    @param[in] a_box     Box where you want to get Vofs. 
    @param[in] a_ebisbox EBISBox
    @details Avoids fetching stuff outside the domain. 
  */
  static Vector<VolIndex> getAllVofsInBox(const Box& a_box, const EBISBox& a_ebisbox);

  /*!
    @brief Modifies a_vofs so it doesn't include cells in a_excludeIVS
    @param[in] a_vofs List of vofs to be modified
    @param[in] a_excludeBox Box to be excluded
    @return Modifed list of vofs which do not include cell in a box
  */
  static void excludeCells(Vector<VolIndex>& a_vofs, const Box& a_excludeBox);

  /*!
    @brief Modifies a_vofs so it doesn't include cells in a_excludeIVS
    @param[in] a_vofs List of vofs to be modified
    @param[in] a_excludeIVS List of cells to be excluded. 
    @return Modifed list of vofs which do not include a_excludeIVS cells. 
  */
  static void excludeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_excludeIVS);

  /*!
    @brief Modifies a_vofs so it only includes cells in a box
    @param[in] a_vofs       List of vofs to be modified
    @param[in] a_includeBox List of cells to be excluded. 
    @return Modifed list of vofs which are guaranteed to lie in a_box
  */
  static void includeCells(Vector<VolIndex>& a_vofs, const Box& a_includeBox);

  /*!
    @brief Modifies a_vofs so it only includes cells in a box
    @param[in] a_vofs       List of vofs to be modified
    @param[in] a_includeIVS List of cells to be included
    @return Modifed list of vofs which are guaranteed to lie in specified grid cells
  */
  static void includeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_includeIVS);

  /*!
    @brief Get all Vofs that are connected to a specific starting Vof. The connection does not have to be direct. 
    @param[in] a_startVof Starting Vof
    @param[in] a_allVofs  All Vofs. This can, but does not have to, include a_startVof.
    @param[in] a_ebisbox  EBISBox. 
    @return List of Vofs that are connected to a_startVof. This list does *not* include a_startVof
  */
  static Vector<VolIndex> getConnectedVofs(const VolIndex& a_startVof, const Vector<VolIndex>& a_allVofs, const EBISBox& a_ebisbox);

  /*!
    @brief Get all connected Vofs in radius which can be reached without crossing a covered face or reaching into a domain boundary. 
    @param[in] a_startVof    Starting Vof
    @param[in] a_ebisbox     EBISBox
    @param[in] a_radius      Search radius
    @param[in] a_excludeIVS Do not add Vofs with these grid indices. 
    @return Returns all Vofs within radius a_radius from a_startVof
  */
  static Vector<VolIndex> getAllConnectedVofsInRadius(const VolIndex&   a_startVof,
						      const EBISBox&    a_ebisbox,
						      const int         a_radius,
						      const IntVectSet& a_excludeIVS);
  
  /*!
    @brief Get all connected Vofs in a box which can be reached without crossing a covered face or reaching into a domain boundary. 
    @param[in] a_startVof    Starting Vof
    @param[in] a_box         Box
    @param[in] a_ebisbox     EBISBox
    @param[in] a_excludeIVS Do not add Vofs with these grid indices. 
    @return Returns all Vofs within radius a_radius from a_startVof
  */
  static Vector<VolIndex> getAllConnectedVofsInBox(const VolIndex&   a_startVof,
						   const EBISBox&    a_ebisbox,
						   const Box&        a_box,
						   const IntVectSet& a_excludeIVS);

  /*!
    @brief Get all Vofs in a quadrant. The quadrant is defined by the normal vector. 
    @param[in] a_startVof     Starting Vof
    @param[in] a_ebisBox      EBISBox
    @param[in] a_normal       Normal vector which defines the quadrant.
    @param[in] a_radius       Radius of quadrant
    @param[in] a_addStartVof  If true, adds a_startVof to the list. 
    @note This creates a quadrant and calls the other version. 
    @return Vector of Vofs in quadrant defined by the normal vector. If a_normal[dir]==0 there is no unique quadrant and the routine returns an empty vector. 
  */
  static Vector<VolIndex> getAllVofsInQuadrant(const VolIndex& a_startVof, const EBISBox& a_ebisbox, const RealVect& a_normal, const int a_radius, const bool a_addStartVof);

  /*!
    @brief Get all Vofs in a quadrant. The quadrant is defined by the normal vector. 
    @param[in] a_startVof     Starting Vof
    @param[in] a_ebisBox      EBISBox
    @param[in] a_normal       Normal vector which defines the quadrant.
    @param[in] a_quadrant     Quadrant
    @param[in] a_addStartVof  If true, adds a_startVof to the list. 
    @return Vector of Vofs in quadrant defined by the normal vector. If a_normal[dir]==0 there is no unique quadrant and the routine returns an empty vector. 
  */
  static Vector<VolIndex> getAllVofsInQuadrant(const VolIndex& a_startVof, const EBISBox& a_ebisbox, const IntVect a_quadrant, const int a_radius, const bool a_addStartVof);

  /*!
    @brief Get all Vofs within a radius, but ignore Vofs that are ''opposite'' to a specified normal direction in space. 
    @param[in] a_startVof    Starting Vof. 
    @param[in] a_ebisbox     EBISBox.
    @param[in] a_cardinal    Specifies cardinal direction in space.
    @param[in] a_side        Specifies direction in space.
    @param[in] a_radius      Radius.
    @param[in] a_addStartVof If true, add a_startVof to the returned vector. 
  */
  static Vector<VolIndex> getAllVofsSymmetric(const VolIndex&      a_startVof,
					      const EBISBox&       a_ebisbox,
					      const int            a_cardinal,
					      const Side::LoHiSide a_side,
					      const int            a_radius,
					      const bool           a_addStartVof);

  /*!
    @brief Get all VoFs that can be found with a monotone path with specified radius
    @param[in] a_startVof    Starting Vof. 
    @param[in] a_ebisbox     EBISBox.
    @param[in] a_radius      Radius.
    @param[in] a_addStartVof If false, a_startVof will not be included in the list. 
    @return Returns list of VoFs that can be found with a monotone path of radius <= a_radius. 
  */
  static Vector<VolIndex> getAllVofsInMonotonePath(VolIndex& a_startVof, const EBISBox& a_ebisbox, const int a_radius, const bool a_addStartVof);

protected:

  /*!
    @brief Get all Vofs in a monotone path. Does not move more than a_radius in any direction. 
    @param[inout] a_vofList    List of vofs that were within the monotone path
    @param[in]    a_startVof   Starting vof. 
    @param[in]    a_ebisbox    EBISBox.
    @param[in]    a_radius     Maximum radius. 
    @param[in]    a_timesMoved Number of times moved in any direction
    @param[in]    a_pathSign   Path sign. 
  */
  static void getAllVofsInMonotonePath(Vector<VolIndex>& a_vofList,
				       const VolIndex&   a_startVof,
				       const EBISBox&    a_ebisbox,
				       const int         a_radius,
				       const IntVect&    a_timesMoved,
				       const IntVect&    a_pathSign);
};

#include <CD_NamespaceFooter.H>

#endif
  
