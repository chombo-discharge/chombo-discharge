/*!
  @file   VoFUtils.H
  @brief  Various functions for getting VoFs near cut-cells
  @author Robert Marskar
  @date   May 2021
*/

#ifndef _VoFUtils_H_
#define _VoFUtils_H_

#include "Stencils.H"
#include "EBISBox.H"
#include "IntVectSet.H"
#include "RealVect.H"
#include "ProblemDomain.H"
#include "VolIndex.H"

namespace VoFUtils {

  /*!
    @brief Check if a quadrant is well-defined for a direction in space.
    @param[in] a_normal The direction in space which defines the quadrant. 
    @return True if the quadrant is well defined and false otherwise. 
  */
  bool isQuadrantWellDefined(const RealVect a_normal);

  /*!
    @brief Returns the coordinate direction
    @param[in] Non-zero normal vector
    @param[in] Returns dirction where ||a_normal[dir]|| == 1. If a_normal[dir] == -1 then dir => -dir
  */
  std::pair<int, Side::LoHiSide> getCardinalDirection(const RealVect a_normal);

  /*!
    @brief Get quadrant defined by direction in space. 
    @param[in] a_normal The direction in space which defines the quadrant. 
    @return Quadrant as IntVect. The IntVect is ill-defined if the quadrant is not well-defined. 
  */
  IntVect getQuadrant(const RealVect a_normal);

  /*!
    @brief Get all VoFs in a radius from a starting VoF
    @param[in] a_startVoF   Starting VoF
    @param[in] a_ebisbox    EBISBox
    @param[in] a_radius     Search radius
  */
  Vector<VolIndex> getAllVoFsInRadius(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const int a_radius, const bool a_addStartVoF);

    /*!
    @brief Get all VoFs within a box.
    @param[in] a_box     Box where you want to get VoFs. 
    @param[in] a_ebisbox EBISBox
    @details Avoids fetching stuff outside the domain. 
  */
  Vector<VolIndex> getAllVoFsInBox(const Box& a_box, const EBISBox& a_ebisbox);

    /*!
    @brief Modifies a_vofs so it doesn't include cells in a_excludeIVS
    @param[in] a_vofs List of vofs to be modified
    @param[in] a_excludeBox Box to be excluded
    @return Modifed list of vofs which do not include cell in a box
  */
  void excludeCells(Vector<VolIndex>& a_vofs, const Box& a_excludeBox);

  /*!
    @brief Modifies a_vofs so it doesn't include cells in a_excludeIVS
    @param[in] a_vofs List of vofs to be modified
    @param[in] a_excludeIVS List of cells to be excluded. 
    @return Modifed list of vofs which do not include a_excludeIVS cells. 
  */
  void excludeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_excludeIVS);

  /*!
    @brief Modifies a_vofs so it only includes cells in a box
    @param[in] a_vofs       List of vofs to be modified
    @param[in] a_includeBox List of cells to be excluded. 
    @return Modifed list of vofs which are guaranteed to lie in a_box
  */
  void includeCells(Vector<VolIndex>& a_vofs, const Box& a_includeBox);

    /*!
    @brief Modifies a_vofs so it only includes cells in a box
    @param[in] a_vofs       List of vofs to be modified
    @param[in] a_includeIVS List of cells to be included
    @return Modifed list of vofs which are guaranteed to lie in specified grid cells
  */
  void includeCells(Vector<VolIndex>& a_vofs, const IntVectSet& a_includeIVS);

  /*!
    @brief Get all VoFs that are connected to a specific starting VoF. The connection does not have to be direct. 
    @param[in] a_startVoF Starting VoF
    @param[in] a_allVoFs  All VoFs. This can, but does not have to, include a_startVoF.
    @param[in] a_ebisbox  EBISBox. 
    @return List of VoFs that are connected to a_startVoF. This list does *not* include a_startVoF
  */
  Vector<VolIndex> getConnectedVoFs(const VolIndex& a_startVoF, const Vector<VolIndex>& a_allVoFs, const EBISBox& a_ebisbox);

  /*!
    @brief Get all connected VoFs in radius which can be reached without crossing a covered face or reaching into a domain boundary. 
    @param[in] a_startVoF    Starting VoF
    @param[in] a_ebisbox     EBISBox
    @param[in] a_radius      Search radius
    @param[in] a_excludeIVS Do not add VoFs with these grid indices. 
    @return Returns all VoFs within radius a_radius from a_startVoF
  */
  Vector<VolIndex> getAllConnectedVoFsInRadius(const VolIndex&   a_startVoF,
					       const EBISBox&    a_ebisbox,
					       const int         a_radius,
					       const IntVectSet& a_excludeIVS);
  
  /*!
    @brief Get all connected VoFs in a box which can be reached without crossing a covered face or reaching into a domain boundary. 
    @param[in] a_startVoF    Starting VoF
    @param[in] a_box         Box
    @param[in] a_ebisbox     EBISBox
    @param[in] a_excludeIVS Do not add VoFs with these grid indices. 
    @return Returns all VoFs within radius a_radius from a_startVoF
  */
  Vector<VolIndex> getAllConnectedVoFsInBox(const VolIndex&   a_startVoF,
					    const EBISBox&    a_ebisbox,
					    const Box&        a_box,
					    const IntVectSet& a_excludeIVS);

  /*!
    @brief Get all VoFs in a quadrant. The quadrant is defined by the normal vector. 
    @param[in] a_startVoF     Starting VoF
    @param[in] a_ebisBox      EBISBox
    @param[in] a_normal       Normal vector which defines the quadrant.
    @param[in] a_radius       Radius of quadrant
    @param[in] a_addStartVoF  If true, adds a_startVoF to the list. 
    @note This creates a quadrant and calls the other version. 
    @return Vector of VoFs in quadrant defined by the normal vector. If a_normal[dir]==0 there is no unique quadrant and the routine returns an empty vector. 
  */
  Vector<VolIndex> getAllVoFsInQuadrant(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const RealVect& a_normal, const int a_radius, const bool a_addStartVoF);

    /*!
    @brief Get all VoFs in a quadrant. The quadrant is defined by the normal vector. 
    @param[in] a_startVoF     Starting VoF
    @param[in] a_ebisBox      EBISBox
    @param[in] a_normal       Normal vector which defines the quadrant.
    @param[in] a_quadrant     Quadrant
    @param[in] a_addStartVoF  If true, adds a_startVoF to the list. 
    @return Vector of VoFs in quadrant defined by the normal vector. If a_normal[dir]==0 there is no unique quadrant and the routine returns an empty vector. 
  */
  Vector<VolIndex> getAllVoFsInQuadrant(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const IntVect a_quadrant, const int a_radius, const bool a_addStartVoF);

  /*!
    @brief Get all VoFs within a radius, but ignore VoFs that are ''opposite'' to a specified normal direction in space. 
  */
  Vector<VolIndex> getAllVoFsSymmetric(const VolIndex&      a_startVoF,
				       const EBISBox&       a_ebisbox,
				       const int            a_cardinal,
				       const Side::LoHiSide a_side,
				       const int            a_radius,
				       const bool           a_addStartVoF);

  /*!
    @brief Get all VoFs within a specific radius from a starting VoF. Can also add the starting VoF to the list
    @param[in] a_startVoF     Starting VoF
    @param[in] a_ebisBox      EBISBox
    @param[in] a_radius       Radius where to get VoFs
    @param[in] a_addStartVoF  If true, adds a_startVoF to the list. 
    @return Vector of VoFs that are within radius a_radius of a_startVoF, possibly including a_startVoF. 
  */
  Vector<VolIndex> getAllVoFsInRadius(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const int a_radius, const bool a_addStartVoF);
};

#endif
  
