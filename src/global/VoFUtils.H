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
    @brief Get quadrant defined by direction in space. 
    @param[in] a_normal The direction in space which defines the quadrant. 
    @return Quadrant as IntVect. The IntVect is ill-defined if the quadrant is not well-defined. 
  */
  IntVect getQuadrant(const RealVect a_normal);

  /*!
    @brief Get all VoFs in a quadrant. The quadrant is defined by the normal vector. 
    @param[in] a_startVoF     Starting VoF
    @param[in] a_ebisBox      EBISBox
    @param[in] a_normal       Normal vector which defines the quadrant.
    @param[in] a_radius       Radius of quadrant
    @param[in] a_addStartVoF  If true, adds a_startVoF to the list. 
    @return Vector of VoFs in quadrant defined by the normal vector. If a_normal[dir]==0 there is no unique quadrant and the routine returns an empty vector. 
  */
  Vector<VolIndex> getAllVoFsInQuadrant(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const RealVect a_normal, const int a_radius, const bool a_addStartVoF);

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
  
