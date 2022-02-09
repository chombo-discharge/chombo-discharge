/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DischargeIO.H
  @brief  Silly, but useful functions that override standard Chombo HDF5 IO.
  @author Robert Marskar
*/

#ifndef CD_DischargeIO_H
#define CD_DischargeIO_H

// Chombo includes
#include <REAL.H>
#include <RealVect.H>
#include <IntVect.H>
#include <Vector.H>
#include <DisjointBoxLayout.H>
#include <LevelData.H>
#include <EBCellFAB.H>
#include <ProblemDomain.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Namespace which encapsulates some changes to standard Chombo HDF5 IO.
*/
namespace DischargeIO {

  /*!
    @brief A shameless copy of Chombo's writeEBHDF5 but including the lower-left corner of the physical domain as well.
    @param[in] a_filename         File name
    @param[in] a_grids            Grids on each AMR level
    @param[in] a_data             Data on the AMR levels
    @param[in] a_variableNames    Grid variable names
    @param[in] a_coarsestDomain   Coarsest grid domain
    @param[in] a_coarsestDx       Coarsest grid resolution
    @param[in] a_dt               Time step
    @param[in] a_time             Time
    @param[in] a_probLo           Lower-left physical corner
    @param[in] a_refinementRatios Refinement ratios
    @param[in] a_numLevels        Number of AMR levels
    @param[in] a_ghostVect        Number of ghost cells to include in the plot files.
    @todo The current implementation doubles the memory in the output data. We should look into more general aliasing functions. 
   */
  void writeEBHDF5(const std::string&                     a_filename,
		   const Vector<DisjointBoxLayout>&       a_grids,
		   const Vector<LevelData<EBCellFAB>* > & a_data,
		   const Vector<std::string>&             a_variableNames,
		   const ProblemDomain&                   a_coarsestDomain,
		   const Real&                            a_coarsestDx,
		   const Real&                            a_dt,
		   const Real&                            a_time,
		   const RealVect&                        a_probLo,
		   const Vector<int>&                     a_refinementRatios,
		   const int&                             a_numLevels,
		   const IntVect                          a_ghostVect= IntVect::Zero);  

}

#include <CD_NamespaceFooter.H>

#endif
