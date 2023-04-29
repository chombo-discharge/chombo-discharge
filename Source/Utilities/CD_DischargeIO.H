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

// Std includes
#include <string>

// Chombo includes
#include <REAL.H>
#include <RealVect.H>
#include <IntVect.H>
#include <Vector.H>
#include <DisjointBoxLayout.H>
#include <LevelData.H>
#include <EBCellFAB.H>
#include <ParticleIO.H>
#include <ProblemDomain.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Namespace which encapsulates chombo-discharge IO functionality. 
*/
namespace DischargeIO {
  /*!
    @brief Number formatting method -- writes big numbers using an input separator. E.g. the number 123456 is written as 123,456
    @param[in] a_number Number to format as string with separator
    @param[in] a_sep    Separator
  */
  std::string
  numberFmt(const long long a_number, char a_sep = ',') noexcept;

  /*!
    @brief Number formatting game. This just prints the number 123456 as 123,456
    @param[in] a_numbers Numbers to format as string with separator
    @param[in] a_sep     Separator
  */
  Vector<std::string>
  numberFmt(const Vector<long long> a_numbers, char a_sep = ',') noexcept;

  /*!
    @brief A shameless copy of Chombo's writeEBHDF5 but including the lower-left corner of the physical domain as well.
    @param[in] a_filename         File name
    @param[in] a_variableNames    Grid variable names
    @param[in] a_grids            Grids on each AMR level
    @param[in] a_data             Data on the AMR levels
    @param[in] a_domains          Grid domains
    @param[in] a_dx               Grid resolutions
    @param[in] a_refinementRatios Refinement ratios
    @param[in] a_dt               Time step
    @param[in] a_time             Time
    @param[in] a_probLo           Lower-left physical corner
    @param[in] a_numLevels        Number of AMR levels
    @param[in] a_numGhost         Number of ghost cells included in the output. 
    @todo The current implementation doubles the memory in the output data. We should look into more general aliasing functions. 
   */
#ifdef CH_USE_HDF5
  void
  writeEBHDF5(const std::string&                   a_filename,
              const Vector<std::string>&           a_variableNames,
              const Vector<DisjointBoxLayout>&     a_grids,
              const Vector<LevelData<EBCellFAB>*>& a_data,
              const Vector<ProblemDomain>&         a_domains,
              const Vector<Real>                   a_dx,
              const Vector<int>                    a_refinementRatios,
              const Real                           a_dt,
              const Real                           a_time,
              const RealVect                       a_probLo,
              const int                            a_numLevels,
              const int                            a_numGhost);
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Write HDF5 header to the input handle
  */
  void
  writeEBHDF5Header(HDF5Handle&                a_handleH5,
                    const int                  a_numLevels,
                    const RealVect&            a_probLo,
                    const Vector<std::string>& a_variableNames) noexcept;

#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Write data to output handle. 
  */
  void
  writeEBHDF5Level(HDF5Handle&                 a_handleH5,
                   const LevelData<EBCellFAB>& a_outputData,
                   const ProblemDomain         a_domain,
                   const Real                  a_dx,
                   const Real                  a_dt,
                   const Real                  a_time,
                   const int                   a_level,
                   const int                   a_refRatio,
                   const int                   a_numGhost) noexcept;
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Debugging function for quickly writing EBAMRCellData to HDF5
    @param[in] a_data Cell-centered EB+AMR data
    @param[in] a_file File name
  */
  void
  writeEBHDF5(const EBAMRCellData& a_data, const std::string& a_file);
#endif

} // namespace DischargeIO

#include <CD_NamespaceFooter.H>

#endif
