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
#include <CD_ParticleContainer.H>
#include <CD_GenericParticle.H>
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
    @details User are not supported to call this. Driver will perform calls to writeEBHDF5Header and writeEBHDF5Level
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
              const int                            a_numGhost) noexcept;
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Write HDF5 header to the input handle
    @param[inout] a_handleH5 Handle to HDF5 file
    @param[in] a_numLevels Number of levels to write
    @param[in] a_probLo Lower-left corner of problem domain
    @param[in] a_variableNames Variable names
  */
  void
  writeEBHDF5Header(HDF5Handle&                a_handleH5,
                    const int                  a_numLevels,
                    const RealVect&            a_probLo,
                    const Vector<std::string>& a_variableNames) noexcept;

#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Write data to output handle. Meant to be called on each level after calling writeEBHDF5 header.
    @param[in] a_handleH5 Handle to HDF5 data
    @param[in] a_outputData Data to write
    @param[in] a_domain Problem domain
    @param[in] a_dx Grid resolution
    @param[in] a_dt Time step
    @param[in] a_time Time
    @param[in] a_level AMR level
    @param[in] a_refRatio Refinement ratio
    @param[in] a_numGhost Number of ghost cells to fill. 
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

  /*!
    @brief Write a particle container to an H5Part file. Good for quick and dirty visualization of particles
    @details Use case is pretty straightforward but the user might need to cast particle types. E.g. call

    writeH5Part<M,N>(a_filename, (const ParticleContainer<GenericParticle<M,N>>&) a_particles, ...)

    Template substitution is not straightforward for this one. 
    @param[in] a_filename  File name
    @param[in] a_particles Particles. Particle type must derive from GenericParticle<M, N>
    @param[in] a_realVars  Variable names for the M real variables
    @param[in] a_vectVars  Variable names for the N vector variables
    @param[in] a_shift     Particle position shift
    @param[in] a_shift     Time 
  */
  template <size_t M, size_t N>
  void
  writeH5Part(const std::string                               a_filename,
              const ParticleContainer<GenericParticle<M, N>>& a_particles,
              const std::vector<std::string>                  a_realVars = std::vector<std::string>(),
              const std::vector<std::string>                  a_vectVars = std::vector<std::string>(),
              const RealVect                                  a_shift    = RealVect::Zero,
              const Real                                      a_time     = 0.0) noexcept;

#ifdef CH_USE_HDF5
  /*!
    @brief Write checkpointed particles to HDF5
    @param[in, out] a_handle Handle to HDF5 files
    @param[in] a_particles Particles to checkpoint
    @param[in] a_dataType String identifier
    @note This is basically a bootleg of Chombo's HDF5 API, but replacing P::size() by P::H5size()
  */
  template <class P>
  void
  writeCheckParticlesToHDF(HDF5Handle& a_handle, const ParticleData<P>& a_particles, const std::string& a_dataType);

  /*!
    @brief Write checkpointed particles to HDF5
    @param[in, out] a_handle Handle to HDF5 files
    @param[in] a_particles Particles to checkpoint
    @param[in] a_dataType String identifier
    @note This is basically a bootleg of Chombo's HDF5 API, but replacing P::size() by P::H5size()    
  */
  template <class P>
  void
  readCheckParticlesFromHDF(HDF5Handle& a_handle, ParticleData<P>& a_particles, const std::string& a_dataType);
#endif
} // namespace DischargeIO

#include <CD_NamespaceFooter.H>

#include <CD_DischargeIOImplem.H>

#endif
