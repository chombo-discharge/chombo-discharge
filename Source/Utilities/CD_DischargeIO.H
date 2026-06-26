/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_DischargeIO.H
  @brief  Silly, but useful functions that override standard Chombo HDF5 IO.
  @author Robert Marskar
*/

#ifndef CD_DISCHARGEIO_H
#define CD_DISCHARGEIO_H

// Std includes
#include <string>
#include <vector>
#include <utility>
#include <functional>

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
#include <CD_ParticleSoA.H>
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Namespace which encapsulates chombo-discharge IO functionality. 
*/
namespace DischargeIO {
  /**
    @brief Number formatting method -- writes big numbers using an input separator. E.g. the number 123456 is written as 123,456
    @param[in] n Number to format as string with separator
    @param[in] a_sep    Separator
    @return Return value
  */
  std::string
  numberFmt(long long n, char a_sep = ',') noexcept;

  /**
    @brief Number formatting game. This just prints the number 123456 as 123,456
    @param[in] a_numbers Numbers to format as string with separator
    @param[in] a_sep     Separator
    @return Return value
  */
  Vector<std::string>
  numberFmt(const Vector<long long>& a_numbers, char a_sep = ',') noexcept;

#ifdef CH_USE_HDF5
  /**
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
  void
  writeEBHDF5(const std::string&                   a_filename,
              const Vector<std::string>&           a_variableNames,
              const Vector<DisjointBoxLayout>&     a_grids,
              const Vector<LevelData<EBCellFAB>*>& a_data,
              const Vector<ProblemDomain>&         a_domains,
              const Vector<Real>&                  a_dx,
              const Vector<int>&                   a_refinementRatios,
              Real                                 a_dt,
              Real                                 a_time,
              const RealVect&                      a_probLo,
              int                                  a_numLevels,
              int                                  a_numGhost) noexcept;
#endif

#ifdef CH_USE_HDF5
  /**
    @brief Write HDF5 header to the input handle
    @param[in,out] a_handleH5 Handle to HDF5 file
    @param[in] a_numLevels Number of levels to write
    @param[in] a_probLo Lower-left corner of problem domain
    @param[in] a_variableNames Variable names
  */
  void
  writeEBHDF5Header(HDF5Handle&                a_handleH5,
                    int                        a_numLevels,
                    const RealVect&            a_probLo,
                    const Vector<std::string>& a_variableNames) noexcept;

#endif

#ifdef CH_USE_HDF5
  /**
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
                   const ProblemDomain&        a_domain,
                   Real                        a_dx,
                   Real                        a_dt,
                   Real                        a_time,
                   int                         a_level,
                   int                         a_refRatio,
                   int                         a_numGhost) noexcept;
#endif

#ifdef CH_USE_HDF5
  /**
    @brief Debugging function for quickly writing EBAMRCellData to HDF5
    @param[in] a_data Cell-centered EB+AMR data
    @param[in] a_file File name
  */
  void
  writeEBHDF5(const EBAMRCellData& a_data, const std::string& a_file);
#endif

  /**
    @brief Accessor returning a scalar payload quantity of SoA particle (a_leaf, a_index) for writeH5Part.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  using H5PartScalar = std::function<Real(const ParticleSoA<P, Traits>&, std::size_t)>;

  /**
    @brief Accessor returning a vector payload quantity of SoA particle (a_leaf, a_index) for writeH5Part.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  using H5PartVect = std::function<RealVect(const ParticleSoA<P, Traits>&, std::size_t)>;

  /**
    @brief Write an SoA particle container to an H5Part file (quick visualization).
    @details Position, id, and weight are taken from the container-owned columns; the caller supplies the
    payload datasets as (name, accessor) pairs so the
    scalar-vs-vector distinction and naming are explicit (the SoA columns tuple interleaves scalar and
    vector-component columns, so there is no implicit reals/vects partition to rely on). Each vector
    dataset is written as name-x/name-y/name-z.
    @param[in] a_filename    File name.
    @param[in] a_particles   SoA particle container.
    @param[in] a_scalarVars  Scalar datasets: (name, accessor) pairs.
    @param[in] a_vectorVars  Vector datasets: (name, accessor) pairs.
    @param[in] a_shift       Particle position shift.
    @param[in] a_time        Physical time.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  writeH5Part(std::string                                                  a_filename,
              const ParticleContainer<P, Traits>&                          a_particles,
              std::vector<std::pair<std::string, H5PartScalar<P, Traits>>> a_scalarVars,
              std::vector<std::pair<std::string, H5PartVect<P, Traits>>>   a_vectorVars,
              RealVect                                                     a_shift,
              Real                                                         a_time) noexcept;

#ifdef CH_USE_HDF5
  /**
    @brief Write checkpointed SoA particles (one level) to HDF5.
    @details SoA analogue of writeCheckParticlesToHDF(ParticleData<P>): the per-particle HDF5 record is
    position + weight + the payload h5-column subset (all doubles), via ParticleSoA::h5LinearizeParticle.
    @param[in,out] a_handle    Handle to HDF5 file.
    @param[in]     a_particles Per-level SoA particles to checkpoint.
    @param[in]     a_dataType  String identifier.
    @tparam P      Payload type of the particles.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  writeCheckParticlesToHDF(HDF5Handle&                               a_handle,
                           const LayoutData<ParticleSoA<P, Traits>>& a_particles,
                           const std::string&                        a_dataType);

  /**
    @brief Read checkpointed SoA particles (one level) from HDF5.
    @details SoA analogue of readCheckParticlesFromHDF(ParticleData<P>): each record is delinearized and
    appended via ParticleSoA::h5DelinearizeAndAppend. The destination layout must match the checkpoint.
    @param[in,out] a_handle    Handle to HDF5 file.
    @param[in,out] a_particles Per-level SoA particles to fill.
    @param[in]     a_dataType  String identifier.
    @tparam P      Payload type of the particles.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  void
  readCheckParticlesFromHDF(HDF5Handle&                         a_handle,
                            LayoutData<ParticleSoA<P, Traits>>& a_particles,
                            const std::string&                  a_dataType);
#endif
} // namespace DischargeIO

#include <CD_NamespaceFooter.H>

#include <CD_DischargeIOImplem.H>

#endif
