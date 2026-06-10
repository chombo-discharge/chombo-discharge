/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MultiFluidIndexSpace.H
  @brief  Multi-fluid index space
  @author Robert Marskar
*/

#ifndef CD_MULTIFLUIDINDEXSPACE_H
#define CD_MULTIFLUIDINDEXSPACE_H

// Chombo includes
#include <GeometryService.H>
#include <EBIndexSpace.H>
#include <RealVect.H>
#include <MFIndexSpace.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Phase names namespace.
*/
namespace phase {
  /**
    @brief Number of supported phases.
  */
  const int numPhases = 2;

  /**
    @brief Enumeration of supported phases.
  */
  enum which_phase
  {
    gas   = 0, ///< Gas phase.
    solid = 1  ///< Solid (dielectric) phase.
  };
} // namespace phase

/**
  @brief Two-fluid index space.
  @details This is basically a wrapper for using two Chombo EBIndexSpaces in the simulation. It calls Chombo for building the EBIS'.
*/
class MultiFluidIndexSpace
{
public:
  /**
    @brief Weak constructor. Must call define later
  */
  MultiFluidIndexSpace();

  /**
    @brief Destructor
  */
  virtual ~MultiFluidIndexSpace();

  /**
    @brief Define function
    @param[in] a_domain          Finest domain
    @param[in] a_origin          Lower-left corner
    @param[in] a_dx              Finest resolution
    @param[in] a_geoservers      Generators for cut-cells
    @param[in] a_distributedData Special flag for Chombo (passing grid generation to the geometryservice)
    @param[in] a_nCellMax        Minimum number of cells in the coarsening
    @param[in] a_maxCoarsenings  Maximum number of coarsenings.
    @param[in] a_fixOnlyFirstPhaseRegNextToMultiValued Special flag for Chombo.
  */
  virtual void
  define(const Box&                      a_domain,
         const RealVect&                 a_origin,
         const Real&                     a_dx,
         const Vector<GeometryService*>& a_geoservers,
         const bool                      a_distributedData,
         int                             a_nCellMax                              = -1,
         int                             a_max_coar                              = -1,
         bool                            a_fixOnlyFirstPhaseRegNextToMultiValued = false);

  /**
    @brief Get a particular EBIndexSpace
    @param[in] a_phase Phase
    @return Const reference to the EBIndexSpace for the given phase.
  */
  const RefCountedPtr<EBIndexSpace>&
  getEBIndexSpace(const phase::which_phase a_phase) const;

  /**
    @brief Get a particular EBIndexSpace
    @param[in] a_phase Phase
    @return Const reference to the EBIndexSpace for the given phase.
  */
  const RefCountedPtr<EBIndexSpace>&
  getEBIndexSpace(const int a_phase) const;

  /**
    @brief Get a particular EBIndexSpace
    @param[in] a_phase Phase
    @return Reference to the EBIndexSpace for the given phase.
  */
  RefCountedPtr<EBIndexSpace>&
  getEBIndexSpace(const phase::which_phase a_phase);

  /**
    @brief Get EBIndexSpace
    @param[in] a_phase Phase
    @return Reference to the EBIndexSpace for the given phase.
  */
  RefCountedPtr<EBIndexSpace>&
  getEBIndexSpace(const int a_phase);

  /**
    @brief Get the interface region. This is defined as the intersection between the irregular cells of the
    gas phase and the other phases.
    @param[in] a_domain Domain
    @return Interface region as an IntVectSet.
  */
  IntVectSet
  interfaceRegion(const ProblemDomain& a_domain) const;

  /**
    @brief Get number of phases
    @return Number of phases.
  */
  static int
  numPhases() ;

protected:
  /**
    @brief EBIndexSpaces
  */
  Vector<RefCountedPtr<EBIndexSpace>> m_ebis;

  /**
    @brief MFIndexSpace
  */
  RefCountedPtr<MFIndexSpace> m_multifluidIndexSpace;
};

#include <CD_NamespaceFooter.H>

#endif
