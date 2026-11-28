/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CdrPlasmaGodunovStorage.H
 * @brief  Declaration of scratch storage for CdrPlasmaGodunovStepper
 * @author Robert Marskar
 */

#ifndef CD_CDRPLASMAGODUNOVSTORAGE_H
#define CD_CDRPLASMAGODUNOVSTORAGE_H

// Our includes
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
namespace CdrPlasma {

/**
 * @brief Per-step storage for one CDR species in CdrPlasmaGodunovStepper.
 * @details The constructor fetches the temporaries from the scratch pool of AmrMesh and the destructor returns them.
 * Objects of this class are created by CdrPlasmaGodunovStepper::advance and do not outlive it, so they never hold
 * temporaries across a regrid.
 */
class CdrPlasmaGodunovStepper::CdrStorage
{
public:
  /**
   * @brief Disallowed constructor. Only strong construction allowed.
   */
  CdrStorage() = delete;

  /**
   * @brief Full constructor. Fetches the temporaries from the scratch pool.
   * @param[in] a_amr   Reference to AmrMesh
   * @param[in] a_realm Realm where storage lives
   * @param[in] a_phase Phase where storage lives (only single-phased data in here)
   */
  CdrStorage(const RefCountedPtr<AmrMesh>& a_amr, const std::string& a_realm, phase::which_phase a_phase);

  /**
   * @brief Destructor. Returns the temporaries to the scratch pool.
   */
  virtual ~CdrStorage();

  CdrStorage(const CdrStorage&) = delete;
  CdrStorage&
  operator=(const CdrStorage&) = delete;
  CdrStorage(CdrStorage&&)     = delete;
  CdrStorage&
  operator=(CdrStorage&&) = delete;

  /**
   * @brief Get scratch storage
   * @return Scratch
   */
  EBAMRCellData&
  getScratch() const
  {
    return *m_scratch;
  }

  /**
   * @brief Get the other scratch storage.
   * @return Scratch2
   */
  EBAMRCellData&
  getScratch2() const
  {
    return *m_scratch2;
  }

  /**
   * @brief Get the extrapolation storage. This is usually used for holding phi + 0.5*dt*S
   * @return Extrap
   */
  EBAMRCellData&
  getExtrap() const
  {
    return *m_cellExtr;
  }

  /**
   * @brief Get storage for holding gradient.
   * @return Gradient
   */
  EBAMRCellData&
  getGradient() const
  {
    return *m_gradient;
  }

  /**
   * @brief Get scratch storage that holds the density on the EB
   * @return Eb state
   */
  EBAMRIVData&
  getEbState() const
  {
    return *m_scratchIV1;
  }

  /**
   * @brief Get the scratch storage holding the normal component of the velocity on the EB
   * @return Eb velo
   */
  EBAMRIVData&
  getEbVelo() const
  {
    return *m_scratchIV2;
  }

  /**
   * @brief Get the scratch storage holding the extrapolated flux on the EB
   * @return Eb flux
   */
  EBAMRIVData&
  getEbFlux() const
  {
    return *m_scratchIV3;
  }

  /**
   * @brief Get the scratch storage holding the extrapolated gradient on the EB
   * @return Eb grad
   */
  EBAMRIVData&
  getEbGrad() const
  {
    return *m_scratchIV4;
  }

  /**
   * @brief Get the scratch storage holding the extrapolated density on domain faces
   * @return Domain state
   */
  EBAMRIFData&
  getDomainState() const
  {
    return *m_scratchIF1;
  }

  /**
   * @brief Get the scratch storage holding the extrapolated velocity (normal component only) on domain faces
   * @return Domain velo
   */
  EBAMRIFData&
  getDomainVelo() const
  {
    return *m_scratchIF2;
  }

  /**
   * @brief Get the scratch storage holding the extrapolated flux on domain faces
   * @return Domain flux
   */
  EBAMRIFData&
  getDomainFlux() const
  {
    return *m_scratchIF3;
  }

  /**
   * @brief Get the scratch storage holding the extrapolated gradient (normal component only) on domain faces
   * @return Domain grad
   */
  EBAMRIFData&
  getDomainGrad() const
  {
    return *m_scratchIF4;
  }

protected:
  /**
   * @brief Cell-centered scratch data
   */
  mutable Scratch<EBAMRCellData> m_scratch;

  /**
   * @brief Cell-centered scratch data
   */
  mutable Scratch<EBAMRCellData> m_scratch2;

  /**
   * @brief Cell-centered scratch data. Used for holding phi + 0.5*dt*S
   */
  mutable Scratch<EBAMRCellData> m_cellExtr;

  /**
   * @brief Cell-centered scratch data for holding the gradient.
   */
  mutable Scratch<EBAMRCellData> m_gradient;

  /**
   * @brief EB-centered scratch data
   */
  mutable Scratch<EBAMRIVData> m_scratchIV1;

  /**
   * @brief EB-centered scratch data
   */
  mutable Scratch<EBAMRIVData> m_scratchIV2;

  /**
   * @brief EB-centered scratch data
   */
  mutable Scratch<EBAMRIVData> m_scratchIV3;

  /**
   * @brief EB-centered scratch data
   */
  mutable Scratch<EBAMRIVData> m_scratchIV4;

  /**
   * @brief Domain-centered scratch data
   */
  mutable Scratch<EBAMRIFData> m_scratchIF1;

  /**
   * @brief Domain-centered scratch data
   */
  mutable Scratch<EBAMRIFData> m_scratchIF2;

  /**
   * @brief Domain-centered scratch data
   */
  mutable Scratch<EBAMRIFData> m_scratchIF3;

  /**
   * @brief Domain-centered scratch data
   */
  mutable Scratch<EBAMRIFData> m_scratchIF4;
};

/**
 * @brief Per-step storage for the electric field in CdrPlasmaGodunovStepper.
 * @details The constructor fetches the temporaries from the scratch pool of AmrMesh and the destructor returns them.
 * Objects of this class are created by CdrPlasmaGodunovStepper::advance and do not outlive it, so they never hold
 * temporaries across a regrid.
 */
class CdrPlasmaGodunovStepper::FieldStorage
{
public:
  /**
   * @brief Disallowed weak constructor
   */
  FieldStorage() = delete;

  /**
   * @brief Full constructor. Fetches the temporaries from the scratch pool.
   * @param[in] a_amr   Reference to AmrMesh
   * @param[in] a_realm Realm where storage lives
   * @param[in] a_phase Phase where storage lives (only single-phased data in here)
   */
  FieldStorage(const RefCountedPtr<AmrMesh>& a_amr, const std::string& a_realm, phase::which_phase a_phase);

  /**
   * @brief Destructor. Returns the temporaries to the scratch pool.
   */
  virtual ~FieldStorage();

  FieldStorage(const FieldStorage&) = delete;
  FieldStorage&
  operator=(const FieldStorage&) = delete;
  FieldStorage(FieldStorage&&)   = delete;
  FieldStorage&
  operator=(FieldStorage&&) = delete;

  /**
   * @brief Get scratch storage that holds the electric field at the cell center
   * @return Electric field cell
   */
  EBAMRCellData&
  getElectricFieldCell() const
  {
    return *m_electricFieldCell;
  }

  /**
   * @brief Get scratch storage that holds the electric field at the EB
   * @return Electric field e b
   */
  EBAMRIVData&
  getElectricFieldEB() const
  {
    return *m_electricFieldEB;
  }

  /**
   * @brief Get scratch storage that holds the electric field at the domainf aces
   * @return Electric field domain
   */
  EBAMRIFData&
  getElectricFieldDomain() const
  {
    return *m_electricFieldDomain;
  }

protected:
  /**
   * @brief Scratch storage for holding electric field at cell center.
   */
  mutable Scratch<EBAMRCellData> m_electricFieldCell;

  /**
   * @brief Scratch storage for holding electric field at the EB.
   */
  mutable Scratch<EBAMRIVData> m_electricFieldEB;

  /**
   * @brief Scratch storage for holding electric field at domain faces.
   */
  mutable Scratch<EBAMRIFData> m_electricFieldDomain;
};

/**
 * @brief Per-step storage for one RTE species in CdrPlasmaGodunovStepper.
 * @details The constructor fetches the temporaries from the scratch pool of AmrMesh and the destructor returns them.
 * Objects of this class are created by CdrPlasmaGodunovStepper::advance and do not outlive it, so they never hold
 * temporaries across a regrid.
 */
class CdrPlasmaGodunovStepper::RtStorage
{
public:
  /**
   * @brief Disallowed weak constructor
   */
  RtStorage() = delete;

  /**
   * @brief Full constructor. Fetches the temporaries from the scratch pool.
   * @param[in] a_amr   Reference to AmrMesh
   * @param[in] a_realm Realm where storage lives
   * @param[in] a_phase Phase where storage lives (only single-phased data in here)
   */
  RtStorage(const RefCountedPtr<AmrMesh>& a_amr, const std::string& a_realm, phase::which_phase a_phase);

  /**
   * @brief Destructor. Returns the temporaries to the scratch pool.
   */
  virtual ~RtStorage();

  RtStorage(const RtStorage&) = delete;
  RtStorage&
  operator=(const RtStorage&) = delete;
  RtStorage(RtStorage&&)      = delete;
  RtStorage&
  operator=(RtStorage&&) = delete;

  /**
   * @brief Get scratch storage that holds the EB flux
   * @return Eb flux
   */
  EBAMRIVData&
  getEbFlux() const
  {
    return *m_scratchIV;
  }

  /**
   * @brief Get scratch storage that holds the domain
   * @return Domain flux
   */
  EBAMRIFData&
  getDomainFlux() const
  {
    return *m_scratchIF;
  }

protected:
  /**
   * @brief Scratch storage on EB
   */
  mutable Scratch<EBAMRIVData> m_scratchIV;

  /**
   * @brief Scratch storage on domain faces
   */
  mutable Scratch<EBAMRIFData> m_scratchIF;
};

} // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
