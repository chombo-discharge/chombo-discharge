/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ComputationalGeometry.H
  @brief  Declaration of base class for defining geometries. 
  @author Robert Marskar
*/

#ifndef CD_ComputationalGeometry_H
#define CD_ComputationalGeometry_H

// Chombo includes
#include <BaseIF.H>
#include <MFIndexSpace.H>
#include <Box.H>
#include <RealVect.H>
#include <ProblemDomain.H>

// Our includes
#include <CD_Dielectric.H>
#include <CD_Electrode.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_RealBox.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Abstract base class for geometries
  @details This class encapsulates computational geometries in chombo-discharge. If you construct this object as-is, you will get a blank geometry. To include EBs one must
  set the electrodes and dielectrics. This is not a pure function, so you can set those objects directly from a ComputationalGeometry object. However, in almost all cases
  one will want to derive from ComputationalGeometry and create a parametrized geometry (that is what $DISCHARGE_HOME/Geometries is for!).
*/
class ComputationalGeometry
{
public:
  /*!
    @brief Constructor. Sets a blank geometry. 
  */
  ComputationalGeometry();

  /*!
    @brief Destructor
  */
  virtual ~ComputationalGeometry();

  /*!
    @brief Get dielectrics
    @return Dielectrics (m_dielectrics)
  */
  const Vector<Dielectric>&
  getDielectrics() const;

  /*!
    @brief Get electrodes
    @return Electrodes (m_electrodes)
  */
  const Vector<Electrode>&
  getElectrodes() const;

  /*!
    @brief Get the background gas permittivity
    @return Background gas permittivity
  */
  Real
  getGasPermittivity() const;

  /*!
    @brief Calls for ComputationalGeometry to use ScanShop rather than Chombo's default geometry generation tool.
    @param[in] a_startDomain Coarse domain where ScanShop begins the load balancing recursion process.
  */
  void
  useScanShop(const ProblemDomain a_startDomain);

  /*!
    @brief Calls for ComputationalGeometry to use Chombo's geometry generation tool
  */
  void
  useChomboShop();

  /*!
    @brief Set dielectrics
    @param[in] a_dielectrics Dielectris
  */
  void
  setDielectrics(const Vector<Dielectric>& a_dielectrics);

  /*!
    @brief Set electrodes
    @param[in] a_electrodes Electrodes
  */
  void
  setElectrodes(const Vector<Electrode>& a_electrodes);

  /*!
    @brief Set the background permittivity
    @param[in] a_eps0 Gas permittivity
  */
  void
  setGasPermittivity(const Real a_eps0);

  /*!
    @brief Get the multifluid index space
  */
  const RefCountedPtr<MultiFluidIndexSpace>&
  getMfIndexSpace() const;

  /*!
    @brief Get the implicit function used to generate the gas-phase EBIS
  */
  const RefCountedPtr<BaseIF>&
  getGasImplicitFunction() const;

  /*!
    @brief Get the implicit function used to generate the solid-phase EBIS
  */
  const RefCountedPtr<BaseIF>&
  getSolidImplicitFunction() const;

  /*!
    @brief Get implicit function
    @param[in] a_phase Phase
  */
  const RefCountedPtr<BaseIF>&
  getImplicitFunction(const phase::which_phase a_phase) const;

  /*!
    @brief Build geometries and the MFIndexSpace.
    @param[in] a_finestDomain Finest domain
    @param[in] a_probLo       Lower-left corner
    @param[in] a_finestDx     Finest grid resolution
    @param[in] a_nCellMax     Patch size
    @param[in] a_numEbGhost   Maximum number of EB ghosts that will be encountered.
    @param[in] a_maxCoarsen   Max coarsenings to run. If = -1 then coarsen all the way down. 
    @details This will build the gas and solid phase domain. The input is the finest-level stuff and you can control the division
    into boxes as well as the maximum number of coarsenings. The computed domain is (a_probLo, a_probLo + a_box*a_dx)
  */
  virtual void
  buildGeometries(const ProblemDomain a_finestDomain,
                  const RealVect      a_probLo,
                  const Real          a_finestDx,
                  const int           a_nCellMax,
                  const int           a_maxGhostEB,
                  const int           a_maxCoarsen = -1);

protected:
  /*!
    @brief Threshold for Vof computation
  */
  constexpr static Real s_thresh = 1.E-15;

  /*!
    @brief Multifluid index spaces
  */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /*!
    @brief Background permittivity
  */
  Real m_eps0;

  /*!
    @brief True if we use the chombo-discharge geometry generation utility.
  */
  bool m_useScanShop;

  /*!
    @brief Grid level where we begin using ScanShop
  */
  ProblemDomain m_scanDomain;

  /*!
    @brief Maximum number of ghost cells that we will ever need
  */
  int m_maxGhostEB;

  /*!
    @brief List of dielectrics
  */
  Vector<Dielectric> m_dielectrics;

  /*!
    @brief List of electrodes
  */
  Vector<Electrode> m_electrodes;

  /*!
    @brief The gas-phase implicit function (i.e. outside electrodes and dielectrics). 
  */
  RefCountedPtr<BaseIF> m_implicitFunctionGas;

  /*!
    @brief The solid-phase implicit function (i.e. the inside of the dielectrics). 
  */
  RefCountedPtr<BaseIF> m_implicitFunctionSolid;

  /*!
    @brief Set up the geometry generation tool for the gas phase.
    @param[inout] a_geoserver     Geometry service object which is later used for making build the EB information.
    @param[in]    a_finestDomain  Finest domain which will be used
    @param[in]    a_probLo        Lower-left corner of simulation domain.
    @param[in]    a_finestDx      Finest resolution which will be used.
  */
  void
  buildGasGeometry(GeometryService*&   a_geoserver,
                   const ProblemDomain a_finestDomain,
                   const RealVect      a_probLo,
                   const Real          a_finestDx);

  /*!
    @brief Set up the geometry generation tool for the solid phase, i.e. the part inside the dielectrics. 
    @param[inout] a_geoserver     Geometry service object which is later used for making build the EB information.
    @param[in]    a_finestDomain  Finest domain which will be used
    @param[in]    a_probLo        Lower-left corner of simulation domain.
    @param[in]    a_finestDx      Finest resolution which will be used.
  */
  void
  buildSolidGeometry(GeometryService*&   a_geoserver,
                     const ProblemDomain a_finestDomain,
                     const RealVect      a_probLo,
                     const Real          a_finestDx);
};

#include <CD_NamespaceFooter.H>

#endif
