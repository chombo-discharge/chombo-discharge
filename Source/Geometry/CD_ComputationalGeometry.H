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
*/
class ComputationalGeometry {
public:

  static bool s_use_new_gshop;
  static ProblemDomain s_ScanDomain;

  /*!
    @brief Constructor
  */
  ComputationalGeometry();

  /*!
    @brief Destructor
  */
  virtual ~ComputationalGeometry();

  /*!
    @brief Minimum refines
  */
  static int s_minRef;

  /*!
    @brief Maximum refines
  */
  static int s_maxRef;

  /*!
    @brief Threshold for Vof computation
  */
  static Real s_thresh;
  
  /*!
    @brief Get dielectrics
  */
  virtual const Vector<Dielectric>& getDielectrics() const;

  /*!
    @brief Get conductors
  */
  virtual const Vector<Electrode>& getElectrodes() const;

  /*!
    @brief Get the background permittivity
  */
  virtual const Real& getGasPermittivity() const;
  
  /*!
    @brief Set electrodes
  */
  virtual void setDielectrics(Vector<Dielectric>& a_dielectrics);

  /*!
    @brief Set electrodes
  */
  virtual void setElectrodes(Vector<Electrode>& a_electrodes);

  /*!
    @brief Get the background permittivity
  */
  virtual void setGasPermittivity(const Real a_eps0);

  /*!
    @brief Get the EBIndexSpaces
  */
  virtual const RefCountedPtr<MultiFluidIndexSpace>& getMfIndexSpace() const;

  /*!
    @brief Get the implicit function used to generate the gas-phase EBIS
  */
  virtual const RefCountedPtr<BaseIF>& getGasImplicitFunction() const;

  /*!
    @brief Get the implicit function used to generate the solid-phase EBIS
  */
  virtual const RefCountedPtr<BaseIF>& getSolidImplicitFunction() const;

  /*!
    @brief Get implicit function
  */
  virtual const RefCountedPtr<BaseIF>& getImplicitFunction(const phase::which_phase a_phase) const;

  /*!
    @brief Build MFIndexSpace
    @details This will build the gas and solid phase domain. The input is the finest-level stuff and you can control the division
    into boxes as well as the maximum number of coarsenings. The computed domain is (a_origin, a_origin + a_box*a_dx)
  */
  virtual void buildGeometries(const ProblemDomain   a_finestDomain,
			       const RealVect        a_origin,
			       const Real            a_finestDx,
			       const int             a_nCellMax = -1,
			       const int             a_maxCoarsen = -1);

  /*!
    @brief Build geometries from files
  */
  virtual void buildGeometriesFromFiles(const std::string&   a_gas_file,
					const std::string&   a_sol_file);

  /*!
    @brief Compute curvature of the implicit function
    @param[in] a_pos    Position where we evaluate curvature
    @param[in] a_diffDx Value to use for finite-difference evaluations. 
  */
  std::pair<Real, Real> getPrincipalCurvatures(const phase::which_phase a_phase, const RealVect pos, const Real a_diffDx) const;
  
protected:

  /*!
    @brief Multifluid index spaces
  */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /*!
    @brief Background permittivity
  */
  Real m_eps0;

  /*!
    @brief dielectrics
  */
  Vector<Dielectric> m_dielectrics;

  /*!
    @brief electrodes
  */
  Vector<Electrode> m_electrodes;

  /*!
    @brief The gas-phase implicit function
  */
  RefCountedPtr<BaseIF> m_gas_if;

  /*!
    @brief The solid-phase implicit function
  */
  RefCountedPtr<BaseIF> m_sol_if;

  /*!
    @brief Geometry server for gas phase
  */
  virtual void buildGasGeoServ(GeometryService*&   a_geoserver,
			       const ProblemDomain a_finestDomain,
			       const RealVect      a_origin,
			       const Real          a_dx);

  /*!
    @brief Geometry server for solid phas
  */
  virtual void buildSolidGeoServ(GeometryService*&   a_geoserver,
				 const ProblemDomain a_finestDomain,
				 const RealVect      a_origin,
				 const Real          a_dx);

#if CH_SPACEDIM==2
  Real curvature2D(const RefCountedPtr<BaseIF>& a_impFunc, const RealVect pos, const Real a_diffDx) const;
#elif CH_SPACEDIM==3
  std::pair<Real, Real> curvature3D(const RefCountedPtr<BaseIF>& a_impFunc, const RealVect pos, const Real a_diffDx) const;
#endif

};

#include <CD_NamespaceFooter.H>

#endif
