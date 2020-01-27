/*!
  @file computational_geometry.H
  @brief Declaration of an abstract base class for defining geometries
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _COMPUTATIONAL_GEOMETRY_
#define _COMPUTATIONAL_GEOMETRY_

#include <BaseIF.H>
#include <MFIndexSpace.H>
#include <Box.H>
#include <RealVect.H>
#include <ProblemDomain.H>

#include "dielectric.H"
#include "electrode.H"
#include "physical_domain.H"
#include "mfis.H"
#include "real_box.H"

/*!
  @brief Abstract base class for geometries
*/
class computational_geometry {
public:

  static bool s_use_new_gshop;
  static ProblemDomain s_ScanDomain;

  /*!
    @brief Constructor
  */
  computational_geometry();

  /*!
    @brief Destructor
  */
  virtual ~computational_geometry();

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
  virtual const Vector<dielectric>& get_dielectrics() const;

  /*!
    @brief Get conductors
  */
  virtual const Vector<electrode>& get_electrodes() const;

  /*!
    @brief Get the background permittivity
  */
  virtual const Real& get_eps0() const;
  
  /*!
    @brief Set electrodes
  */
  virtual void set_dielectrics(Vector<dielectric>& a_dielectrics);

  /*!
    @brief Set electrodes
  */
  virtual void set_electrodes(Vector<electrode>& a_electrodes);

  /*!
    @brief Get the background permittivity
  */
  virtual void set_eps0(const Real a_eps0);

  /*!
    @brief Get the EBIndexSpaces
  */
  virtual const RefCountedPtr<mfis>& get_mfis() const;

  /*!
    @brief Get the implicit function used to generate the gas-phase EBIS
  */
  virtual const RefCountedPtr<BaseIF>& get_gas_if() const;

  /*!
    @brief Get the implicit function used to generate the solid-phase EBIS
  */
  virtual const RefCountedPtr<BaseIF>& get_sol_if() const;

  /*!
    @brief Build MFIndexSpace
    @details This will build the gas and solid phase domain. The input is the finest-level stuff and you can control the division
    into boxes as well as the maximum number of coarsenings. The computed domain is (a_origin, a_origin + a_box*a_dx)
  */
  virtual void build_geometries(const physical_domain& a_physDom,
				const ProblemDomain&   a_finestDomain,
				const Real&            a_finestDx,
  				const int&             a_nCellMax = -1,
  				const int&             a_maxCoarsen = -1);

  /*!
    @brief Build geometries from files
  */
  virtual void build_geo_from_files(const std::string&   a_gas_file,
				    const std::string&   a_sol_file);

protected:

  /*!
    @brief Multifluid index spaces
  */
  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief Background permittivity
  */
  Real m_eps0;

  /*!
    @brief dielectrics
  */
  Vector<dielectric> m_dielectrics;

  /*!
    @brief electrodes
  */
  Vector<electrode> m_electrodes;

  /*!
    @brief The gas-phase implicit function
  */
  RefCountedPtr<BaseIF> m_gas_if;

  /*!
    @brief The solid-phase implicit function
  */
  RefCountedPtr<BaseIF> m_sol_if;

  // Voxels for gas phase
  Vector<real_box> m_regular_voxels_gas;
  Vector<real_box> m_covered_voxels_gas;
  Vector<real_box> m_bounded_voxels_gas;

  // Voxels for solid phase
  Vector<real_box> m_regular_voxels_sol;
  Vector<real_box> m_covered_voxels_sol;
  Vector<real_box> m_bounded_voxels_sol;

  /*!
    @brief Geometry server for gas phase
  */
  virtual void build_gas_geoserv(GeometryService*&    a_geoserver,
				 const ProblemDomain& a_finestDomain,
				 const RealVect&      a_origin,
				 const Real&          a_dx);

  /*!
    @brief Geometry server for solid phas
  */
  virtual void build_solid_geoserv(GeometryService*&     a_geoserver,
				   const ProblemDomain& a_finestDomain,
				   const RealVect&      a_origin,
				   const Real&          a_dx);

};
#endif
