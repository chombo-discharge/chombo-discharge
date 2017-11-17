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
#include "wall_bc.H"

/*!
  @brief Global stuff for phases
*/
namespace Phase {
  const int num_phases = 2;
  
  enum WhichPhase {
    Gas = 0,
    Solid = 1
  };
};


/*!
  @brief Abstract base class for geometries
*/
class computational_geometry {
public:

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
    @brief Set Dirichlet wall BC
  */
  virtual void set_dirichlet_wall_bc(const int a_dir, Side::LoHiSide a_side, const bool a_live);

  /*!
    @brief Set Neumann wall BC
  */
  virtual void set_neumann_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value);

  /*!
    @brief Get Poisson BC at specified domain wall
  */
  virtual const wall_bc& get_wall_bc(const int a_dir, Side::LoHiSide a_side) const;

  /*!
    @brief Get the background permittivity
  */
  virtual void set_eps0(const Real a_eps0);

  /*!
    @brief Get the EBIndexSpaces
  */
  virtual const RefCountedPtr<MFIndexSpace> get_ebis_mf() const;

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

protected:

  /*!
    @brief MF index space
  */
  RefCountedPtr<MFIndexSpace> m_ebis_mf;

  /*!
    @brief Background permittivity
  */
  Real m_eps0;

  /*!
    @brief Wall boundary conditions
  */
  Vector<RefCountedPtr<wall_bc> > m_wallbc;

  /*!
    @brief dielectrics
  */
  Vector<dielectric> m_dielectrics;

  /*!
    @brief electrodes
  */
  Vector<electrode> m_electrodes;

  /*!
    @brief Geometry server for gas phase
  */
  virtual void build_gas_geoserv(GeometryService*&     a_geoserver,
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

  /*!
    @brief Sanity check
  */
  virtual void sanity_check();

  /*!
    @brief Internal bc map
  */
  int map_bc(const int a_dir, const Side::LoHiSide a_side) const;

};
#endif
