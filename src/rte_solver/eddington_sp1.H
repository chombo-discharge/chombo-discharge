/*!
  @file eddington_sp1.H
  @brief Declaration of a first order Eddington solver
  @author Robert Marskar
  @date Jan. 2018
*/

#ifndef _EDDINGTON_SP1_
#define _EDDINGTON_SP1_

#include "rte_solver.H"
#include "ebconductivityopfactory.H"
#include "larsen_coefs.H"
#include "robinconductivityebbcfactory.H"
#include "robinconductivitydomainbcfactory.H"
#include "conductivitydomainbc_wrapper_factory.H"
#include "wall_bc.H"

#include <AMRMultiGrid.H>
#include <AMRTGA.H>
#include <BiCGStabSolver.H>
#include <EBBackwardEuler.H>
#include <EBSimpleSolver.H>
#include <EBConductivityOpFactory.H>

#include <random>

#include "CD_NamespaceHeader.H"

/*!
  @brief Radiative tranfer equation solver in the SP1 (diffusion) approximation. Also includes stochastically drawn photons. 
*/
class eddington_sp1 : public rte_solver {
protected:

  /*!
    @brief Silly enum for switching between relaxation methods.
  */
  enum class relax{
    jacobi,
    gauss_seidel,
    gsrb_fast
  };

  /*!
    @brief Silly enum for switching between multigrid cycling
  */
  enum class amrmg{
    full,
    vcycle,
    fcycle,
  };
  
public:

  /*!
    @brief Constructor
  */
  eddington_sp1();

  /*!
    @brief Destructor
  */
  virtual ~eddington_sp1();

  /*!
    @brief Advance RTE onto state a_phi
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zerophi = false);

  /*!
    @brief Parse class options
  */
  virtual void parseOptions();

  /*!
    @brief Parse class options
  */
  virtual void parseRuntimeOptions();
  
  /*!
    @brief Allocate internal storage
  */
  virtual void allocateInternals();

  /*!
    @brief Allocate wall bc
  */
  virtual void allocate_wall_bc();
  
  /*!
    @brief Cache state
  */
  virtual void preRegrid(const int a_base, const int a_oldFinestLevel);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocateInternals();

  /*!
    @brief Regrid function for this class
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Register operators
  */
  virtual void registerOperators();

  /*!
    @brief Set reflection coefficients. These must be obtained by integrals over the reflectivity.
  */
  virtual void set_reflection_coefficients(const Real r1, const Real r2);

  /*!
    @brief Compute the boundary flux. For Eddington, the boundary flux is = c*phi/2
  */
  virtual void compute_boundary_flux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi);

  /*!
    @brief Compute the domain flux. For Eddington, the domain flux is = c*phi/2
  */
  virtual void compute_domain_flux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_phi);

  /*!
    @brief Compute the flux
    @details For diffusive models, the flux will be something like grad(a_phi). For higher-order models, the flux 
    will be contained in a_phi (somehow);
  */
  virtual void computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi);

  /*!
    @brief Get isotropic part. 
  */
  virtual void compute_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi);

  /*!
    @brief Write plot file
  */
  virtual void writePlotFile();

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Return required number of ghost cells
  */
  virtual int queryGhost() const;

protected:
  
  /*!
    @brief Relaxation type for gmg
  */
  relax m_gmg_relax_type;

  /*!
    @brief GMG multigrid type
  */
  amrmg m_gmg_type;

  /*!
    @brief Needs setup
  */
  bool m_needs_setup;

  /*!
    @brief Has deeper MG levels or not
  */
  bool m_has_mg_stuff;

  /*!
    @brief Use tga for transient solevs
  */
  bool m_use_tga;

  /*!
    @brief Use the source term to draw stochastic photons
  */
  bool m_stochastic_photons;

  /*!
    @brief Verbosity for geometric multigrid
  */
  int m_gmg_verbosity;

  /*!
    @brief Precoarsening for multigrid
  */
  int m_gmg_coarsen;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_gmg_pre_smooth;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_gmg_post_smooth;

  /*!
    @brief Number of smoothing before bottom solver
  */
  int m_gmg_bot_smooth;

  /*!
    @brief Maximum number of iterations
  */
  int m_gmg_max_iter;

  /*!
    @brief Minimum number of iterations
  */
  int m_gmg_min_iter;

  /*!
    @brief Bottom solver type
  */
  int m_bottomsolver;

  /*!
    @brief Number of smoothing for bottom solver
  */
  int m_numsmooth;

  /*!
    @brief Set bottom drop depth
  */
  int m_bottom_drop;

  /*!
    @brief RNG seed
  */
  int m_seed;

  /*!
    @brief 
  */
  Real m_gmg_eps;
  
  /*!
    @brief 
  */
  Real m_gmg_hang;

  /*!
    @brief Reflection coefficient
  */
  Real m_r1;

  /*!
    @brief Reflection coefficient
  */
  Real m_r2;

  /*!
    @brief RNG
  */
  std::mt19937_64* m_rng;

  /*!
    @brief Geometric multigrid solver
  */
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > m_gmg_solver;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB> > m_bicgstab;

  /*!
    @brief TGA solver
  */
  RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > > m_tgasolver;

  /*!
    @brief Euler solver
  */
  RefCountedPtr<EBBackwardEuler> m_eulersolver;

  /*!
    @brief Operator factory
  */
  RefCountedPtr<ebconductivityopfactory> m_opfact;

  /*!
    @brief Larsen coefficients for photon outflow bc
  */
  RefCountedPtr<larsen_coefs> m_robinco;

  /*!
    @brief Photon outflow boundary conditions
  */
  RefCountedPtr<robinconductivityebbcfactory> m_ebfact;

  /*!
    @brief Photon outflow boundary conditions
  */
  RefCountedPtr<robinconductivitydomainbcfactory> m_domfact;

  /*!
    @brief Wall boundary conditions
  */
  Vector<RefCountedPtr<wall_bc> > m_wallbc;

  /*!
    @brief MG levelsgrids
  */
  Vector<EBLevelGrid> m_mg_levelgrids;

  /*!
    @brief multi-fluid simple solver
  */
  EBSimpleSolver m_simple_solver;
  
  /*!
    @brief Residue for stationary solver
  */
  EBAMRCellData m_resid;

  /*!
    @brief a-coefficient
  */
  EBAMRCellData m_aCoefficient;

  /*!
    @brief b-coefficient
  */
  EBAMRFluxData m_bco;

  /*!
    @brief b-coefficient
  */
  EBAMRIVData m_bco_irreg;

  /*!
    @brief Set up geometric multigrid
  */
  virtual void setup_gmg();

  /*!
    @brief Set kappa
  */
  virtual void set_coefficients();

  /*!
    @brief Set the a-coefficient and b-coefficient
  */
  virtual void set_aco_and_bco();

  /*!
    @brief Set both aco and bco (irregular version)
  */
  virtual void set_aco_and_bco_box(EBCellFAB&       a_aco,
				   BaseIVFAB<Real>& a_bco_irreg,
				   const Box        a_cellbox,
				   const RealVect   a_origin,
				   const Real       a_dx,
				   const int        a_lvl,
				   const DataIndex& a_dit);

  /*!
    @brief Define the GMG levels
  */
  virtual void define_mg_levels();
  
  /*!
    @brief Set the operator factory
  */
  virtual void setup_operator_factory();

  /*!
    @brief Set the GMG solver
  */
  virtual void setup_multigrid();

  /*!
    @brief Set Neumann wall boundary condition
  */
  virtual void set_neumann_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value);

  /*!
    @brief Set Robin wall BC
  */
  virtual void set_robin_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value);

  /*!
    @brief Set up the TGA solver
  */
  virtual void setup_tga();
  
  /*!
    @brief Set up the Euler solver
  */
  virtual void setup_euler();

  // Parse options functions
  virtual void parseDomainBc();    // Parses domain BC options
  virtual void parse_stationary();   // Parse stationary solver
  virtual void parsePlotVariables();    // Parses plot variables
  virtual void parse_gmg_settings(); // Parses solver parameters for geometric multigrid
  virtual void parse_reflection();   // Parse reflection coefficient
};
#include "CD_NamespaceFooter.H"
#endif
