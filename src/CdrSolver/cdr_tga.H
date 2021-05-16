/*!
  @file   cdr_tga.H
  @brief  Extensions of CdrSolver which use ebconductivityop and TGA for diffusion parts
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _CDR_TGA_
#define _CDR_TGA_

#include "CD_CdrSolver.H"
#include "ebconductivityopfactory.H"

#include <AMRMultiGrid.H>
#include <AMRTGA.H>
#include <BiCGStabSolver.H>
#include <EBBackwardEuler.H>
#include <EBSimpleSolver.H>

#include <random>
#include <time.h>
#include <chrono>

#include "CD_NamespaceHeader.H"
  
/*!
  @brief Extension class of CdrSolver that uses ebconductivityop for diffusion part. Can also solve for stochastic diffusion
*/
class cdr_tga : public CdrSolver {
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
  cdr_tga();

  /*!
    @brief Destructor
  */
  virtual ~cdr_tga();

  /*!
    @brief Parse class options
  */
  virtual void parseOptions() = 0;

  /*!
    @brief Parse class options
  */
  virtual void parseRuntimeOptions() = 0;

  /*!
    @brief Compute div(nv - D*grad(n)) explicitly
  */
  virtual void computeDivJ(EBAMRCellData& a_divJ,
			    EBAMRCellData& a_phi,
			    const Real     a_extrapDt,
			    const bool     a_ebFlux = true);

  /*!
    @brief Compute div(v*phi) explicitly
  */
  virtual void computeDivF(EBAMRCellData& a_divF,
			    EBAMRCellData& a_phi,
			    const Real     a_extrapDt,
			    const bool     a_ebFlux = true);

  /*!
    @brief Compute div(D*grad(phi)) explicitly
  */
  virtual void computeDivD(EBAMRCellData& a_divD,
			    EBAMRCellData& a_phi,
			    const bool     a_ebFlux = false);

  /*!
    @brief Set implicit diffusion solver solver parameters. These are normally good for this class. 
  */
  virtual void set_gmg_solver_parameters(relax      a_relax_type  = relax::gsrb_fast,
					 amrmg      a_gmg_type    = amrmg::vcycle,
					 const int  a_verbosity   = 10,
					 const int  a_pre_smooth  = 2,
					 const int  a_post_smooth = 2,
					 const int  a_bot_smooth  = 2,
					 const int  a_max_iter    = 16,
					 const int  a_min_iter    = 5,
					 const Real a_eps         = 1.E-8,
					 const Real a_hang        = 0.2);

  /*!
    @brief Set the bottom solver.  
    @param[in] a_whichsolver Solver type. 0 = smoother, 1 = BiCGStab
  */
  virtual void set_bottom_solver(const int a_whichsolver);

  /*!
    @brief Set the depth for bottom drop
  */
  virtual void set_bottom_drop(const int a_bottom_drop);

  /*!
    @brief Set number of smoothing for the bottom solver (simplified solver only)
  */
  virtual void set_botsolver_smooth(const int a_smooth);

  /*!
    @brief Use tga for temporal integration
  */
  virtual void set_tga(const bool a_use_tga);

  /*!
    @brief Query number of ghost cells
  */
  virtual int queryGhost() const = 0;

  /*!
    @brief Implicit diffusion TGA advance with source terms
  */
  virtual void advanceTGA(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt);

  /*!
    @brief Implicit diffusion TGA advance with source terms
  */
  virtual void advanceTGA(EBAMRCellData&       a_newPhi,
			   const EBAMRCellData& a_oldPhi,
			   const EBAMRCellData& a_source,
			   const Real           a_dt);

  /*!
    @brief Implicit diffusion Euler advance without source terms
  */
  virtual void advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt);

  /*!
    @brief Implicit diffusion Euler advance with source terms
  */
  virtual void advanceEuler(EBAMRCellData&       a_newPhi,
			     const EBAMRCellData& a_oldPhi,
			     const EBAMRCellData& a_source,
			     const Real           a_dt);

  /*!
    @brief Write plot data
  */
  virtual void writePlotData(EBAMRCellData& a_output, int& a_comp);

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
    @brief Geometric multigrid solver
  */
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > m_gmg_solver;

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
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB> > m_bicgstab;

  /*!
    @brief Simple solver
  */
  EBSimpleSolver m_simple_solver;

  /*!
    @brief MG levelsgrids
  */
  Vector<EBLevelGrid> m_mg_levelgrids;

  /*!
    @brief Verbosity for geometric multigrid
  */
  int m_gmg_verbosity;

  /*!
    @brief Pre-coarsening of GMG levels
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
    @brief Use tga
  */
  bool m_use_tga;

  /*!
    @brief Has deeper MG levels or not
  */
  bool m_has_mg_stuff;

  /*!
    @brief 
  */
  Real m_gmg_eps;
  
  /*!
    @brief 
  */
  Real m_gmg_hang;

  /*!
    @brief Advection-only extrapolation to faces
  */
  virtual void advect_to_faces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt) = 0;

  /*!
    @brief Set up geometric multigrid
  */
  virtual void setup_gmg();

  /*!
    @brief Define the GMG levels
  */
  virtual void define_mg_levels();

  /*!
    @brief Setup the operator factory
  */
  virtual void setup_operator_factory();

  /*!
    @brief Setup multigrid
  */
  virtual void setup_multigrid();

  /*!
    @brief Setup the TGA solver
  */
  virtual void setup_tga();

  /*!
    @brief Setup the Euler solver
  */
  virtual void setup_euler();

  /*!
    @brief Parse solver settings for geometric multigrid
  */
  virtual void parse_gmg_settings();
};
#include "CD_NamespaceFooter.H"
#endif
