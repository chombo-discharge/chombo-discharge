/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrTGA.H
  @brief  Extensions of CdrSolver which use EbHelmholtzOp and TGA for diffusion parts
  @author Robert Marskar
*/

#ifndef CD_CdrTGA_H
#define CD_CdrTGA_H

// Std includes
#include <random>
#include <time.h>
#include <chrono>

// Chombo includes
#include <AMRMultiGrid.H>
#include <AMRTGA.H>
#include <BiCGStabSolver.H>
#include <EBBackwardEuler.H>
#include <EBSimpleSolver.H>

// Our includes
#include <CD_EbHelmholtzOpFactory.H>
#include <CD_CdrSolver.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Extension class of CdrSolver that uses multigrid for diffusion part. Can also solve for stochastic diffusion.
*/
class CdrTGA : public CdrSolver {
protected:

  /*!
    @brief Enum class for supported bottom solvers in multigrid.
  */
  enum class BottomSolver {
    Simple,
    BiCGStab,
    GMRES
  };

  /*!
    @brief Enum for supported relaxation methods in multigrid. 
  */
  enum class RelaxationMethod {
    Jacobi,
    GaussSeidel,
    GSRB,
  };

  /*!
    @brief Enum for multigrid cycle types. 
  */
  enum class MultigridType {
    FAS,
    VCycle,
    FCycle,
  };
  
public:

  /*!
    @brief Constructor
  */
  CdrTGA();

  /*!
    @brief Destructor
  */
  virtual ~CdrTGA();

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
			   const bool     a_ebFlux = true) override;

  /*!
    @brief Compute div(v*phi) explicitly
  */
  virtual void computeDivF(EBAMRCellData& a_divF,
			   EBAMRCellData& a_phi,
			   const CdrBc    a_domainBc,
			   const Real     a_extrapDt,
			   const bool     a_ebFlux = true) override;

  /*!
    @brief Compute div(D*grad(phi)) explicitly
    @param[out] a_divF      Divergence term, including redistribution magic. 
    @param[in]  a_phi       Cell-centered state
    @param[in]  a_domainBc  How to set domain fluxes
    @param[in]  a_useEbFlux If true, the embedded boundary flux will be included in div(D)
    @note a_phi is non-const because ghost cells will be interpolated in this routine. Valid data in a_phi is not touched. 
  */
  virtual void computeDivD(EBAMRCellData& a_divD,
			   EBAMRCellData& a_phi,
			   const CdrBc    a_domainBc,			   
			   const bool     a_ebFlux = false) override;

  /*!
    @brief Query number of ghost cells
  */
  virtual int queryGhost() const = 0;

  /*!
    @brief Implicit diffusion TGA advance without source term.
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t
    @param[in]    a_dt     Time step
    @note This version sets the source term to zero. 
  */
  virtual void advanceTGA(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt) override;

  /*!
    @brief Implicit diffusion TGA advance with source term.
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t
    @param[in]    a_source Source (at t+dt/2)
    @param[in]    a_dt     Time step
    @note Source term should be centered at t+dt/2 for second order convergence. 
  */
  virtual void advanceTGA(EBAMRCellData&       a_newPhi,
			  const EBAMRCellData& a_oldPhi,
			  const EBAMRCellData& a_source,
			  const Real           a_dt) override;

  /*!
    @brief Implicit diffusion Euler advance without source terms
  */
  virtual void advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt) override;

  /*!
    @brief Implicit diffusion Euler advance with source terms
  */
  virtual void advanceEuler(EBAMRCellData&       a_newPhi,
			    const EBAMRCellData& a_oldPhi,
			    const EBAMRCellData& a_source,
			    const Real           a_dt) override;

  /*!
    @brief Write plot data
  */
  virtual void writePlotData(EBAMRCellData& a_output, int& a_comp);

protected:

  /*!
    @brief Relaxation type for gmg
  */
  RelaxationMethod m_multigridRelaxMethod;

  /*!
    @brief GMG multigrid type
  */
  MultigridType m_multigridType;
  
  /*!
    @brief Geometric multigrid solver
  */
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > m_multigridSolver;

  /*!
    @brief TGA solver
  */
  RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > > m_tgaSolver;

  /*!
    @brief Euler solver
  */
  RefCountedPtr<EBBackwardEuler> m_eulerSolver;

  /*!
    @brief Operator factory
  */
  RefCountedPtr<EbHelmholtzOpFactory> m_operatorFactory;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB> > m_bicgstab;

  /*!
    @brief Simple solver
  */
  EBSimpleSolver m_simpleSolver;

  /*!
    @brief MG levelsgrids
  */
  Vector<EBLevelGrid> m_mg_levelgrids;

  /*!
    @brief Verbosity for geometric multigrid
  */
  int m_multigridVerbosity;

  /*!
    @brief Pre-coarsening of GMG levels
  */
  int m_numCoarseningsBeforeAggregation;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_multigridPreSmooth;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_multigridPostSmooth;

  /*!
    @brief Number of smoothing before bottom solver
  */
  int m_multigridBottomSmooth;

  /*!
    @brief Maximum number of iterations
  */
  int m_multigridMaxIterations;

  /*!
    @brief Minimum number of iterations
  */
  int m_multigridMinIterations;

  /*!
    @brief Bottom solver type
  */
  BottomSolver m_bottomSolver;

  /*!
    @brief Number of smoothing for bottom solver
  */
  int m_numSmoothingsForSimpleSolver;

  /*!
    @brief Set bottom drop depth
  */
  int m_numCellsBottomDrop;

  /*!
    @brief Has deeper MG levels or not
  */
  bool m_hasDeeperMultigridLevels;

  /*!
    @brief Multigrid tolerance
  */
  Real m_multigridTolerance;
  
  /*!
    @brief Multigrid hand
  */
  Real m_multigridHang;

  /*!
    @brief Advection-only extrapolation to faces
  */
  virtual void advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt) = 0;

  /*!
    @brief Set up geometric multigrid
  */
  virtual void setupMultigrid();

  /*!
    @brief Define the GMG levels
  */
  virtual void defineDeeperMultigridLevels();

  /*!
    @brief Setup the operator factory
  */
  virtual void setupOperatorFactory();

  /*!
    @brief Setup multigrid
  */
  virtual void setupMultigridSolver();

  /*!
    @brief Setup the TGA solver
  */
  virtual void setupTGA();

  /*!
    @brief Setup the Euler solver
  */
  virtual void setupEuler();

  /*!
    @brief Parse solver settings for geometric multigrid
  */
  virtual void parseMultigridSettings();
};

#include <CD_NamespaceFooter.H>

#endif
