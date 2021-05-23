/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */
/*!
  @file   CD_EbHelmholtzOp.H
  @brief  Declares Helmholtz operator which fits in Chombo3's multigrid templates. 
  @author Robert Marskar
*/

#ifndef CD_EbHelmholtzOp_H
#define CD_EbHelmholtzOp_H

// Std includes
#include <map>

// Chombo includes
#include <REAL.H>
#include <Box.H>
#include <FArrayBox.H>
#include <Vector.H>
#include <RefCountedPtr.H>
#include <AMRMultiGrid.H>
#include <EBAMRPoissonOp.H>
#include <EBIndexSpace.H>
#include <EBCellFAB.H>
#include <EBCellFactory.H>
#include <EBStencil.H>
#include <EBLevelDataOps.H>
#include <BaseEBBC.H>
#include <BaseDomainBC.H>
#include <CFIVS.H>
#include <EBFluxRegister.H>
#include <EBFastFR.H>
#include <EBMGAverage.H>
#include <EBMGInterp.H>
#include <PolyGeom.H>
#include <EBQuadCFInterp.H>
#include <EBLevelGrid.H>
#include <AMRTGA.H>
#include <AMRPoissonOp.H>
#include <CFRegion.H>
#include <ConductivityBaseDomainBC.H>

// Our includes
#include <CD_NamespaceHeader.H>

#if CH_SPACEDIM==2
#define EBAMRPOVC_NUMSTEN 4
#elif CH_SPACEDIM==3
#define EBAMRPOVC_NUMSTEN 8
#endif

/*!
  @brief Class which implements an operator that can solve the Helmholtz equation using Chombo's AMRLevelOp interface
*/
class EbHelmholtzOp : public LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB > {
public:

  /*!
    @brief Weak constructor is disallowed. 
  */
  EbHelmholtzOp() = delete;

  /*!
    @brief Copy constructor is disallowed
  */
  EbHelmholtzOp(const EbHelmholtzOp& a_other) = delete;

  /*!
    @brief Constructor for EbHelmholtzOp's. 
    @note Ghost cell arguments are there for caching reasons. Once you set them, an error is thrown if you send in data that does not match.
  */
  EbHelmholtzOp(const EBLevelGrid &                                  a_eblgFine,
		const EBLevelGrid &                                  a_eblg,
		const EBLevelGrid &                                  a_eblgCoar,
		const EBLevelGrid &                                  a_eblgCoarMG,
		const RefCountedPtr<EBQuadCFInterp>&                 a_quadCFI,
		const RefCountedPtr<EBFluxRegister>&                 a_fastFR,
		const RefCountedPtr<ConductivityBaseDomainBC>&       a_domainBC,
		const RefCountedPtr<ConductivityBaseEBBC>&           a_ebBC,
		const Real&                                          a_dx,
		const Real&                                          a_dxCoar,
		const int&                                           a_refToFine,
		const int&                                           a_refToCoar,
		const bool&                                          a_hasFine,
		const bool&                                          a_hasCoar,
		const bool&                                          a_hasMGObjects,
		const bool&                                          a_layoutChanged,
		const Real&                                          a_alpha,
		const Real&                                          a_beta,
		const RefCountedPtr<LevelData<EBCellFAB> >&          a_acoef,
		const RefCountedPtr<LevelData<EBFluxFAB> >&          a_bcoef,
		const RefCountedPtr<LevelData<BaseIVFAB<Real> > >&   a_bcoIrreg,
		const IntVect&                                       a_ghostCellsPhi,
		const IntVect&                                       a_ghostCellsRHS,
		const int&                                           a_relaxType);

  /*!
    @brief Destructor
  */
  ~EbHelmholtzOp();

  /*!
    @brief Assignment constructor is disallowed
  */
  void operator=(const EbHelmholtzOp& a_opin) = delete;

  /*!
    @brief Weights by kappa. time and tga have their demands.
    @param[in] a_rhs Data to be scaled by volume fraction. 
  */
  void kappaScale(LevelData<EBCellFAB> & a_rhs);

  /*!
    @brief Computes a_residual = a_rhs - L(a_phiFine, a_phi)
    @details Assumes no coarser AMR level. 
    @param[out] a_residual      Residual
    @param[in]  a_phiFine       Fine-level solution
    @param[in]  a_phi           Current-level solution. 
    @param[in]  a_rhs           Right-hand side
    @param[in]  a_homogeneousBc Use homogeneous BCs or not. 
    @param[in]  a_finerOp       Fine-grid operator
    @details Assumes no coaser AMR level
  */
  void AMRResidualNC(LevelData<EBCellFAB>&              a_residual,
		     const LevelData<EBCellFAB>&        a_phiFine,
		     const LevelData<EBCellFAB>&        a_phi,
		     const LevelData<EBCellFAB>&        a_rhs,
		     bool                               a_homogeneousBC,
		     AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp) override;


  /*!
    @brief Applies the AMR operator with no coaser AMR level
    @details Assumes no coarser AMR level. 
    @param[out] a_LofPhi        Operator applied to phi
    @param[in]  a_phiFine       Fine-level solution
    @param[in]  a_phi           Current-level solution. 
    @param[in]  a_homogeneousBc Use homogeneous BCs or not. 
    @param[in]  a_finerOp       Fine-grid operator
  */
  void AMROperatorNC(LevelData<EBCellFAB>&       a_LofPhi,
		     const LevelData<EBCellFAB>& a_phiFine,
		     const LevelData<EBCellFAB>& a_phi,
		     bool a_homogeneousBC,
		     AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp) override;


  /*!
    @brief Set alpha and beta coefficients. Need when e.g. using backward euler discretization. 
    @param[in] a_alpha alpha-coefficient
    @param[in] a_beta  beta-coefficient
    @details This also recomputes relaxation weights.
  */
  void setAlphaAndBeta(const Real& a_alpha, const Real& a_beta) override;

  /*!
    @brief Do diagonal scaling
    @param[inout] a_rhs Data to be scaled diagonally
    @param[in]    a_kappaWeighted If true, also multiply by volume fraction.
  */
  void diagonalScale(LevelData<EBCellFAB> & a_rhs, bool a_kappaWeighted) override;

  /*!
    @brief Divide by m_acoefficient
    @param[in] a_rhs Data to be divided.
  */
  void divideByIdentityCoef(LevelData<EBCellFAB> & a_rhs) override;

  /*!
    @brief Part of the LevelTGA interface
  */
  void fillGrad(const LevelData<EBCellFAB>& a_phi) override;

  /*!
    @brief Compute flux on face centroids
    @param[out] a_flux  Flux
    @param[in]  a_data  Data
    @param[in]  a_grid  Box 
    @param[in]  a_dit   DataIndex
    @param[in]  a_scale Flux scaling
    @details This iterates over the sides and calls the other (directional) version. 
  */
  void getFlux(EBFluxFAB&                    a_flux,
	       const LevelData<EBCellFAB>&   a_data,
	       const Box&                    a_grid,
	       const DataIndex&              a_dit,
	       Real                          a_scale);

  /*!
    @brief Compute flux on face centroids
    @param[out] a_fluxCentroid Flux ont he centroids
    @param[in]  a_phi          Data
    @param[in]  a_ghostedBox   Ghosted box
    @param[in]  a_fabBox       Grid box
    @param[in]  a_domain       Problem domain
    @param[in]  a_ebisBox      EBIS box
    @param[in]  a_dx           Current grid resolution.
    @param[in]  a_datInd       Current dataindex
    @param[in]  a_idir         Flux "direction"
    @details This computes face-centers first, then interpolates to centroids. 
  */
  void getFlux(EBFaceFAB&                    a_fluxCentroid,
	       const EBCellFAB&              a_phi,
	       const Box&                    a_ghostedBox,
	       const Box&                    a_fabBox,
	       const ProblemDomain&          a_domain,
	       const EBISBox&                a_ebisBox,
	       const Real&                   a_dx,
	       const DataIndex&              a_datInd,
	       const int&                    a_idir);


  /*!
    @brief Set time for this operator
    @param[in] a_oldTime Old time
    @param[in] a_mu      Mu
    @param[in] a_dt      Time step. 
    @note  This does nothing as time is disabled for EbHelmholtzOp.
  */
  void setTime(Real a_oldTime, Real a_mu, Real a_dt) override;

  /*!
    @brief Compute the residual
    @param[out] a_residual          Residual
    @param[in]  a_phi               Cell-centered data
    @param[in]  a_rhs               Right-hand side
    @param[in]  a_homogeneousPhysBC Use homogeneous BC or not
  */
  void residual(LevelData<EBCellFAB>&       a_residual,
		const LevelData<EBCellFAB>& a_phi,
		const LevelData<EBCellFAB>& a_rhs,
		bool                        a_homogeneousPhysBC=false) override;

  /*!
    @brief Precondition before bottom solver.
    @param[out] a_opPhi Preconditioned soln
    @param[out] a_phi   Data
    @details This just does a few relaxations.
  */
  void preCond(LevelData<EBCellFAB>& a_opPhi, const LevelData<EBCellFAB>& a_phi) override;

  /*!
    @brief Apply operator to solution
    @param[out] a_opPhi             Operator applied to soluion. 
    @param[in]  a_phi               Cell-centered data on this level
    @param[in]  a_phiCoarse         Coarser data
    @param[in]  a_homogeneousPhysBC Homogeneous physical BCs or not
    @param[in]  a_homogeneousCFBC   Homogeneous physical coarse-fine boundary conditions or not
    @param[in]  a_doExchange        Do an exchange operation or not. 
    @details This function assumes that coarse-fine boundary condtions have been dealt with.
  */
  void applyOp(LevelData<EBCellFAB>&             a_opPhi,
	       const LevelData<EBCellFAB>&       a_phi,
	       const LevelData<EBCellFAB>* const a_phiCoarse,
	       const bool&                       a_homogeneousPhysBC,
	       const bool&                       a_homogeneousCFBC,
	       const bool&                       a_doExchange = true);

  /*!
    @brief Apply operator to solution
    @param[out] a_opPhi             Operator applied to soluion. 
    @param[in]  a_phi               Cell-centered data on this level
    @param[in]  a_phiCoarse         Coarser data
    @param[in]  a_homogeneousPhysBC Homogeneous physical BCs or not
    @param[in]  a_homogeneousCFBC   Homogeneous physical coarse-fine boundary conditions or not
    @param[in]  a_dit               Data iterator
    @param[in]  a_doExchange        Do an exchange operation or not. 
    @details This function assumes that coarse-fine boundary condtions have been dealt with.
  */
  void applyOp(LevelData<EBCellFAB>&             a_opPhi,
	       const LevelData<EBCellFAB>&       a_phi,
	       const LevelData<EBCellFAB>* const a_phiCoarse,
	       const bool&                       a_homogeneousPhysBC,
	       const bool&                       a_homogeneousCFBC,
	       DataIterator&                     a_dit,
	       const bool&                       a_doExchange = true);

  /*!
    @brief Applies operator but without boundary conditions
    @param[in] a_opPhi L(phi)
    @param[in] a_phi   Cell-centered data
    @details Function called by LevelTGA
  */
  void applyOpNoBoundary(LevelData<EBCellFAB>& a_opPhi, const LevelData<EBCellFAB>& a_phi) override;

  /*!
    @brief Apply operator
    @param[out] a_opPhi Operator
    @param[in]  a_phi   Data
    @param[in]  a_homogeneousPhysBC Homogeneous physical BCs or not. 
    @details Coarse is null because that's all we can do. 
  */
  void applyOp(LevelData<EBCellFAB>& a_opPhi, const LevelData<EBCellFAB>& a_phi, bool a_homogeneousPhysBC) override;

  /*!
    @brief Create a clone of data (does not copy)
    @param[out] a_lhs Clone
    @param[out] a_rhs Stuff to clone
  */
  void create(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override;

  /*!
    @brief Create a coarsened version of data
    @param[out] a_lhs Coarsened data structure
    @param[in]  a_rhs Data tructure
    @param[in]  a_refRat Refinement ratio.
  */
  void createCoarsened(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int& a_refRat) override;

  /*!
    @brief Set lhs = rhs
    @param[out] a_lhs Left hand side
    @param[in]  a_rhs Data
  */
  void assign(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override;

  /*!
    @brief Compute dot product
    @param[in] a_1 Data
    @param[in] a_2 Data
    @return Dot product, includes kappa. 
  */
  Real dotProduct(const LevelData<EBCellFAB>& a_1, const LevelData<EBCellFAB>& a_2) override;

  /*!
    @brief Increment data
    @param[inout] a_lhs   Data to be incremented
    @param[inout] a_x     Incrementation
    @param[inout] a_scale Incrementation scale
  */
  void incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_x, Real a_scale) override;

  /*!
    @brief Compute a*x + b*y
    @param[out] a_lhs Left hand side
    @param[in]  a_x   Data
    @param[in]  a_y   Data
    @param[in]  a_a   Scale
    @param[in]  a_b   Scale
  */
  void axby(LevelData<EBCellFAB>&       a_lhs,
	    const LevelData<EBCellFAB>& a_x,
	    const LevelData<EBCellFAB>& a_y,
	    Real                        a_a,
	    Real                        a_b) override;

  /*!
    @brief Scale data
    @param[inout] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  void scale(LevelData<EBCellFAB>& a_lhs, const Real& a_scale) override;

  /*!
    @brief Compute normal of data
    @param[in] a_rhs Data 
    @param[in] a_ord Norm order
  */
  Real norm(const LevelData<EBCellFAB>& a_rhs, int a_ord) override;

  /*!
    @brief Utility function which computes norm on local data
    @param[in] a_rhs Data
  */
  Real localMaxNorm(const LevelData<EBCellFAB>& a_rhs) override;

  /*!
    @brief Set data to zero
    @param[inout] a_lhs Data
  */
  void setToZero(LevelData<EBCellFAB>& a_lhs) override;

  /*!
    @brief Set constant value in data
    @param[inout] a_lhs  Data
    @param[in]    a_value Value
  */
  void setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value); 

  /*!
    @brief Create coarser data from fine data
    @param[out] a_coarse  Coarse leveldata
    @param[in]  a_fine    Fine leveldata
    @param[in]  a_ghosted Dummy argument. Does nothing
  */
  void createCoarser(LevelData<EBCellFAB>&       a_coarse,
		     const LevelData<EBCellFAB>& a_fine,
		     bool                        a_ghosted) override;


  /*!
    @brief Relax data
    @param[inout] a_e          Err
    @param[in]    a_residual   Resid
    @param[in]    a_iterations Number of iterations
  */
  void relax(LevelData<EBCellFAB>&       a_e,
	     const LevelData<EBCellFAB>& a_residual,
	     int                         a_iterations);


  /*!
    @brief Relax data using gauss-seidel
    @param[inout] a_e          Err
    @param[in]    a_residual   Resid
    @param[in]    a_iterations Number of iterations
  */
  void relaxGauSai(LevelData<EBCellFAB>&       a_e,
		   const LevelData<EBCellFAB>& a_residual,
		   int                         a_iterations);
  /*!
    @brief Relax data using Jacobi iteration
    @param[inout] a_e          Err
    @param[in]    a_residual   Resid
    @param[in]    a_iterations Number of iterations
  */
  void relaxPoiJac(LevelData<EBCellFAB>&       a_e,
		   const LevelData<EBCellFAB>& a_residual,
		   int                         a_iterations);

  /*!
    @brief Relax data using red-black Gauss-Seidel
    @param[inout] a_e          Err
    @param[in]    a_residual   Resid
    @param[in]    a_iterations Number of iterations
  */
  void relaxGSRBFast(LevelData<EBCellFAB>&       a_e,
		     const LevelData<EBCellFAB>& a_residual,
		     int                         a_iterations);

  /*!
    @brief Restrict the residual
    @param[out] a_resCoarse Residual on coarse level
    @param[inout] a_phiFine Fine data
    @param[inout] a_rhsFine Fine rhs
    @details Calculate restricted residual: a_resCoarse[2h] = I[h->2h] (a_rhsFine[h] - L[h](a_phiFine[h]))
  */
  void restrictResidual(LevelData<EBCellFAB>& a_resCoarse, LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_rhsFine) override;

  /*!
    @brief Prolongation function
    @param[inout] a_phiThisLevel Correction on this level
    @param[in]    a_correctCoarse Coarse-level correction
    @details Correct the fine solution based on coarse correction: a_phiThisLevel += I[2h->h] (a_correctCoarse)
  */
  void prolongIncrement(LevelData<EBCellFAB>& a_phiThisLevel, const LevelData<EBCellFAB>& a_correctCoarse) override;

  /*!
    @brief Refinement ratio between this level and coarser level.
    @details Returns 1 when there are no coarser AMRLevelOp objects 
  */
  int refToCoarser();

  /*!
    @brief Refinement ratio between this level and coarser level.
    @details Returns 1 when there are no coarser AMRLevelOp objects
  */
  int refToFiner();

  /*!
    @brief Compute the residual
    @param[out] a_residual     Residual
    @param[in]  a_phiFine       Fine data
    @param[in]  a_phi           This data
    @param[in]  a_phiCoar       Coarse data
    @param[in]  a_rhs           right-hand side
    @param[in]  a_homogeneousBC Use homogeneous BCs or not
    @param[in]  a_finerOp       Finer operator (if it exists
    @details Computes a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) 
  */
  void AMRResidual(LevelData<EBCellFAB>&              a_residual,
		   const LevelData<EBCellFAB>&        a_phiFine,
		   const LevelData<EBCellFAB>&        a_phi,
		   const LevelData<EBCellFAB>&        a_phiCoarse,
		   const LevelData<EBCellFAB>&        a_rhs,
		   bool                               a_homogeneousBC,
		   AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp) override;

  /*!
    @brief Compute residual on finest level
    @param[out] a_residual      Residual
    @param[in]  a_phi           Data
    @param[in]  a_phiCoarse     Coarse data
    @param[in]  a_rhs           Right-hand side
    @param[in]  a_homogeneousBC Homogeneous BC or not
    
    @details Computes a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
  void AMRResidualNF(LevelData<EBCellFAB>& a_residual,
		     const LevelData<EBCellFAB>& a_phi,
		     const LevelData<EBCellFAB>& a_phiCoarse,
		     const LevelData<EBCellFAB>& a_rhs,
		     bool a_homogeneousBC) override;

  /*!
    @brief AMR operator
    @param[out] a_LofPhi        L(phi)
    @param[in]  a_phiFine        Fine phi
    @param[in]  a_phi           Phi
    @param[in]  a_phiCoarse     Coarse phi
    @param[in]  a_homogeneousBC Homogeneous BC or not
  */
  void AMROperator(LevelData<EBCellFAB>& a_LofPhi,
		   const LevelData<EBCellFAB>& a_phiFine,
		   const LevelData<EBCellFAB>& a_phi,
		   const LevelData<EBCellFAB>& a_phiCoarse,
		   bool a_homogeneousBC,
		   AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp) override;

  /*!
    @brief Apply linear operator assuming no finer level
    @param[out] a_LofPhi        L(phi)
    @param[in]  a_phi           Phi
    @param[in]  a_phiCoarse     Coarse phi
    @param[in]  a_homogeneousBC Homogeneous BC or not
  */
  void AMROperatorNF(LevelData<EBCellFAB>& a_LofPhi,
		     const LevelData<EBCellFAB>& a_phi,
		     const LevelData<EBCellFAB>& a_phiCoarse,
		     bool a_homogeneousBC) override;


  /*!
    @brief Restrict residual to coarse
    @param[out] a_resCoars         Coarse residual
    @param[out] a_residual         Residual
    @param[out] a_correction       Correction
    @param[out] a_coarseCorrection Coarse correction
    @details Computes resCoarse = I[h-2h] (a_residual - L(a_correction, a_coarseCorrection)) 
  */
  void AMRRestrict(LevelData<EBCellFAB>& a_resCoarse,
		   const LevelData<EBCellFAB>& a_residual,
		   const LevelData<EBCellFAB>& a_correction,
		   const LevelData<EBCellFAB>& a_coarseCorrection, 
		   bool a_skip_res = false ) override;


  /*!
    @brief Prolong ocrection
    @param[inout] a_correction       Correction
    @param[out]   a_coarseCorrection Coarse correction
    @details Compute a_correction += I[2h->h](a_coarseCorrection) 
  */
  void AMRProlong(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_coarseCorrection) override;

  /*!
    @brief Update AMR Residual
    @param[in]  a_residual         Residual
    @param[out] a_correction       Correction
    @param[out] a_coarseCorrection Coarse correction
    @details a_residual = a_residual - L(a_correction, a_coarseCorrection) 
  */
  void AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
				 const LevelData<EBCellFAB>& a_correction,
				 const LevelData<EBCellFAB>& a_coarseCorrection) override;

  /*!
    @brief Reflux solution
    @param[inout] a_residual Residual
    @param[in]    a_phiFine  Fine phi
    @param[in]    a_phi      Phi
    @param[in]    a_finerOp  Finer operator
  */
  void reflux(LevelData<EBCellFAB>& a_residual,
	      const LevelData<EBCellFAB>& a_phiFine,
	      const LevelData<EBCellFAB>& a_phi,
	      AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  /*!
    @brief Do red-black color magic
  */
  void gsrbColor(LevelData<EBCellFAB>&       a_phi,
		 const LevelData<EBCellFAB>& a_lph,
		 const LevelData<EBCellFAB>& a_rhs,
		 const IntVect&              a_color);

  /*!
    @brief Do red-black color magic
  */
  void gsrbColor(EBCellFAB&       a_phi,
		 const EBCellFAB& a_lph,
		 const EBCellFAB& a_rhs,
		 const Box&       a_box,
		 const DataIndex& a_dit,
		 const IntVect&   a_color);

  /*!
    @brief Get divergence stencil in point
  */
  void getDivFStencil(VoFStencil&      a_vofStencil,
		      const VolIndex&  a_vof,
		      const DataIndex& a_dit);

  /*!
    @brief Get flux stencil for face
  */
  void getFluxStencil(VoFStencil&      a_fluxStencil,
		      const FaceIndex& a_face,
		      const DataIndex& a_dit);

  /*!
    @brief Get flux stencil for face
  */
  void getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
				  const FaceIndex& a_face,
				  const DataIndex& a_dit);

  /*!
    @brief Increment regular operator
  */
  void incrOpRegularDir(EBCellFAB&             a_lhs,
			const EBCellFAB&       a_phi,
			const bool&            a_homogeneous,
			const int&             a_dir,
			const DataIndex&       a_datInd);

  /*!
    @brief Apply operator in irregular cells
  */
  void applyOpIrregular(EBCellFAB&             a_lhs,
			const EBCellFAB&       a_phi,
			const bool&            a_homogeneous,
			const DataIndex&       a_datInd);

  /*!
    @brief Get EBLevelGrid
  */
  const EBLevelGrid& getEBLG() const {
    return m_eblg;
  }
  
protected:
  
  /*!
    @brief For turning on/off BCs
  */
  static bool s_turnOffBCs;

  /*!
    @brief Get safety factor
    @details Returns 1.0
  */
  Real getSafety() const;

  /*!
    @brief Define all stencils
  */
  void defineStencils();

  /*!
    @brief For computing regular amrop
  */
  void incrOpRegularAllDirs(Box * a_loBox,
			    Box * a_hiBox,
			    int * a_hasLo,
			    int * a_hasHi,
			    Box & a_curDblBox,
			    Box & a_curPhiBox,
			    int a_nComps,
			    BaseFab<Real> & a_curOpPhiFAB,
			    const BaseFab<Real> & a_curPhiFAB,
			    bool a_homogeneousPhysBC,
			    const DataIndex& a_dit);
  /*!
    @brief Apply domain flux
  */
  void applyDomainFlux(Box * a_loBox,
		       Box * a_hiBox,
		       int * a_hasLo,
		       int * a_hasHi,
		       Box & a_dblBox,
		       int a_nComps,
		       BaseFab<Real> & a_phiFAB,
		       bool a_homogeneousPhysBC,
		       const DataIndex& a_dit);

  /*!
    @brief Gauss-seidel in irregular cells
  */
  void GSColorAllIrregular(EBCellFAB&                        a_phi,
			   const EBCellFAB&                  a_rhs,
			   const int&                        a_icolor,
			   const DataIndex&                  a_dit);


  /*!
    @brief Calculate diagonal weight
  */
  void calculateAlphaWeight();

  /*!
    @brief Calcualte relaxation coefficient
  */
  void calculateRelaxationCoefficient();

  /*!
    @brief Define stencils
  */
  void defineColorStencils(Box a_ideBoxLo[SpaceDim], Box a_ideBoxHi[SpaceDim]);
  
  /*!
    @brief Stencils for doing EBCF stuff
  */
  void defineEBCFStencils();

  /*!
    @brief Get flux, knowing that there is an EBCF.
  */
  void getFluxEBCF(EBFaceFAB&                    a_flux,
		   const EBCellFAB&              a_phi,
		   const Box&                    a_ghostedBox,
		   Vector<FaceIndex>&            a_faceitEBCF,
		   Vector<VoFStencil>&           a_stenEBCF);

  /*!
    @brief Get flux
  */
  void getFluxRegOnly(EBFaceFAB&                    a_fluxCentroid,
		      const EBCellFAB&              a_phi,
		      const Box&                    a_ghostedBox,
		      const Real&                   a_dx,
		      const DataIndex&              a_datInd,
		      const int&                    a_idir);

  // Stuff to make EBCF go faster
  LayoutData< Vector<FaceIndex>  >   m_faceitCoar[2*SpaceDim];
  LayoutData< Vector<VoFStencil> >  m_stencilCoar[2*SpaceDim];

  //stuff to make relaxation go faster
  //not for the faint of heart
  LayoutData<VoFIterator >         m_vofItIrregColorDomLo[EBAMRPOVC_NUMSTEN][SpaceDim];
  LayoutData<VoFIterator >         m_vofItIrregColorDomHi[EBAMRPOVC_NUMSTEN][SpaceDim];
  LayoutData<BaseIVFAB<Real> >     m_cacheEBxDomainFluxLo[EBAMRPOVC_NUMSTEN][SpaceDim];
  LayoutData<BaseIVFAB<Real> >     m_cacheEBxDomainFluxHi[EBAMRPOVC_NUMSTEN][SpaceDim];
  LayoutData<RefCountedPtr<EBStencil> >  m_colorEBStencil[EBAMRPOVC_NUMSTEN];

  int                             m_relaxType;
  const IntVect                   m_ghostCellsPhi;
  const IntVect                   m_ghostCellsRHS;

  RefCountedPtr<EBQuadCFInterp>   m_quadCFIWithCoar;
  RefCountedPtr<EBFluxRegister>   m_ext_fastFR;

  EBLevelGrid                     m_eblg;
  EBLevelGrid                     m_eblgFine;
  EBLevelGrid                     m_eblgCoar;
  EBLevelGrid                     m_eblgCoarMG;
  EBLevelGrid                     m_eblgCoarsenedFine;

  RefCountedPtr<ConductivityBaseDomainBC>     m_domainBC;
  RefCountedPtr<ConductivityBaseEBBC>         m_ebBC;

  Real                            m_dxFine;
  Real                            m_dx;
  Real                            m_dxCoar;

  //! "Current" (time-interpolated) value of the a coefficient. For a
  //! time-independent a coefficient, this is where the coefficient lives.
  RefCountedPtr<LevelData<EBCellFAB> >          m_aCoef;

  //! "Current" (time-interpolated) value of the b coefficient. For a
  //! time-independent a coefficient, this is where the coefficient lives.
  RefCountedPtr<LevelData<EBFluxFAB> >          m_bcoef;

  //! "Current" (time-interpolated) value of the b coefficient on irregular
  //! cells.
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >   m_bcoIrreg;

  Real                            m_alpha;
  Real                            m_beta;
  //weights that get multiplied by alpha
  LayoutData<BaseIVFAB<Real> >       m_alphaDiagWeight;
  //weights that get multiplied by beta
  LayoutData<BaseIVFAB<Real> >       m_betaDiagWeight;
  int                             m_refToFine;
  int                             m_refToCoar;
  bool                            m_hasEBCF;
  bool                            m_hasFine;
  bool                            m_hasInterpAve;
  bool                            m_hasCoar;

  //restriction object
  EBMGAverage                    m_ebAverage;
  //prolongation object
  EBMGInterp                     m_ebInterp;

  //stencils for operator evaluation
  LayoutData<RefCountedPtr<EBStencil> >  m_opEBStencil;
  //stencils for operator evaluation on gauss-seidel colors

  //! Multigrid relaxation coefficient
  LevelData<EBCellFAB>       m_relCoef;

  //cache the vofiterators
  //for irregular cell iteration (includes buffer around multivalued cells)
  LayoutData<VoFIterator >                     m_vofIterIrreg;
  LayoutData<VoFIterator >                     m_vofIterMulti;
  //for domain boundary conditions at ir regular cells
  LayoutData<VoFIterator >                     m_vofIterDomLo[CH_SPACEDIM];
  LayoutData<VoFIterator >                     m_vofIterDomHi[CH_SPACEDIM];

  // Coarse-fine stencils for homogeneous CFInterp
  LayoutData<CFIVS> m_loCFIVS[SpaceDim];
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

  //flux register with finer level
  RefCountedPtr<EBFluxRegister>       m_fastFR;
  //special mg objects for when we do not have
  //a coarser level or when the refinement to coarse
  //is greater than two
  //flag for when we need special MG objects
  bool                        m_hasMGObjects;
  bool                        m_layoutChanged;
  //stuff below is only defined if m_hasMGObjects==true
  EBMGAverage                 m_ebAverageMG;
  EBMGInterp                  m_ebInterpMG;
  DisjointBoxLayout           m_dblCoarMG;
  EBISLayout                  m_ebislCoarMG;
  ProblemDomain               m_domainCoarMG;

  Vector<IntVect> m_colors;

private:

  /*!
    @brief Increment flux register
  */
  void incrementFRCoar(EBFluxRegister&             a_fluxReg,
		       const LevelData<EBCellFAB>& a_phiFine,
		       const LevelData<EBCellFAB>& a_phi);

  /*!
    @brief Increment flux register
  */
  void incrementFRFine(EBFluxRegister&             a_fluxReg,
		       const LevelData<EBCellFAB>& a_phiFine,
		       const LevelData<EBCellFAB>& a_phi,
		       AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  /*!
    @brief Compute flux
  */
  void getFlux(FArrayBox&                    a_flux,
	       const FArrayBox&              a_phi,
	       const Box&                    a_faceBox,
	       const int&                    a_idir,
	       const Real&                   a_dx,
	       const DataIndex&              a_datInd);


  /*!
    @brief Interpolate ghost cells
  */
  void applyCFBCs(LevelData<EBCellFAB>&             a_phi,
		  const LevelData<EBCellFAB>* const a_phiCoarse,
		  bool a_homogeneousCFBC);

  /*!
    @brief Homogeneous CFBCs
  */
  void applyHomogeneousCFBCs(LevelData<EBCellFAB>&   a_phi);

  /*!
    @brief Homogeneous CFBCs
  */
  void applyHomogeneousCFBCs(EBCellFAB&            a_phi,
			     const DataIndex&      a_datInd,
			     int                   a_idir,
			     Side::LoHiSide        a_hiorlo);
};

#include <CD_NamespaceFooter.H>

#endif
