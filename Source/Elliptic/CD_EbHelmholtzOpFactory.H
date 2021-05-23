/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbHelmholtzOpFactory.H
  @brief  Declaration of factory class for EbHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EbHelmholtzOpFactory_H
#define CD_EbHelmholtzOpFactory_H

// Std includes
#include <map>

// Chombo includes
#include <REAL.H>
#include <Box.H>
#include <FArrayBox.H>
#include <Vector.H>
#include <RefCountedPtr.H>
#include <AMRMultiGrid.H>
#include <EBIndexSpace.H>
#include <EBCellFAB.H>
#include <EBCellFactory.H>
#include <EBLevelDataOps.H>
#include <BaseEBBC.H>
#include <BaseDomainBC.H>
#include <CFIVS.H>
#include <EBFluxRegister.H>
#include <EBMGAverage.H>
#include <EBMGInterp.H>
#include <EBCoarsen.H>
#include <PolyGeom.H>
#include <EBAMRPoissonOp.H>
#include <EBLevelGrid.H>

// Our includes
#include <CD_EbHelmholtzOp.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Factory class for making EbHelmholtzOps. Follows Chombo interface. 
*/
class EbHelmholtzOpFactory: public AMRLevelOpFactory<LevelData<EBCellFAB> >{
public:

  /*!
    @brief Disallowed constructor
  */
  EbHelmholtzOpFactory() = delete;

  /*!
    @brief Disallowed constructor
  */
  EbHelmholtzOpFactory(const EbHelmholtzOpFactory& a_opin) = delete;

  /*!
    @brief Constructs a factory that builds EbHelmholtzOps with time-independent coefficients
  */
  EbHelmholtzOpFactory(const Vector<EBLevelGrid>&                                  a_eblgs,
		       const Vector<RefCountedPtr<EBQuadCFInterp> >&               a_quadCFI,
		       const Vector<RefCountedPtr<EBFluxRegister> >&               a_fastFR,
		       const Real&                                                 a_alpha,
		       const Real                                         &        a_beta,
		       const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
		       const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_bcoef,
		       const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_bcoefIrreg,
		       const Real&                                                 a_dxCoarse,
		       const Vector<int>&                                          a_refRatio,
		       const RefCountedPtr<BaseDomainBCFactory>&                   a_domainBCFactory,
		       const RefCountedPtr<BaseEBBCFactory>    &                   a_ebBcFactory,
		       const IntVect&                                              a_ghostCellsPhi,
		       const IntVect&                                              a_ghostCellsRhs,
		       const int&                                                  a_relaxType,
		       const int&                                                  a_botdrop = 8,
		       const int                                                   a_numLevels = -1,
		       const Vector<EBLevelGrid>&                                  m_mg_eblgs = Vector<EBLevelGrid>(0));


  /*!
    @brief Destructor
  */
  virtual ~EbHelmholtzOpFactory();

  /*!
    @brief Disallowed assignment operator
  */
  void operator=(const EbHelmholtzOpFactory& a_opin) = delete;

  /*!
    @brief Create multigrid operator
  */
  EbHelmholtzOp* MGnewOp(const ProblemDomain& a_FineindexSpace,
			 int                  a_depth,
			 bool                 a_homoOnly = true) override;

  /*!
    @brief Create AMR operator
  */
  EbHelmholtzOp* createOperator(const EBLevelGrid&             a_eblgMGLevel,
				const EBLevelGrid&             a_eblgCoarMG,
				const bool&                    a_hasMGObjects,
				const RealVect&                a_dxMGLevel,
				const RealVect&                a_dxCoar,
				const int&                     a_whichLevel);
  /*!
    @brief Reclaim (delete) operator
  */
  void reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim);

  /*!
    @brief Create operator for domain
  */
  EbHelmholtzOp* AMRnewOp(const ProblemDomain& a_domain) override;

  /*!
    @brief Delete operator
  */
  void AMRreclaim(EbHelmholtzOp* a_reclaim);

  /*! 
    @brief Get refinement ratio  
    @details Returns 1 when there are no coarser AMRLevelOp objects
  */
  int refToFiner(const ProblemDomain& a_domain) const;

  /*!
    @brief Set maximum box size for aggregration
    @param[in] a_maxBoxSize Maximum box size 
  */
  static void setMaxBoxSize(int a_maxBoxSize){
    s_maxBoxSize = a_maxBoxSize;
  }

  /*! 
    @brief Reset the coefficients for the conductivity operator (time-independent).
    @param[in] a_alpha      Alpha
    @param[in] a_beta       Beta 
    @param[in] a_acoef      a-coefficient
    @param[in] a_bcoef      b-coefficient
    @param[in] a_bcoefIrreg Irregular b-coefficient
  */
  void resetCoefficients(const Real&                                                 a_alpha,
			 const Real&                                                 a_beta,
			 const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
			 const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_bcoef,
			 const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_bcoefIrreg){
    m_alpha        = a_alpha      ;
    m_beta         = a_beta       ;
    m_aCoef = a_acoef      ;
    m_bcoef        = a_bcoef      ;
    m_bcoefIrreg   = a_bcoefIrreg ;
  }

protected:

  static int s_maxBoxSize;

  Vector< Vector<EBLevelGrid> >                                   m_eblgsMG;

  //! Time-independent A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_aCoefMG;

  //! Beginning-of-step (time-dependent) A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_aCoefMG0;

  //! End-of-step (time-dependent) A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_aCoefMG1;

  Vector< Vector< RefCountedPtr<LevelData<EBFluxFAB> > > >        m_bcoefMG;
  Vector< Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > > > m_bcoefIrregMG;

  int                                                  m_relaxType;
  int                                                  m_botdrop;
  Vector<EBLevelGrid>                                  m_eblgs;
  Vector<RefCountedPtr<EBQuadCFInterp> >               m_quadCFI;
  Vector<RefCountedPtr<EBFluxRegister> >               m_fastFR;
  Real                                                 m_alpha;
  Real                                                 m_beta;

  //! Time-independent A coefficients.
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >        m_aCoef;


  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        m_bcoef;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_bcoefIrreg;
  Real                                                 m_dxCoarse;
  Vector<int>                                          m_refRatio;
  RefCountedPtr<BaseDomainBCFactory>                   m_domainBCFactory;
  RefCountedPtr<BaseEBBCFactory>                       m_ebBCFactory;
  IntVect                                              m_ghostCellsPhi;
  IntVect                                              m_ghostCellsRhs;
  int                                                  m_numLevels;
  std::vector< bool  >                                 m_hasMGObjects;
  Vector<Real>                                         m_dx;


  void coarsenCoefficients(LevelData<EBCellFAB>               & a_acoefCoar,
			   LevelData<EBFluxFAB>               & a_bcoefCoar,
			   LevelData<BaseIVFAB<Real> >        & a_bcoefCoarIrreg,
			   const EBLevelGrid                  & a_eblgFine,
			   const EBLevelGrid                  & a_eblgCoar,
			   const LevelData<EBCellFAB>         & a_acoefFine,
			   const LevelData<EBFluxFAB>         & a_bcoefFine,
			   const LevelData<BaseIVFAB<Real> >  & a_bcoefFineIrreg,
			   const int                          & a_refToDepth);
};

#include <CD_NamespaceFooter.H>

#endif
