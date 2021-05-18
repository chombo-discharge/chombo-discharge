/* chombo-discharge
 * This is a copy of Chombos EBConductivityOpFactory and retains Chombo copyrights. 
 * Please refer to Copyright.txt in the Chombo root directory.
 */

/*!
  @file   CD_EbHelmholtzOpFactory.H
  @brief  Declares a copy of Chombos EBConductivityOpFactory
  @author Robert Marskar
*/

#ifndef CD_EbHelmholtzOpFactory_H
#define CD_EbHelmholtzOpFactory_H

// Chombo includes
#include <REAL.H>
#include <Box.H>
#include <FArrayBox.H>
#include <Vector.H>
#include <map>
#include <RefCountedPtr.H>
#include <AMRMultiGrid.H>
#include <EBIndexSpace.H>
#include <EBCellFAB.H>
#include <EBCellFactory.H>
#include <CD_EbHelmholtzOp.H>
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
#include <CD_NamespaceHeader.H>
  
//! \class EbHelmholtzOpFactory
//! Factory class to generate EbHelmholtzOps.  This follows the
//! AMRLevelOpFactory interface.
class EbHelmholtzOpFactory: public AMRLevelOpFactory<LevelData<EBCellFAB> >
{
public:

  //! Constructs a factory that builds EbHelmholtzOps with time-independent
  //! A and B coefficients.
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


  //! Destructor.
  virtual ~EbHelmholtzOpFactory();

  ///
  virtual void setData( Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_data)
  {
    m_data = a_data;
    m_dataBased = true;
  }

  ///
  virtual EbHelmholtzOp*
  MGnewOp(const ProblemDomain& a_FineindexSpace,
	  int                  a_depth,
	  bool                 a_homoOnly = true);

  EbHelmholtzOp* createOperator(const EBLevelGrid&             a_eblgMGLevel,
				   const EBLevelGrid&             a_eblgCoarMG,
				   const bool&                    a_hasMGObjects,
				   const RealVect&                a_dxMGLevel,
				   const RealVect&                a_dxCoar,
				   const int&                     a_whichLevel);
  ///
  virtual void reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim);

  ///
  virtual EbHelmholtzOp*
  AMRnewOp(const ProblemDomain& a_FineindexSpace);

  ///
  virtual void AMRreclaim(EbHelmholtzOp* a_reclaim);

  ///
  /** Refinement ratio between this level and coarser level.
      Returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToFiner(const ProblemDomain& a_domain) const;

  ///
  /**
     testRef is the size of the coarsest domain allowed in multigrid. If testRef=2,
     then the coarsest domain in multigrid will be 2x2(x2)
  **/
  static void setTestRef(int a_testRef)
  {
    s_testRef = a_testRef;
  }

  static void setMaxBoxSize(int a_maxBoxSize)
  {
    s_maxBoxSize = a_maxBoxSize;
  }

  static int s_testRef;
  static int s_maxBoxSize;

  //! Reset the coefficients for the conductivity operator (time-independent).
  void resetCoefficients(const Real&                                                 a_alpha,
			 const Real                                         &        a_beta,
			 const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
			 const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_bcoef,
			 const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_bcoefIrreg)
  {
    m_alpha       = a_alpha      ;
    m_beta        = a_beta       ;
    m_aCoefficientef       = a_acoef      ;
    m_bcoef       = a_bcoef      ;
    m_bcoefIrreg  = a_bcoefIrreg ;
  }


protected:

  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_data;
  bool m_dataBased;
  Vector< Vector<EBLevelGrid> >                                   m_eblgsMG;

  //! Time-independent A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_aCoefficientefMG;

  //! Beginning-of-step (time-dependent) A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_aCoefficientefMG0;

  //! End-of-step (time-dependent) A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_aCoefficientefMG1;

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
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >        m_aCoefficientef;


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
private:
  ///weak construction bad
  EbHelmholtzOpFactory()
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  EbHelmholtzOpFactory(const EbHelmholtzOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  void operator=(const EbHelmholtzOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

#include <CD_NamespaceFooter.H>

#endif
