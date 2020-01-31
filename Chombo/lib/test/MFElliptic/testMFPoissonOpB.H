#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TESTMFPOISSONOPB_H_
#define _TESTMFPOISSONOPB_H_

#include "TestUtilities.H"

#include "MFIndexSpace.H"
#include "MFCellFAB.H"
#include "MFAliasFactory.H"
#include "MFPoissonOp.H"
#include "DirichletPoissonDomainBC.H"
#include "MFGeometry.H"
#include "EBAMRIO.H"
#include "ParmParse.H"
// #include "DebugDump.H"
// #include "CH_Attach.H"
#include <cxxtest/TestSuite.h>
#include <cxxtest/GlobalFixture.h>
#include "UsingNamespace.H"

class ChomboWorld : public CxxTest::GlobalFixture
{
public:
  bool setUpWorld()
  {
#ifdef CH_MPI
    MPI_Init(_argcptr, _argvptr);
    // registerDebugger();
    pout() << "Set up for MPI run." << endl;
#endif
    return true;
  }
  bool tearDownWorld()
  {
#ifdef CH_MPI
    pout() << "Finalizing MPI." << endl;
    MPI_Finalize();
#endif
    return true;
  }
};

static ChomboWorld chomboWorld;

#define NUM_COEFFICIENTS 5
static Real s_coefficient_values[NUM_COEFFICIENTS] =
{
  1e-4, 1e-1, 1e0, 1e1, 1e4
};
#ifdef CH_USE_FLOAT
static Real s_eps = 1e-4;
#else
static Real s_eps = 1e-10;
#endif

class TestMFPoissonOpB: public CxxTest::TestSuite
{
public:

  bool setUpWorld()
  {
    return true;
  }

  bool tearDownWorld()
  {
    return true;
  }

  /*
    ----------------------------------------------------------------------
    Called before each test.
    Reads in parameters used by the various tests.
    ----------------------------------------------------------------------
  */
  void setUp()
  {
    const char* in_file = "testB.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);

    // get domain info from inputs
    m_refRatio.resize(3,2);
    readDomainInfo(m_nlevels,
                   m_domBox,
                   m_dx,
                   m_origin,
                   m_maxBoxSize);
    // set derived domain info
    m_pDomain.define(m_domBox);
    m_dxVect = RealVect::Unit;
    m_dxVect *= m_dx;

    // define data holders
    m_phi.resize(m_nlevels, NULL);
    m_rhs.resize(m_nlevels, NULL);
    m_resid.resize(m_nlevels, NULL);

    // read geometry info
    readSphereInfo(m_center, m_radius);
    // read plane geometry info. normal is not used.
    RealVect dummyVect;
    readPlaneInfo(dummyVect, m_point);
    // initialization of MFIndexSpace; will be done again later
    m_mfis = RefCountedPtr<MFIndexSpace>(new MFIndexSpace());

    m_ncomp = 1;
    m_comps.resize(2, m_ncomp);
    m_ghost = 3;
  }

  /*
    ----------------------------------------------------------------------
    Called after each test.
    ----------------------------------------------------------------------
  */
  void tearDown()
  {
  }

  /*
    ----------------------------------------------------------------------
    Test of single-level MFPoissonOp using Quadratic IBCs and sphere
    geometry.
    ----------------------------------------------------------------------
  */
  void testTruncationErrorForSphereWithQuadraticBCs()
  {
    // diagnostics
    pout() << "Sphere test..." << endl;

    // define sphere geometry and MFIndexSpace
    {
      // first get domain and dx to finest level
      ProblemDomain finestDomain(m_domBox);
      RealVect finestDx = m_dxVect;
      for (int ilev=0; ilev<(m_nlevels-1); ++ilev)
        {
          finestDomain.refine(m_refRatio[ilev]);
          finestDx /= m_refRatio[ilev];
        }
      m_ifunc = MFGeometry::defineSphereGeometry(m_mfis,
                                                 finestDomain,
                                                 finestDx,
                                                 m_origin,
                                                 m_center,
                                                 m_radius);
    }

    //make grids
    IntVectSet tags(m_domBox);
    tags.grow(1);
    makeHierarchy(m_dblv, m_pDomain, tags, m_refRatio, m_nlevels, m_maxBoxSize);

    // define data holders (2)
    {
      ProblemDomain domLev(m_domBox);
      for (int ilev=0; ilev<m_nlevels; ++ilev)
        {
          MFCellFactory factory(*m_mfis, m_dblv[ilev], domLev.domainBox(), m_comps, m_ghost);
          m_phi[ilev]   = new LevelData<MFCellFAB>(m_dblv[ilev], m_ncomp, m_ghost*IntVect::Unit, factory);
          m_rhs[ilev]   = new LevelData<MFCellFAB>(m_dblv[ilev], m_ncomp,         IntVect::Zero, factory);
          m_resid[ilev] = new LevelData<MFCellFAB>(m_dblv[ilev], m_ncomp,         IntVect::Zero, factory);
          domLev.refine(m_refRatio[ilev]);
        }
    }

    // set values for all coefficients
    Real gD = 0.;
    Real gN = 0.;
    Vector<Real> alpha(2), beta(2);

    // For when I implement a range of alpha coefficients
    /*
    for (int ia0=0; ia0<NUM_COEFFICIENTS; ia0++)
      {
        alpha[0] = s_coefficient_values[ia0];
        for (int ia1=0; ia1<NUM_COEFFICIENTS; ia1++)
          {
            alpha[1] = s_coefficient_values[ia1];
    */
    alpha[0] = 0.;
    alpha[1] = 0.;

    for (int ib0=0; ib0<NUM_COEFFICIENTS; ib0++)
      {
        beta[0] = s_coefficient_values[ib0];
        for (int ib1=0; ib1<NUM_COEFFICIENTS; ib1++)
          {
            beta[1] = s_coefficient_values[ib1];

            // BC setting and data holder initialization
            Vector<RefCountedPtr<BaseBCValue> > bcv;
            // create BCValue for inside (phase=0)
            Real iConstant = 1.0;
            RealVect iLinearCoefficients = RealVect::Zero;
            RealVect iXXCoefficients = RealVect::Unit/beta[0];
            RealVect iXYCoefficients = RealVect::Zero;
            RefCountedPtr<BaseBCValue> ival( new TriquadraticDirichletBC(iConstant,
                                                                         iLinearCoefficients,
                                                                         iXXCoefficients,
                                                                         iXYCoefficients) );
            // add this BCValue to bcv for passing to initialization routine
            bcv.push_back(ival);

            // create BCValue for outside (phase=1), using inside values
            RealVect oLinearCoefficients = iLinearCoefficients*(beta[0]/beta[1]);
            RealVect oXXCoefficients = iXXCoefficients*(beta[0]/beta[1]);
            RealVect oXYCoefficients = iXYCoefficients*(beta[0]/beta[1]);
            RealVect boundaryPoint = m_center + m_radius*BASISREALV(0);
            Real iBoundaryValue = ival->value(boundaryPoint, RealVect::Zero, 0.0, 0);
            TriquadraticDirichletBC oHomogeneousBCValue(0.,
                                                        oLinearCoefficients,
                                                        oXXCoefficients,
                                                        oXYCoefficients);
            Real oHomogeneousBoundaryValue = oHomogeneousBCValue.value(boundaryPoint,
                                                                       RealVect::Zero,
                                                                       0.0,
                                                                       0);
            Real oConstant = gD + iBoundaryValue - oHomogeneousBoundaryValue;
            RefCountedPtr<BaseBCValue> oval(new TriquadraticDirichletBC(oConstant,
                                                                        oLinearCoefficients,
                                                                        oXXCoefficients,
                                                                        oXYCoefficients) );
            // add this BCValue to bcv for passing to initialization routine
            bcv.push_back(oval);

            // create vector containing the (constant) value of the RHS in each phase
            Vector<Real> rhsVals(2);
            rhsVals[0] = 2.*(beta[0])*iXXCoefficients.sum();
            rhsVals[1] = 2.*(beta[1])*oXXCoefficients.sum();

            // initialize ICs and BCs
            initializeQuadraticDirichletIBC(m_bc, bcv, rhsVals);

            // eps needs to be scaled for these coefficients
            Real scaledEps = Max(beta[0], beta[1])*s_eps;

            // diagnostics
            pout() << "  Test results for alpha=[" << alpha[0] << ", " << alpha[1] << "] "
                   << "beta=[" << beta[0] << ", " << beta[1] << "]" << endl;
            pout() << "           Scaled EPS = " << scaledEps << endl;

            // if we pass all level tests, then passTest is true
            bool passTest = true;
            // loop over levels, running the test
            RealVect dxLev = m_dxVect;
            ProblemDomain domLev(m_domBox);
            for (int ilev=0; ilev<m_nlevels; ++ilev)
              {
                Real r0norm, r1norm, r2norm;
                runLevelTest(r0norm, r1norm, r2norm,
                             domLev, dxLev, ilev,
                             alpha, beta, gD, gN);

                // cxxTest assertion for testing
                // TS_ASSERT( (r0norm < scaledEps) && (r1norm < scaledEps) && (r2norm < scaledEps) );

                // did we pass the test on this level? If not, write out diagnostics.
                bool passTestLevel = (r0norm < scaledEps) &&
                  (r1norm < scaledEps) &&
                  (r2norm < scaledEps);
                // update passTest with result on this level
                passTest = passTest && passTestLevel;

                // diagnostics
                if (!passTestLevel)
                  {
                    pout() << "    level " << ilev
                           << "; L0norm=" << r0norm
                           << "  L1norm=" << r1norm
                           << "  L2norm=" << r2norm << endl;
                  }

                // update domain and dx for next level
                domLev.refine(m_refRatio[ilev]);
                dxLev /= m_refRatio[ilev];
              }

            // Temporary cxxTest assertion, for all levels
            TS_ASSERT(passTest);

            // write out HDF5 file if we failed any of the level tests
            if (!passTest)
              {
                int id = 0;
                writeDataOnTestFailure("sphere", id, alpha, beta);
              }

          } // end loop over beta1
      } // end loop over beta0

    // clean up
    for (int ilev=0; ilev<m_nlevels; ilev++)
      {
        delete m_phi[ilev];
        delete m_rhs[ilev];
        delete m_resid[ilev];
      }
  }

  /*
    ----------------------------------------------------------------------
    Test of single-level MFPoissonOp using Quadratic IBCs and plane
    geometry.
    ----------------------------------------------------------------------
  */
  void testTruncationErrorForPlaneWithQuadraticBCs()
  {
    // diagnostics
    pout() << "Plane test..." << endl;

    //make grids
    IntVectSet tags(m_domBox);
    tags.grow(1);
    makeHierarchy(m_dblv, m_pDomain, tags, m_refRatio, m_nlevels, m_maxBoxSize);

    // loop over directions; normals will be along axes
    for (int idir=0; idir<SpaceDim; ++idir)
      {
        RealVect normal = BASISREALV(idir);

        // define plane geometry and MFIndexSpace
        {
          // first get domain and dx to finest level
          ProblemDomain finestDomain(m_domBox);
          RealVect finestDx = m_dxVect;
          for (int ilev=0; ilev<(m_nlevels-1); ++ilev)
            {
              finestDomain.refine(m_refRatio[ilev]);
              finestDx /= m_refRatio[ilev];
            }
          m_ifunc = MFGeometry::definePlaneGeometry(m_mfis,
                                                    finestDomain,
                                                    finestDx,
                                                    m_origin,
                                                    normal,
                                                    m_point);
        }

        // define data holders (2)
        {
          ProblemDomain domLev(m_domBox);
          for (int ilev=0; ilev<m_nlevels; ++ilev)
            {
              MFCellFactory factory(*m_mfis, m_dblv[ilev], domLev.domainBox(), m_comps, m_ghost);
              m_phi[ilev]   = new LevelData<MFCellFAB>(m_dblv[ilev], m_ncomp, m_ghost*IntVect::Unit, factory);
              m_rhs[ilev]   = new LevelData<MFCellFAB>(m_dblv[ilev], m_ncomp,         IntVect::Zero, factory);
              m_resid[ilev] = new LevelData<MFCellFAB>(m_dblv[ilev], m_ncomp,         IntVect::Zero, factory);
              domLev.refine(m_refRatio[ilev]);
            }
        }

        // set values for all coefficients
        Real gD = 0.;
        Real gN = 0.;
        Vector<Real> alpha(2), beta(2);

        // For when I implement a range of alpha coefficients
        /*
        for (int ia0=0; ia0<NUM_COEFFICIENTS; ia0++)
          {
            alpha[0] = s_coefficient_values[ia0];
            for (int ia1=0; ia1<NUM_COEFFICIENTS; ia1++)
              {
                alpha[1] = s_coefficient_values[ia1];
        */
        alpha[0] = 0.;
        alpha[1] = 0.;

        for (int ib0=0; ib0<NUM_COEFFICIENTS; ib0++)
          {
            beta[0] = s_coefficient_values[ib0];
            for (int ib1=0; ib1<NUM_COEFFICIENTS; ib1++)
              {
                beta[1] = s_coefficient_values[ib1];

                // BC setting and data holder initialization
                Vector< RefCountedPtr<BaseBCValue> > bcv;
                // create BCValue for inside (phase=0)
                Real iConstant = 1.0;
                RealVect iLinearCoefficients = RealVect::Zero;
                RealVect iXXCoefficients = normal/beta[0];
                RealVect iXYCoefficients = RealVect::Zero;
                RefCountedPtr<BaseBCValue> ival( new TriquadraticDirichletBC(iConstant,
                                                                             iLinearCoefficients,
                                                                             iXXCoefficients,
                                                                             iXYCoefficients) );
                // add this BCValue to bcv for passing to initialization routine
                bcv.push_back(ival);

                // create BCValue for outside (phase=1), using inside values
                RealVect oLinearCoefficients = iLinearCoefficients*(beta[0]/beta[1]);
                RealVect oXXCoefficients = iXXCoefficients*(beta[0]/beta[1]);
                RealVect oXYCoefficients = iXYCoefficients*(beta[0]/beta[1]);
                RealVect boundaryPoint = m_point;
                Real iBoundaryValue = ival->value(boundaryPoint, RealVect::Zero, 0.0, 0);
                TriquadraticDirichletBC oHomogeneousBCValue(0.,
                                                            oLinearCoefficients,
                                                            oXXCoefficients,
                                                            oXYCoefficients);
                Real oHomogeneousBoundaryValue = oHomogeneousBCValue.value(boundaryPoint,
                                                                           RealVect::Zero,
                                                                           0.0,
                                                                           0);
                Real oConstant = gD + iBoundaryValue - oHomogeneousBoundaryValue;
                RefCountedPtr<BaseBCValue> oval( new TriquadraticDirichletBC(oConstant,
                                                                             oLinearCoefficients,
                                                                             oXXCoefficients,
                                                                             oXYCoefficients) );
                // add this BCValue to bcv for passing to initialization routine
                bcv.push_back(oval);

                // create vector containing the (constant) value of the RHS in each phase
                Vector<Real> rhsVals(2);
                rhsVals[0] = 2.*(beta[0])*iXXCoefficients.sum();
                rhsVals[1] = 2.*(beta[1])*oXXCoefficients.sum();

                // initialize ICs and BCs
                initializeQuadraticDirichletIBC(m_bc, bcv, rhsVals);

                // eps needs to be scaled for these coefficients
                Real scaledEps = Max(beta[0], beta[1])*s_eps;

                // diagnostics
                pout() << "  Test results for alpha=[" << alpha[0] << ", " << alpha[1] << "] "
                       << "beta=[" << beta[0] << ", " << beta[1] << "]" << endl;
                pout() << "           Scaled EPS = " << scaledEps << endl;

                // if we pass all level tests, then passTest is true
                bool passTest = true;
                // loop over levels, running the test
                RealVect dxLev = m_dxVect;
                ProblemDomain domLev(m_domBox);
                for (int ilev=0; ilev<m_nlevels; ++ilev)
                  {
                    Real r0norm, r1norm, r2norm;
                    runLevelTest(r0norm, r1norm, r2norm,
                                 domLev, dxLev, ilev,
                                 alpha, beta, gD, gN);

                    // cxxTest assertion for testing
                    // TS_ASSERT( (r0norm < scaledEps) && (r1norm < scaledEps) && (r2norm < scaledEps) );

                    // did we pass the test on this level? If not, write out diagnostics.
                    bool passTestLevel = (r0norm < scaledEps) &&
                      (r1norm < scaledEps) &&
                      (r2norm < scaledEps);
                    // update passTest with result on this level
                    passTest = passTest && passTestLevel;

                    // diagnostics
                    if (!passTestLevel)
                      {
                        pout() << "    level " << ilev
                               << "; L0norm=" << r0norm
                               << "  L1norm=" << r1norm
                               << "  L2norm=" << r2norm << endl;
                      }

                    // update domain and dx for next level
                    domLev.refine(m_refRatio[ilev]);
                    dxLev /= m_refRatio[ilev];
                  }

                // Temporary cxxTest assertion, for all levels
                TS_ASSERT(passTest);

                // write out HDF5 file if we failed any of the level tests
                if (!passTest)
                  {
                    writeDataOnTestFailure("plane", idir, alpha, beta);
                  }

              } // end loop over beta1
          } // end loop over beta0

        // clean up
        for (int ilev=0; ilev<m_nlevels; ilev++)
          {
            delete m_phi[ilev];
            delete m_rhs[ilev];
            delete m_resid[ilev];
          }
      } // end loop over normal directions
  }

private:

  /*
    ----------------------------------------------------------------------
    Initializes the domain BCs, initial phi, and rhs=L(phi).
    ----------------------------------------------------------------------
  */
  void initializeQuadraticDirichletIBC(      Vector< RefCountedPtr<BaseDomainBC> >&  a_bcs,
                                       const Vector< RefCountedPtr<BaseBCValue> >&   a_bcVals,
                                       const Vector<Real>&                           a_rhsVals)
  {
    a_bcs.clear();
    RefCountedPtr<DirichletPoissonDomainBC> insideBC( new DirichletPoissonDomainBC() );
    RefCountedPtr<DirichletPoissonDomainBC> outsideBC( new DirichletPoissonDomainBC() );
    a_bcs.push_back(insideBC);
    a_bcs.push_back(outsideBC);
    insideBC->setFunction(a_bcVals[0]);
    outsideBC->setFunction(a_bcVals[1]);

    //initialize data holders
    ProblemDomain domLev(m_domBox);
    RealVect dxLev = m_dxVect;
    for (int ilev=0; ilev<m_nlevels; ++ilev)
      {
        LevelData<EBCellFAB> phase0, phase1;
        aliasMF(phase0, 0, *m_phi[ilev]);
        aliasMF(phase1, 1, *m_phi[ilev]);
        setValue(phase0, *(a_bcVals[0]), domLev.domainBox(), dxLev, m_origin);
        setValue(phase1, *(a_bcVals[1]), domLev.domainBox(), dxLev, m_origin);
        aliasMF(phase0, 0, *m_rhs[ilev]);
        aliasMF(phase1, 1, *m_rhs[ilev]);
        setValue(phase0, ConstantDirichletBC(a_rhsVals[0]),
                 domLev.domainBox(), dxLev, m_origin, true);
        setValue(phase1, ConstantDirichletBC(a_rhsVals[1]),
                 domLev.domainBox(), dxLev, m_origin, true);
        // update domain and dx for next level
        domLev.refine(m_refRatio[ilev]);
        dxLev /= m_refRatio[ilev];
      }
  }

  /*
    ----------------------------------------------------------------------
    Runs the test at a given level, returning the norms in the
    corresponding inputs.
    ----------------------------------------------------------------------
  */
  void runLevelTest(      Real&               a_norm0,
                          Real&               a_norm1,
                          Real&               a_norm2,
                    const ProblemDomain&      a_domLev,
                    const RealVect&           a_dx,
                    const int                 a_ilev,
                    const Vector<Real>&       a_alpha,
                    const Vector<Real>&       a_beta,
                    const Real                a_gD,
                    const Real                a_gN)
  {
    MFPoissonOp op;
    DisjointBoxLayout dumdbl;
    op.define(*m_mfis, m_ncomp, m_dblv[a_ilev], dumdbl, false, true, dumdbl, dumdbl,
              a_dx, 2, 2, a_domLev, m_bc, m_ghost*IntVect::Unit, IntVect::Zero,
              false, false, a_alpha, a_beta);

    op.setJump(a_gD, a_gN);
    op.setToZero(*m_resid[a_ilev]);

    op.residual(*m_resid[a_ilev], *m_phi[a_ilev], *m_rhs[a_ilev]);
    zeroAdjacentToBoundary(*m_resid[a_ilev], a_domLev);
    a_norm0 = op.norm(*m_resid[a_ilev], 0);
    a_norm1 = op.norm(*m_resid[a_ilev], 1);
    a_norm2 = op.norm(*m_resid[a_ilev], 2);
  }

  /*
    ----------------------------------------------------------------------
    Writes out phi and residual to HDF5 file.
    ----------------------------------------------------------------------
  */
  void writeDataOnTestFailure(const char *   a_prefix,
                              const int&     a_id,
                              Vector<Real>&  a_alpha,
                              Vector<Real>&  a_beta)
  {
#ifdef CH_USE_HDF5
    ProblemDomain domLev(m_domBox);
    Vector< LevelData<FArrayBox>* > levelSet(m_nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > phaseA(m_nlevels, NULL), phaseB(m_nlevels,NULL);
    for (int ilev=0; ilev<m_nlevels; ++ilev)
      {
        MFCellFactory factory(*m_mfis, m_dblv[ilev], domLev.domainBox(), m_comps, m_ghost);
        levelSet[ilev] = new LevelData<FArrayBox>(m_dblv[ilev], m_ncomp, IntVect::Zero);
        phaseA[ilev]   = new LevelData<EBCellFAB>();
        phaseB[ilev]   = new LevelData<EBCellFAB>();
        domLev.refine(m_refRatio[ilev]);
      }

    // initialize level set (used for writing to HDF)
    MFGeometry::makeLevelSet(levelSet,
                             *m_ifunc,
                             m_origin,
                             m_dxVect,
                             m_refRatio);

    // write residual and phi out to HDF
    aliasMF(phaseA, 0, m_resid);
    aliasMF(phaseB, 1, m_resid);
    char iter_str[80];
    // For when I implement a range of alpha coefficients
    /*
    int a0 = (int)log10(a_alpha[0]);
    int a1 = (int)log10(a_alpha[1]);
    */
    int a0=0;
    int a1=0;
    int b0 = (int)log10( a_beta[0]);
    int b1 = (int)log10( a_beta[1]);
    sprintf(iter_str, "%s_residual_alpha%+1d%+1d_beta%+1d%+1d.%03d.%dd.hdf5",
            a_prefix, a0, a1, b0, b1, a_id, SpaceDim);
    Vector<std::string> names(2); names[0]="residual0"; names[1]="residual1";
    int step = 0;
    Vector<Real> coveredVal(m_ncomp, 1.2345678e90);
    writeEBHDF5(iter_str, m_dblv, &phaseA, &phaseB, &levelSet, names, m_pDomain,
                m_dx, 1, step, m_refRatio, m_nlevels, true, coveredVal);

    aliasMF(phaseA, 0, m_phi);
    aliasMF(phaseB, 1, m_phi);
    sprintf(iter_str, "%s_phi_alpha%+1d%+1d_beta%+1d%+1d.%03d.%dd.hdf5",
            a_prefix, a0, a1, b0, b1, a_id, SpaceDim);
    names[0]="phi0"; names[1]="phi1";
    writeEBHDF5(iter_str, m_dblv, &phaseA, &phaseB, &levelSet, names, m_pDomain,
                m_dx, 1, step, m_refRatio, m_nlevels, true, coveredVal);

    for (int ilev=0; ilev<m_nlevels; ilev++)
      {
        delete phaseA[ilev];
        delete phaseB[ilev];
        delete levelSet[ilev];
      }
#endif
  }

private:
  Vector<int> m_refRatio;
  int m_nlevels;
  int m_maxBoxSize;
  Box m_domBox;
  RealVect m_origin;
  Real m_dx;
  RealVect m_dxVect;
  ProblemDomain m_pDomain;
  Vector< LevelData<MFCellFAB>* > m_phi, m_rhs, m_resid;
  RefCountedPtr<MFIndexSpace> m_mfis;
  int m_ncomp;
  Vector<int> m_comps;
  int m_ghost;
  Vector<DisjointBoxLayout> m_dblv;
  Real m_radius;
  RealVect m_center;
  RealVect m_point;
  RefCountedPtr<BaseIF> m_ifunc;
  Vector<RefCountedPtr<BaseDomainBC> > m_bc;
};

#endif
