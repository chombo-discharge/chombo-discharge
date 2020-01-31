#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Sat, July 18, 1999

#include "LevelMG.H"
#include "LevelMGF_F.H"
#include "NamespaceHeader.H"

void LevelMG::setDefaultValues()
{
  // set defaults
  m_numBottomGSRB = 4;
  m_numSmoothUp = 4;
  m_numSmoothDown = 4;

  m_lCoarsePtr = NULL;
  m_levelopPtr = NULL;

  m_dxLevel = -1;

  m_isDefined = false;
}

bool LevelMG::isDefined() const
{
  return m_isDefined;
}

void LevelMG::clearMemory()
{
  if (m_lCoarsePtr != NULL)
    {
      delete m_lCoarsePtr;
      m_lCoarsePtr = NULL;
    }

  if (m_levelopPtr != NULL)
    {
      delete m_levelopPtr;
      m_levelopPtr = NULL;
    }

  m_isDefined = false;
}

// default constructor
LevelMG::LevelMG()
{
  setDefaultValues();
}

// Constructor with DisjointBoxLayout, number of coarser levels
// this calls the levelmg big define function
LevelMG::LevelMG(const DisjointBoxLayout& a_ba,
                 const DisjointBoxLayout* a_baseBaPtr,
                 Real                     a_dxLevel,
                 int                      a_refRatio,
                 const Box&               a_domain,
                 int                      a_nCoarserLevels,
                 const LevelOp* const     a_opin,
                 int                      a_ncomp)
{
  setDefaultValues();

  ProblemDomain physdomain(a_domain);

  define(a_ba,a_baseBaPtr, a_dxLevel, a_refRatio, physdomain,
         a_nCoarserLevels, a_opin, a_ncomp);
}

// Constructor with DisjointBoxLayout, number of coarser levels
// this calls the levelmg big define function
LevelMG::LevelMG(const DisjointBoxLayout& a_ba,
                 const DisjointBoxLayout* a_baseBaPtr,
                 Real                     a_dxLevel,
                 int                      a_refRatio,
                 const ProblemDomain&     a_domain,
                 int                      a_nCoarserLevels,
                 const LevelOp* const     a_opin,
                 int                      a_ncomp)
{
  setDefaultValues();

  define(a_ba,a_baseBaPtr, a_dxLevel, a_refRatio, a_domain,
         a_nCoarserLevels, a_opin, a_ncomp);
}

// corresponding define function
void LevelMG::define(const DisjointBoxLayout& a_ba,
                     const DisjointBoxLayout* a_baseBaPtr,
                     Real                     a_dxLevel,
                     int                      a_refRatio,
                     const  Box&              a_domain,
                     int                      a_nCoarserLevels,
                     const LevelOp* const     a_opin,
                     int                      a_ncomp)
{
  ProblemDomain physdomain(a_domain);

  define(a_ba, a_baseBaPtr, a_dxLevel, a_refRatio, a_domain,
         a_nCoarserLevels, a_opin, a_ncomp);
}

// corresponding define function
void LevelMG::define(const DisjointBoxLayout& a_ba,
                     const DisjointBoxLayout* a_baseBaPtr,
                     Real                     a_dxLevel,
                     int                      a_refRatio,
                     const  ProblemDomain&    a_domain,
                     int                      a_nCoarserLevels,
                     const LevelOp* const     a_opin,
                     int                      a_ncomp)
{
  CH_assert(a_opin != NULL);
  CH_assert(a_nCoarserLevels >= 0);
  CH_assert(!a_domain.isEmpty());
  CH_assert(a_ba.checkPeriodic(a_domain));

  //either there are no coarser levels
  //or we have a legitimate ref ratio
  CH_assert((a_refRatio > 0) || (a_baseBaPtr == NULL));
  CH_assert(a_dxLevel > 0);

  clearMemory();

  m_ba = a_ba;
  m_baseBaPtr = a_baseBaPtr;

  m_dxLevel = a_dxLevel;
  m_domain  = a_domain;

  m_nCoarserLevels = a_nCoarserLevels;

  m_isDefined = true;
  m_resid.define(m_ba, a_ncomp, IntVect::Zero);
  m_levelopPtr = a_opin->new_levelop();

  // no need to do inhomgeneous CFInterp here
  // all levelmg operations (smooth, etc) use
  // only homogeneous cf interpolation
  bool homogeneousOnly = true;

  m_levelopPtr->define(m_ba,  m_baseBaPtr,
                       m_dxLevel,  a_refRatio,
                       m_domain, homogeneousOnly, a_ncomp);
  // recursively define MG hierarchy
  if (m_nCoarserLevels != 0)
    {

      //this constructor is called at the finest mg level.
      //this means that it might need an averaging operator
      //but it does not need an interpolation operator
      m_refToCoar = 2;

      if (m_baCoarsened.isClosed())
      {
        m_baCoarsened = DisjointBoxLayout();
      }

      coarsen(m_baCoarsened, m_ba, m_refToCoar);

      m_crseCorr .define(m_baCoarsened, a_ncomp, IntVect::Unit);
      m_crseResid.define(m_baCoarsened, a_ncomp, IntVect::Zero);

      m_averageOp.define(m_ba, m_baCoarsened, a_ncomp, m_refToCoar);

      m_lCoarsePtr = new LevelMG(*this, m_refToCoar, a_opin);

      if (m_lCoarsePtr == NULL)
        {
          MayDay::Error("Out of Memory in LevelMG::define");
        }
    }
}

// --------------------------------------------------------------
//  constructor for coarsened version of levelMG
// refcoarse is the ratio from coarsened thing to this thing
// --------------------------------------------------------------
LevelMG::LevelMG(const LevelMG& a_level,
                 int            a_refCoarse,
                 const LevelOp* a_opin)
{
  // set pointers to NULL so define does not try to delete them
  setDefaultValues();

  define(a_level, a_refCoarse, a_opin);
}

// --------------------------------------------------------------
// constructor for coarsened version of levelMG
// refcoarse is the ratio from coarsened thing to this thing
// --------------------------------------------------------------
void LevelMG::define(const LevelMG&       a_level,
                     int                  a_refCoarse,
                     const LevelOp* const a_opin)
{
  CH_assert(a_refCoarse == 2);
  CH_assert(a_opin != NULL);

  clearMemory();

  m_isDefined = true;
  m_ba = a_level.m_baCoarsened;

  Real dxcoarse = a_refCoarse*a_level.m_dxLevel;
  int ncomp = a_level.m_resid.nComp();

  ProblemDomain dprobc = coarsen(a_level.m_domain, a_refCoarse);

  m_nCoarserLevels = a_level.m_nCoarserLevels - 1;

  CH_assert(m_nCoarserLevels >= 0);

  m_baseBaPtr = a_level.m_baseBaPtr;
  m_dxLevel = dxcoarse;
  m_domain = dprobc;

  m_resid.define(m_ba, ncomp, IntVect::Zero);

  m_levelopPtr = a_level.m_levelopPtr->new_levelop();

  m_levelopPtr->define(a_level.m_levelopPtr, a_refCoarse);

  m_numSmoothDown = a_level.m_numSmoothDown;
  m_numSmoothUp = a_level.m_numSmoothUp;

  // recursively define MG hierarchy
  if (m_nCoarserLevels != 0)
    {
      m_refToCoar = 2;

      coarsen(m_baCoarsened, m_ba, m_refToCoar);

      m_crseCorr .define(m_baCoarsened, ncomp, IntVect::Unit);
      m_crseResid.define(m_baCoarsened, ncomp, IntVect::Zero);

      m_averageOp.define(m_ba, m_baCoarsened, ncomp, m_refToCoar);

      m_lCoarsePtr = new LevelMG(*this, m_refToCoar, a_opin);
      if (m_lCoarsePtr == NULL)
        {
          MayDay::Error("Out of Memory in LevelMG::define");
        }
    }
  else
    {
      m_lCoarsePtr = NULL;
    }
}

// ------------------------------------------------------------
// destructor
// ------------------------------------------------------------
LevelMG::~LevelMG()
{
  clearMemory();
}

// -------------------------------------------------------------
void
LevelMG::mgRelax(LevelData<FArrayBox>&       a_soln,
                 const LevelData<FArrayBox>& a_rhs,
                 bool                        a_bottomSolveFlag)
{
  CH_assert(isDefined());

  int ncomp = a_rhs.nComp();

  CH_assert(a_soln.nComp() == ncomp);
  CH_assert (ncomp == m_resid.nComp());
  CH_assert( a_soln.ghostVect() >= IntVect::Unit);

  // now, either coarsen or call bottom solver
  if (m_lCoarsePtr != NULL)
    {
      // first, smooth on this level
      for (int iter=0; iter< m_numSmoothDown; iter++)
        {
          m_levelopPtr->smooth(a_soln,a_rhs);
        }

      DataIterator ditCoar = m_crseCorr.dataIterator();
      for (ditCoar.reset(); ditCoar.ok(); ++ditCoar)
        {
          m_crseCorr[ditCoar()].setVal(0.0);
        }

      // define residual on this grid
      // using homogeneous cf bcs
      m_levelopPtr->applyOpH(a_soln,m_resid);

      DataIterator dit = m_resid.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          m_resid[dit()] -= a_rhs[dit()];
          m_resid[dit()].negate();
        }

      // coarsen residual and physical boundary condition
      CH_assert(m_averageOp.isDefined());
      m_averageOp.averageToCoarse(m_crseResid, m_resid);

      // relax recursively on coarsened level
      m_lCoarsePtr->mgRelax(m_crseCorr, m_crseResid, a_bottomSolveFlag);

      // interpolate correction back up and modify solution
      crseCorrect(a_soln, m_crseCorr, m_refToCoar);

      // smooth again on this level
      for (int iter=0; iter < m_numSmoothUp; iter++)
        {
          m_levelopPtr->smooth(a_soln,a_rhs);
        }
    }
  else
  if (m_domain.domainBox().numPts() == 1)
    {
      // grid is 1x1, just need to smooth
      for (int iter = 0; iter < m_numBottomGSRB; iter++)
        {
          m_levelopPtr->smooth(a_soln,a_rhs);
        }
    }
  else
  if (a_bottomSolveFlag)
    {
      // first, smooth on this level
      for (int iter = 0; iter < m_numBottomGSRB; iter++)
        {
          m_levelopPtr->smooth(a_soln,a_rhs);
        }

      // then call bottom smoother
      m_levelopPtr->bottomSmoother(a_soln,a_rhs);
    }
  else
    {
      // if we get here, then there is no coarser level, but
      // we're not at the bottom.  This means that we're at the bottom
      // of a mini-V-cycle between AMR levels.  In this case, we just
      // call smooth on this level.  Do this max(numSmoothUp,numSmoothDown)
      // in case one or the other is zero
      int numSmooth = std::max(m_numSmoothUp, m_numSmoothDown);

      // do numSmoothDown iterations
      for (int iter = 0; iter < numSmooth; iter++)
        {
          m_levelopPtr->smooth(a_soln,a_rhs);
        }
    }
}

void LevelMG::crseCorrect(LevelData<FArrayBox>&       a_fine,
                          const LevelData<FArrayBox>& a_crse,
                          int                         a_refRat)
{
  DataIterator crseIt = a_crse.dataIterator();

  const DisjointBoxLayout& crseBoxes = a_crse.getBoxes();
  const DisjointBoxLayout& fineBoxes = a_fine.getBoxes();

  for (crseIt.reset(); crseIt.ok(); ++crseIt)
    {
      const Box& crseBox = crseBoxes[crseIt()];
      const FArrayBox& crseFab = a_crse[crseIt()];

      FArrayBox& fineFab = a_fine[crseIt()];
      Box fineBox = fineBoxes[crseIt()];

      fineBox.coarsen(a_refRat);
      fineBox &= crseBox;

      // c set refinement box
      Box nrefbox(IntVect::Zero,
                  (a_refRat-1)*IntVect::Unit);

      if (!fineBox.isEmpty())
        {
          FORT_INTERPMG(CHF_FRA(fineFab),
                        CHF_FRA(crseFab),
                        CHF_BOX(fineBox),
                        CHF_CONST_INT(a_refRat),
                        CHF_BOX(nrefbox));
         }
    }
}

LevelOp* LevelMG::levelOpPtr()
{
  return m_levelopPtr;
}

LevelMG* LevelMG::lCoarsePtr()
{
  return m_lCoarsePtr;
}

/// Set number of smoothing steps on the way up V-cycle
void
LevelMG::setNumSmoothUp(int a_numSmoothUp)
{
  m_numSmoothUp = a_numSmoothUp;
  // may need to pass this down the chain
  if (lCoarsePtr() != NULL)
    {
      lCoarsePtr()->setNumSmoothUp(a_numSmoothUp);
    }
}

/// Set number of smoothing steps on the way down V-cycle
void
LevelMG::setNumSmoothDown(int a_numSmoothDown)
{
  m_numSmoothDown = a_numSmoothDown;
  // may need to pass this down the chain
  if (lCoarsePtr() != NULL)
    {
      lCoarsePtr()->setNumSmoothDown(a_numSmoothDown);
    }
}

void
LevelMG::setConvergenceMetric(Real a_metric, int a_comp)
{
  // first set on this level, then pass down to coarser level
  CH_assert (m_levelopPtr != NULL);
  m_levelopPtr->setConvergenceMetric(a_metric, a_comp);

  if (m_lCoarsePtr != NULL)
    {
      m_lCoarsePtr->setConvergenceMetric(a_metric, a_comp);
    }
}
#include "NamespaceFooter.H"
