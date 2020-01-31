#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LoadBalance.H"
#include "BaseIVFactory.H"
#include "CH_Timer.H"
#include "EBArith.H"
#include "MFPoissonOpFactory.H"
#include "MFPoissonOp.H"
#include "NamespaceHeader.H"

int MFPoissonOpFactory::s_testRef = 4;
int MFPoissonOpFactory::s_relaxType = 0;

////
MFPoissonOpFactory::~MFPoissonOpFactory()
{
}
///
MFPoissonOpFactory::
MFPoissonOpFactory(const RefCountedPtr<MFIndexSpace>&            a_mfis,
                   const Vector<DisjointBoxLayout>&              a_dilboVec,
                   const Vector<int>&                            a_refRatio,
                   const ProblemDomain&                          a_domainCoar,
                   const RealVect&                               a_dxCoar,
                   const RealVect&                               a_origin,
                   const Vector<RefCountedPtr<BaseDomainBC> >&   a_bc,
                   const Vector<Real>&                           a_alpha,
                   const Vector<Real>&                           a_beta,
                   const int&                                    a_ncomp,
                   const IntVect&                                a_ghostCellsPhi,
                   const IntVect&                                a_ghostCellsRHS,
                   int a_numLevels)
  : m_ghostCellsPhi( a_ghostCellsPhi ),
    m_ghostCellsRHS( a_ghostCellsRHS )
{
  m_ncomp = a_ncomp,
  m_mfis = a_mfis;
  m_bc = a_bc;
  CH_assert(a_dilboVec.size() <= a_refRatio.size());
  if (a_numLevels > 0)
    {
      m_numLevels = a_numLevels;
    }
  else
    {
      m_numLevels = a_dilboVec.size();
    }

  m_dilboVec         = a_dilboVec;
  m_refRatioVec     = a_refRatio;
  m_origin          = a_origin;
  m_dxVec.resize(m_numLevels);
  m_domainVec.resize(m_numLevels);

  m_dxVec[0] = a_dxCoar;
  m_domainVec[0] = a_domainCoar;
  for (int ilev = 1; ilev < m_numLevels; ilev++)
    {
      m_dxVec[ilev] = m_dxVec[ilev-1];
      m_dxVec[ilev] /= m_refRatioVec[ilev-1];

      m_domainVec[ilev] = m_domainVec[ilev-1];
      m_domainVec[ilev].refine(m_refRatioVec[ilev-1]);
    }

  m_alpha = a_alpha;
  m_beta = a_beta;

  // set the Neumann and Dirichlet jumps to be scalar constants
  m_scalarGD     = 0.0;
  m_scalarGN     = 0.0;
  m_isScalarJump = true;
  m_analyticJump = false;

  m_dilboVecMG.resize(m_numLevels);
  m_domainVecMG.resize(m_numLevels);
  m_hasMGObjects.resize(m_numLevels);
  m_layoutChanged.resize(m_numLevels);
  m_layoutChangedMG.resize(m_numLevels);
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      if ((ilev==0) || (m_refRatioVec[ilev] > 2))
        {
          m_hasMGObjects[ilev] = true;

          int mgRef = 2;
          m_dilboVecMG[ilev].resize(0);
          m_domainVecMG[ilev].resize(0);
          m_dilboVecMG[ilev].push_back(m_dilboVec[ilev]);
          m_domainVecMG[ilev].push_back(m_domainVec[ilev]);
          m_layoutChangedMG[ilev].resize(0);
          m_layoutChangedMG[ilev].push_back(m_layoutChanged[ilev]);

          bool hasCoarser = true;
          bool atAMRLevel = true;
          ProblemDomain curDomain = m_domainVec[ilev];
          while (hasCoarser)
            {
              int imgsize = m_dilboVecMG[ilev].size();
              const DisjointBoxLayout& dilboFine=  m_dilboVecMG[ilev][imgsize-1];

              DisjointBoxLayout dblCoarMG;
              ProblemDomain  domainCoarMG;
              int maxBoxSize = 32;
              bool layoutChanged;
              hasCoarser = EBArith::getCoarserLayouts(dblCoarMG,
                                                      domainCoarMG,
                                                      dilboFine,
                                                      curDomain,
                                                      mgRef,
                                                      maxBoxSize,
                                                      layoutChanged,
                                                      s_testRef);
              if ((atAMRLevel) && !hasCoarser)
                {
                  m_hasMGObjects[ilev] = false;
                }

              if (atAMRLevel)
                {
                  m_layoutChanged[ilev] = layoutChanged;
                  atAMRLevel= false;
                }

              if (hasCoarser)
                {
                  m_dilboVecMG[ilev].push_back(dblCoarMG);
                  m_domainVecMG[ilev].push_back(domainCoarMG);
                  m_layoutChangedMG[ilev].push_back(layoutChanged);
                  curDomain.coarsen(mgRef);
                }
            }
        }
      else
        {
          m_hasMGObjects[ilev] = false;
        }
    }
}
///
void MFPoissonOpFactory::setJump(const Real& a_gD,
                                 const Real& a_gN)
{
  CH_assert(m_ncomp == 1);
  m_scalarGD = a_gD;
  m_scalarGN = a_gN;
  m_isScalarJump = true;
  m_analyticJump = false;
}
///
void MFPoissonOpFactory::setJump(const RealVect& a_gD,
                                 const RealVect& a_gN)
{
  CH_assert(m_ncomp == SpaceDim);
  m_vectorGD = a_gD;
  m_vectorGN = a_gN;
  m_isScalarJump = false;
  m_analyticJump = false;
}
///
void MFPoissonOpFactory::setJump(const Vector< RefCountedPtr<BaseBCValue> >& a_phiValVect,
                                 const Vector< RefCountedPtr<BaseBCValue> >& a_flxValVect)
{
  m_phiValVect = a_phiValVect;
  m_flxValVect = a_flxValVect;
  m_isScalarJump = (m_ncomp == 1);
  m_analyticJump = true;
}
///
MGLevelOp<LevelData<MFCellFAB> >*
MFPoissonOpFactory::
MGnewOp(const ProblemDomain& a_domainFine,
        int                  a_depth,
        bool                 a_homoOnly)
{

  //find out if there is a real starting point here.
  int whichlev;
  bool found = false;

  for (int ilev = 0; ilev < m_numLevels && !found; ilev++)
    {
      if (a_domainFine == m_domainVec[ilev])
        {
          found = true;
          whichlev = ilev ;
        }
    }
  if (!found)
    {
      MayDay::Error("No corresponding AMRLevel to starting point of MGnewOp");
    }

  //multigrid operator.  coarn by two from depth  no
  //coarse or finer to worry about.
  DisjointBoxLayout dilboMGLevel;
  ProblemDomain domainMGLevel;
  DisjointBoxLayout dilboCoarMG;
  RealVect      dxMGLevel;
  RealVect dxCoar = RealVect::Unit;
  dxCoar *= -1.0;
  if (whichlev > 0)
    {
      dxCoar = m_dxVec[whichlev-1];
    }
  bool hasCoarMGObjects = false;
  bool coarLayoutChanged = true;
  int img = 0;
  if (a_depth == 0)
    {
      dilboMGLevel    = m_dilboVec[whichlev];
      domainMGLevel    = m_domainVec[whichlev];
      dxMGLevel      = m_dxVec[whichlev];

      hasCoarMGObjects = m_hasMGObjects[whichlev];
      if (hasCoarMGObjects)
        {
          dilboCoarMG = m_dilboVecMG[whichlev][1];
          coarLayoutChanged = m_layoutChangedMG[whichlev][1];
        }
    }
  else
    {
      int icoar = 1;
      for (int idep = 0; idep < a_depth; idep++)
        {
          icoar *= 2;
        }
      const ProblemDomain domainFine = m_domainVec[whichlev];
      ProblemDomain domainBoxMGLevel = coarsen(domainFine, icoar);
      bool foundMGLevel = false;
      int numMGLevels = m_dilboVecMG[whichlev].size();
      for (; img < numMGLevels; img++)
        {
          if (m_domainVecMG[whichlev][img] == domainBoxMGLevel)
            {
              dilboMGLevel = m_dilboVecMG[whichlev][img];
              domainMGLevel = m_domainVecMG[whichlev][img];
              foundMGLevel = true;

              hasCoarMGObjects = ((img+1) < (numMGLevels));
              if (hasCoarMGObjects)
                {
                  dilboCoarMG = m_dilboVecMG[whichlev][img+1];
                  coarLayoutChanged = m_layoutChangedMG[whichlev][img+1];
                }
              break;
            }
        }
      bool coarsenable = foundMGLevel;

      dxMGLevel = m_dxVec[whichlev];
      dxMGLevel *= Real(icoar);

      if (!coarsenable)
        {
          //not coarsenable.
          //return null
          return NULL;
        }
    }

  //creates coarse and finer info and bcs and all that
  MFPoissonOp* op = createOperator(dilboMGLevel,  dilboCoarMG, domainMGLevel, hasCoarMGObjects, coarLayoutChanged,
                                   dxMGLevel, dxCoar, whichlev, img);
  return op;
}
/////
AMRLevelOp<LevelData<MFCellFAB> >*
MFPoissonOpFactory::
AMRnewOp(const ProblemDomain& a_domainFine)
{
  //figure out which level we are at.
  int whichlev;
  bool found = false;

  for (int ilev = 0; ilev < m_numLevels && !found; ilev++)
    {
      if (a_domainFine == m_domainVec[ilev])
        {
          found = true;
          whichlev = ilev ;
        }
    }
  if (!found)
    {
      MayDay::Error("No corresponding AMRLevel to starting point of AMRnewOp");
    }

  RealVect dxCoar = RealVect::Unit;
  dxCoar *= -1.0;
  if (whichlev > 0)
    {
      dxCoar = m_dxVec[whichlev-1];
    }
  int img = 0;
  //creates coarse and finer info and bcs and all that
  DisjointBoxLayout      dilboMGLevel = m_dilboVec[whichlev];
  ProblemDomain domainMGLevel = m_domainVec[whichlev];
  RealVect           dxMGLevel =   m_dxVec[whichlev];

  bool hasCoarMGObjects = m_hasMGObjects[whichlev];
  bool coarLayoutChanged = m_layoutChanged[whichlev];
  DisjointBoxLayout dilboCoarMG;
  if (hasCoarMGObjects)
    {
      dilboCoarMG = m_dilboVecMG[whichlev][1];
    }

  MFPoissonOp* op = createOperator(dilboMGLevel, dilboCoarMG, domainMGLevel, hasCoarMGObjects, coarLayoutChanged,
                                   dxMGLevel, dxCoar, whichlev, img);
  return op;
}
//////
MFPoissonOp*
MFPoissonOpFactory::createOperator(const DisjointBoxLayout&       a_dilboMGLevel,
                                   const DisjointBoxLayout&       a_dilboCoarMG,
                                   const ProblemDomain&           a_domainMGLevel,
                                   const bool&                    a_hasMGObjects,
                                   const bool&                    a_layoutChanged,
                                   const RealVect&                a_dxMGLevel,
                                   const RealVect&                a_dxCoar,
                                   const int&                     a_whichLevel,
                                   const int&                     a_mgLevel)
{

  //fine and coarse stuff undefined because this is a MG level
  DisjointBoxLayout  dilboFine,  dilboCoar;

  int refToFine   = 1;
  int refToCoar = 1;
  bool hasFine   = ((a_whichLevel < m_numLevels - 1)  && (a_domainMGLevel == m_domainVec[a_whichLevel]));
  bool hasCoar   = ((a_whichLevel > 0) && (a_domainMGLevel == m_domainVec[a_whichLevel]));
  if (hasFine)
    {
      refToFine = m_refRatioVec[a_whichLevel];
      dilboFine  = m_dilboVec[a_whichLevel+1];
    }
  if (hasCoar)
    {
      refToCoar  =  m_refRatioVec[a_whichLevel-1];
      dilboCoar  =     m_dilboVec[a_whichLevel-1];
    }

  MFPoissonOp* op = new MFPoissonOp();
  op->define(*m_mfis, m_ncomp, a_dilboMGLevel, a_dilboCoarMG, a_hasMGObjects, a_layoutChanged,
             dilboFine, dilboCoar, a_dxMGLevel, refToCoar, refToFine, a_domainMGLevel,
             m_bc, m_ghostCellsPhi, m_ghostCellsRHS, hasCoar, hasFine, m_alpha, m_beta);

  // set the relaxation type for the operator
  op->m_relax = s_relaxType;

  // set the jump conditions for the operator
  if (m_analyticJump)
  { // jump was passed as a BaseIVFab
    op->setJump(m_phiValVect, m_flxValVect);
  }
  else
  {
    if (m_isScalarJump)
      { // jump was passed as a Real
        op->setJump(m_scalarGD, m_scalarGN);
      }
    else
      { // jump was passed as a RealVect
        op->setJump(m_vectorGN, m_vectorGN);
      }
  }

  return op;
}
///
int
MFPoissonOpFactory::
refToFiner(const ProblemDomain& a_domain) const
{
  int retval = -1;
  bool found = false;
  for (int ilev = 0; ilev < m_domainVec.size(); ilev++)
    {
      if (m_domainVec[ilev] == a_domain)
        {
          retval = m_refRatioVec[ilev];
          found = true;
        }
    }
  if (!found)
    {
      MayDay::Error("Domain not found in AMR hierarchy");
    }
  return retval;
}
/****/
void
MFPoissonOpFactory::
reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim)
{
  delete a_reclaim;
}
/****/
void
MFPoissonOpFactory::
AMRreclaim(MFPoissonOp* a_reclaim)
{
  delete a_reclaim;
}
#include "NamespaceFooter.H"
