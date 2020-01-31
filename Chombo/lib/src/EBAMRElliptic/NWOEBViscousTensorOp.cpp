#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NWOEBViscousTensorOp.H"
#include "EBAMRPoissonOpF_F.H"
#include "EBAMRPoissonOp.H"
#include "InterpF_F.H"
#include "EBArith.H"
#include "EBViscousTensorOpF_F.H"
#include "ViscousTensorOpF_F.H"
#include "NWOViscousTensorOp.H"
#include "AMRPoissonOpF_F.H"
#include "EBAMRPoissonOpF_F.H"
#include "EBCoarseAverage.H"
#include "ParmParse.H"
#include "EBLevelDataOpsF_F.H"
#include "EBCellFactory.H"
#include "EBFluxFactory.H"
#include "BaseIVFactory.H"
#include "EBAlias.H"
#include "BCFunc.H"
#include "DebugOut.H"
#include "NWOEBQuadCFInterp.H"
#include "NamespaceHeader.H"
bool NWOEBViscousTensorOp::s_doLazyRelax = true;
                                         
int NWOEBViscousTensorOp::s_whichLev = -1;
int NWOEBViscousTensorOp::s_step = -1;
bool NWOEBViscousTensorOp::s_turnOffBCs = false; //really needs to default to false


bool NWOEBViscousTensorOp::s_forceNoEBCF = true;

//-----------------------------------------------------------------------
NWOEBViscousTensorOp::
NWOEBViscousTensorOp(const EBLevelGrid &                                a_eblgFine,
                     const EBLevelGrid &                                a_eblg,
                     const EBLevelGrid &                                a_eblgCoar,
                     const EBLevelGrid &                                a_eblgCoarMG,
                     const Real&                                        a_alpha,
                     const Real&                                        a_beta,
                     const RefCountedPtr<LevelData<EBCellFAB> >&        a_acoef,
                     const RefCountedPtr<LevelData<EBFluxFAB> >&        a_eta,
                     const RefCountedPtr<LevelData<EBFluxFAB> >&        a_lambda,
                     const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_etaIrreg,
                     const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_lambdaIrreg,
                     const Real&                                        a_dx,
                     const Real&                                        a_dxCoar,
                     const int&                                         a_refToFine,
                     const int&                                         a_refToCoar,
                     const RefCountedPtr<ViscousBaseDomainBC>&          a_domainBC,
                     const RefCountedPtr<ViscousBaseEBBC>&              a_ebBC,
                     const bool&                                        a_hasMGObjects,
                     const IntVect&                                     a_ghostPhi,
                     const IntVect&                                     a_ghostRHS):
LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB>(false), // is time-independent
  m_eblgFine(a_eblgFine),
  m_eblg(a_eblg),
  m_eblgCoar(a_eblgCoar),
  m_eblgCoarMG(a_eblgCoarMG),
  m_alpha(a_alpha),
  m_beta(a_beta),
  m_alphaDiagWeight(),
  m_betaDiagWeight(),
  m_acoef(a_acoef),
  m_eta(a_eta),
  m_lambda(a_lambda),
  m_etaIrreg(a_etaIrreg),
  m_lambdaIrreg(a_lambdaIrreg),
  m_fastFR(),
  m_dx(a_dx),
  m_dxCoar(a_dxCoar),
  m_hasFine(a_eblgFine.isDefined()),
  m_hasCoar(a_eblgCoar.isDefined()),
  m_refToFine(a_refToFine),
  m_refToCoar(a_refToCoar),
  m_hasMGObjects(a_hasMGObjects),
  m_ghostPhi(a_ghostPhi),
  m_ghostRHS(a_ghostRHS),
  m_opEBStencil(),
  m_relCoef(),
  m_vofIterIrreg(),
  m_vofIterMulti(),
  m_vofIterDomLo(),
  m_vofIterDomHi(),
  m_interpWithCoarser(),
  m_ebAverage(),
  m_ebAverageMG(),
  m_ebInterp(),
  m_ebInterpMG(),
  m_domainBC(a_domainBC),
  m_ebBC(a_ebBC),
  m_colors()
{
  CH_TIME("nwoebvto::define");
  m_operatorForTimingOnly = false;
  EBArith::getMultiColors(m_colors);
  if (m_hasCoar)
    {
      CH_TIME("nwoebvto::coarsefine");
  
      ProblemDomain domainCoar = coarsen(m_eblg.getDomain(), m_refToCoar);
      m_interpWithCoarser = RefCountedPtr<NWOEBQuadCFInterp>(new NWOEBQuadCFInterp(m_eblg.getDBL(),  m_eblgCoar.getDBL(),
                                                                                   m_eblg.getEBISL(),m_eblgCoar.getEBISL(),
                                                                                   domainCoar, m_refToCoar, SpaceDim, m_dx,
                                                                                   m_ghostPhi, *m_eblg.getCFIVS()));
      //if this fails, then the AMR grids violate proper nesting.
      ProblemDomain domainCoarsenedFine;
      DisjointBoxLayout dblCoarsenedFine;

      int maxBoxSize = 32;
      bool dumbool;
      bool hasCoarser = EBAMRPoissonOp::getCoarserLayouts(dblCoarsenedFine,
                                                          domainCoarsenedFine,
                                                          m_eblg.getDBL(),
                                                          m_eblg.getEBISL(),
                                                          m_eblg.getDomain(),
                                                          m_refToCoar,
                                                          m_eblg.getEBIS(),
                                                          maxBoxSize, dumbool);

      //should follow from coarsenable
      if (hasCoarser)
        {
          //      EBLevelGrid eblgCoarsenedFine(dblCoarsenedFine, domainCoarsenedFine, 4, Chombo_EBIS::instance());
          EBLevelGrid eblgCoarsenedFine(dblCoarsenedFine, domainCoarsenedFine, 4, m_eblg.getEBIS() );
          m_ebInterp.define( m_eblg.getDBL(),     m_eblgCoar.getDBL(),
                             m_eblg.getEBISL(), m_eblgCoar.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, SpaceDim,
                             m_eblg.getEBIS(),     a_ghostPhi, true, true);
          m_ebAverage.define(m_eblg.getDBL(),   eblgCoarsenedFine.getDBL(),
                             m_eblg.getEBISL(), eblgCoarsenedFine.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, SpaceDim,
                             m_eblg.getEBIS(), a_ghostRHS);

        }
    }
  if (m_hasMGObjects)
    {
      int mgRef = 2;
      m_eblgCoarMG = a_eblgCoarMG;

      m_ebInterpMG.define( m_eblg.getDBL(),   m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain(), mgRef, SpaceDim,
                           m_eblg.getEBIS(),   a_ghostPhi, true, true);
      m_ebAverageMG.define(m_eblg.getDBL(),   m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain() , mgRef, SpaceDim,
                           m_eblg.getEBIS(),   a_ghostRHS);

    }
  /**/
  {
    CH_TIME("ivs calc");
    m_ivsIrregCCFlux.define(m_eblg.getDBL());
    DataIterator dit = m_eblg.getDBL().dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
      {
        IntVectSet& ivsIrreg = m_ivsIrregCCFlux[dit[mybox]];
        const ProblemDomain& domain = m_eblg.getDomain();
        Box ccFluxBox =  m_eblg.getDBL()[dit[mybox]];
        ccFluxBox.grow(1);
        ccFluxBox &= domain;

        ivsIrreg = m_eblg.getEBISL()[dit[mybox]].getIrregIVS(ccFluxBox);
      }
  }
  /**/

  m_exchangeCopier.define(m_eblg.getDBL(), m_eblg.getDBL(), m_ghostPhi,  true);
  
  defineStencils();
}
//-----------------------------------------------------------------------
NWOEBViscousTensorOp::
NWOEBViscousTensorOp(const EBLevelGrid &                                a_eblg,
                     const Real&                                        a_alpha,
                     const Real&                                        a_beta,
                     const RefCountedPtr<LevelData<EBCellFAB> >&        a_acoef,
                     const RefCountedPtr<LevelData<EBFluxFAB> >&        a_eta,
                     const RefCountedPtr<LevelData<EBFluxFAB> >&        a_lambda,
                     const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_etaIrreg,
                     const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_lambdaIrreg,
                     const Real&                                        a_dx,
                     const RefCountedPtr<ViscousBaseDomainBC>&          a_domainBC,
                     const RefCountedPtr<ViscousBaseEBBC>&              a_ebBC,
                     const IntVect&                                     a_ghostPhi,
                     const IntVect&                                     a_ghostRHS):
LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB>(false) // is time-independent
{
  CH_TIME("nwoebvto::defineforLoadBalance");
  m_operatorForTimingOnly = true;
  m_acoef       = a_acoef;
  m_eta         = a_eta;
  m_lambda      = a_lambda;
  m_etaIrreg    = a_etaIrreg;
  m_lambdaIrreg = a_lambdaIrreg;
  m_eblg        = a_eblg;
  m_ghostPhi    = a_ghostPhi;
  m_ghostRHS    = a_ghostRHS;
  m_domainBC    = a_domainBC;
  m_ebBC = a_ebBC;
  m_alpha = a_alpha;
  m_beta = a_beta;
  m_dx   =  a_dx;
  //none of these values matter as this definition is only for timing
  m_hasFine = false;
  m_hasCoar = false;
  m_refToFine= -1;
  m_refToCoar= -1;
  m_hasMGObjects = false;

  EBArith::getMultiColors(m_colors);

  m_exchangeCopier.define(m_eblg.getDBL(), m_eblg.getDBL(), m_ghostPhi,  true);

  defineStencilsForTimingOnly();
}
/**************/

void 
NWOEBViscousTensorOp::
averageCellToFace(EBFaceFAB           &      a_fluxData,
                  const EBCellFAB     &      a_cellData,
                  const Box           &      a_grid,
                  const EBISBox       &      a_ebisBox,
                  const ProblemDomain &      a_domain,
                  const DataIndex& a_dit,
                  int isrc, int idst, int inco,
                  bool a_interpolateToCentroid)
{
  CH_TIME("nwoebvto::averagecelltoface::averagecelltoface");
  CH_assert(!m_operatorForTimingOnly);

  int idir = a_fluxData.direction();
  EBFaceFAB cellCenteredFlux;
  Box ccFluxBox = a_grid;
  if (a_interpolateToCentroid)
    {
      ccFluxBox.grow(1);
      ccFluxBox &= a_domain;
    }

  cellCenteredFlux.define(a_ebisBox, ccFluxBox, idir, a_fluxData.nComp());
  cellCenteredFlux.setVal(0.);

  faceCenteredAverageCellsToFaces(cellCenteredFlux, a_cellData, ccFluxBox, a_ebisBox, a_domain, a_dit, isrc, idst, inco);
  //first copy (this does all the regular cells)
  Interval srcInt(isrc, isrc+inco-1);
  Interval dstInt(idst, idst+inco-1);
  a_fluxData.copy(a_grid, dstInt, a_grid,  cellCenteredFlux, srcInt);

  if (a_interpolateToCentroid)
    {
      CH_TIME("interpolation_to centroids");
      //if required, do the fancy interpolation to centroids.
      IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(a_grid);
      IntVectSet cfivs;
      FaceStop::WhichFaces stopCritGrid =  FaceStop::SurroundingWithBoundary;
      for (FaceIterator faceit(ivsIrreg, a_ebisBox.getEBGraph(), idir, stopCritGrid); faceit.ok(); ++faceit)
        {
          FaceStencil sten = EBArith::getInterpStencil(faceit(), cfivs, a_ebisBox, a_domain);
          for (int icomp = 0; icomp < inco; icomp++)
            {
              sten.setAllVariables(isrc+icomp);
              a_fluxData(faceit(), idst+icomp) = applyFaceStencil(sten, cellCenteredFlux, 0);
            }
        }
    }
}
//-----------------------------------------------------------------------
void 
NWOEBViscousTensorOp::
averageCellToFace(LevelData<EBFluxFAB>&         a_fluxData,
                  const LevelData<EBCellFAB>&   a_cellData,
                  const DisjointBoxLayout&      a_grids,
                  const EBISLayout&             a_ebisl,
                  const ProblemDomain&          a_domain,
                  int isrc, int idst, int inco,
                  bool a_interpolateToCentroid)
{

  CH_TIME("EBLDO::averagecelltoface::ldflux");
  CH_assert(!m_operatorForTimingOnly);
  DataIterator dit = a_grids.dataIterator();
  int nbox=dit.size();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //not thread safe
      //#pragma omp parallel for  
      for (int mybox=0;mybox<nbox; mybox++)
        {
          fillVelGhost(a_cellData[dit[mybox]], dit[mybox], false);
          averageCellToFace(a_fluxData[dit[mybox]][idir],
                            a_cellData[dit[mybox]],
                            a_grids[dit[mybox]],
                            a_ebisl[dit[mybox]],
                            a_domain, dit[mybox],
                            isrc, idst, inco,
                            a_interpolateToCentroid);
        }
    }

}
/****/
void 
NWOEBViscousTensorOp::
faceCenteredAverageCellsToFaces(EBFaceFAB           &      a_faceData,
                                const EBCellFAB     &      a_cellData,
                                const Box           &      ccFluxBox,
                                const EBISBox       &      a_ebisBox,
                                const ProblemDomain &      a_domain,
                                const DataIndex& a_dit,
                                int isrc, int idst, int inco)
{
  CH_TIME("nwoebvto::faceCenteredAverageCellsToFaces");
  CH_assert(!m_operatorForTimingOnly);
  //get the surrounding faces box
  Box faceBox = ccFluxBox;
  int idir = a_faceData.direction();
  faceBox.surroundingNodes(idir);
  //do regular cells in fortran for the faces in faceBox
  {
    CH_TIME("fortran calls");
    for (int icomp = 0; icomp < inco; icomp++)
      {
        int faceComp= idst + icomp;
        int cellComp= isrc + icomp;

        FORT_EBAVECELLTOFACE(CHF_FRA1(      a_faceData.getSingleValuedFAB(), faceComp),
                             CHF_CONST_FRA1(a_cellData.getSingleValuedFAB(), cellComp),
                             CHF_CONST_INT(idir),
                             CHF_BOX(faceBox));
      }
  
  }

  //fix up irregular faces and boundary faces
  IntVectSet& ivsIrreg = m_ivsIrregCCFlux[a_dit];
  {
    CH_TIME("faceiterator_foo");
    FaceStop::WhichFaces stopCritGrid =  FaceStop::SurroundingWithBoundary;
    for (FaceIterator faceit(ivsIrreg, a_ebisBox.getEBGraph(), idir, stopCritGrid); faceit.ok(); ++faceit)
      {
        const FaceIndex& face = faceit();
        if (!face.isBoundary())
          {
            for (int icomp = 0; icomp < inco; icomp++)
              {
                a_faceData(face, idst + icomp) = 0.5*(a_cellData(face.getVoF(Side::Hi), isrc + icomp) +
                                                      a_cellData(face.getVoF(Side::Lo), isrc + icomp));
              }
          }
        else
          {

            for (int icomp = 0; icomp < inco; icomp++)
              {
                VolIndex whichVoF;
                const Box& domainBox = a_domain.domainBox();
                if (domainBox.contains(face.gridIndex(Side::Lo)))
                  {
                    whichVoF = face.getVoF(Side::Lo);
                  }
                else if (domainBox.contains(face.gridIndex(Side::Hi)))
                  {
                    whichVoF = face.getVoF(Side::Hi);
                  }
                else
                  {
                    MayDay::Error("face and domain inconsistent:  logic bust in average cells to faces");
                  }

                a_faceData(face, idst + icomp) = a_cellData(whichVoF, isrc + icomp);
              }
          }
      }
  }

}

//-----------------------------------------------------------------------
void
NWOEBViscousTensorOp::
getKappaDivSigmaU(LevelData<EBCellFAB>& a_divSigmaU,
                  const LevelData<EBCellFAB>& a_velocity,
                  const LevelData<EBCellFAB>* a_veloCoar,
                  int a_level)
{
  CH_TIME("nwoebvto::getkappadivsigmau");
  CH_assert(!m_operatorForTimingOnly);
  LevelData<EBCellFAB>& velcast =   (LevelData<EBCellFAB>&) a_velocity;
  if (a_level > 0)
    {
      CH_TIME("nwoebvto::getkappadivsigmau::cfinterp");
      m_interpWithCoarser->coarseFineInterp(velcast, *a_veloCoar, 0, 0, SpaceDim);
    }
  velcast.exchange();

  EBFluxFactory fact(m_eblg.getEBISL());
  LevelData<EBFluxFAB>       velFluxCenter(m_eblg.getDBL(), SpaceDim, IntVect::Unit, fact);
  LevelData<EBFluxFAB>     velFluxCentroid(m_eblg.getDBL(), SpaceDim, IntVect::Unit, fact);
  LevelData<EBFluxFAB>       sigma(m_eblg.getDBL(), SpaceDim, IntVect::Zero, fact);
  LevelData<EBFluxFAB> velDotSigma(m_eblg.getDBL(),        1,  IntVect::Zero, fact);

  averageCellToFace(velFluxCenter, a_velocity,
                    m_eblg.getDBL(), m_eblg.getEBISL(),
                    m_eblg.getDomain(), 0, 0, SpaceDim, true);

  velFluxCenter.exchange();
  Interval inter(0, SpaceDim-1);
  velFluxCenter.copyTo(inter, velFluxCentroid, inter);
  {
    CH_TIME("nwoebvto::getkappadivsigmau::interpolatetocentroid");
    EBArith::interpolateFluxToCentroids(velFluxCentroid,
                                        velFluxCenter,
                                        m_eblg.getDBL(),
                                        m_eblg.getEBISL(),
                                        m_eblg.getDomain());
  }

  //so now we have velocity at faces.
  //now we get the viscous flux at faces, dot it with the velocity and take the divergence.
  int ibox=0;
  {
    CH_TIME("nwoebvto::getkappadivsigmau::getfluxbit");

    DataIterator dit = m_eblg.getDBL().dataIterator();
    int nbox=dit.size();
    //not thread safe
    //#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
      {
        EBFluxFAB& thisFlux = sigma[dit[mybox]];
        thisFlux.setVal(0.);
        getFlux(thisFlux, a_velocity, m_eblg.getDBL()[dit[mybox]], dit[mybox], 1.0);

        ibox++;
      }
  }

  getVelDotSigma(velDotSigma, velFluxCentroid, sigma);

  //now we can take the divergence of vel dot sigma and add it to the output
  CH_assert(velDotSigma.ghostVect() == IntVect::Zero);
  CH_assert(a_divSigmaU.ghostVect() == m_ghostPhi);
  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
#pragma omp parallel for
  for (int mybox=0;mybox<nbox; mybox++)
    {
      BaseIVFAB<Real> dummy; //no boundary flux here
      m_divergenceStencil[dit[mybox]]->divergence(a_divSigmaU[dit[mybox]], velDotSigma[dit[mybox]], dummy, 0, false);
    }
}

//-----------------------------------------------------------------------
void
NWOEBViscousTensorOp::
getCCSigma(LevelData<EBCellFAB>      & a_sigma,
           const LevelData<EBCellFAB>& a_gradU,
           const LevelData<EBCellFAB>& a_etaCell,
           const LevelData<EBCellFAB>& a_lambdaCell
           )
{
  CH_TIME("nwoebvto::getCCSigma");
  CH_assert(!m_operatorForTimingOnly);
  //Now we loop through all the boxes and do the appropriate multiplication of
  //cofficients with gradients to get sigma
  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
  {
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
      {
        const Box& grid = m_eblg.getDBL()[dit[mybox]];
        a_sigma[dit[mybox]].setVal(0); //needed because this is done through increment
        //first get the eta*(grad u + grad u^T)
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            for (int jdir = 0; jdir < SpaceDim; jdir++)
              {
                int gradcomp = TensorCFInterp::gradIndex(idir, jdir);
                int trancomp = TensorCFInterp::gradIndex(jdir, idir);
                FORT_INCRCCSIGMAETA(CHF_FRA1(        a_sigma[dit[mybox]].getSingleValuedFAB(), gradcomp),
                                    CHF_CONST_FRA(   a_gradU[dit[mybox]].getSingleValuedFAB()),
                                    CHF_CONST_FRA1(a_etaCell[dit[mybox]].getSingleValuedFAB(), 0),
                                    CHF_BOX(grid),
                                    CHF_INT(gradcomp),
                                    CHF_INT(trancomp));
                if (idir == jdir)
                  {
                    //now add the lambda*divU bits
                    for (int divDir = 0; divDir < SpaceDim; divDir++)
                      {
                        int diagcomp = TensorCFInterp::gradIndex(divDir, divDir);
                        FORT_INCRCCSIGMALAMBDA(CHF_FRA1(           a_sigma[dit[mybox]].getSingleValuedFAB(), gradcomp),
                                               CHF_CONST_FRA(      a_gradU[dit[mybox]].getSingleValuedFAB()),
                                               CHF_CONST_FRA1(a_lambdaCell[dit[mybox]].getSingleValuedFAB(), 0),
                                               CHF_BOX(grid),
                                               CHF_INT(diagcomp));
                      }
                  }
              }
            
          }
        for (m_vofIterIrreg[dit[mybox]].reset(); m_vofIterIrreg[dit[mybox]].ok(); ++m_vofIterIrreg[dit[mybox]])
          {
            const VolIndex& vof = m_vofIterIrreg[dit[mybox]]();
            for (int idir = 0; idir < SpaceDim; idir++)
              {
                for (int jdir = 0; jdir < SpaceDim; jdir++)
                  {
                    Real irregSigma = 0;
                    
                    int gradcomp = TensorCFInterp::gradIndex(idir, jdir);
                    int trancomp = TensorCFInterp::gradIndex(jdir, idir);
                    // put in eta(grad u + grad u ^T)
                    irregSigma += a_etaCell[dit[mybox]](vof, 0)*(a_gradU[dit[mybox]](vof, gradcomp) + a_gradU[dit[mybox]](vof, trancomp));
                    
                    if (idir == jdir)
                      {
                        //now add the lambda*divU bits
                        for (int divDir = 0; divDir < SpaceDim; divDir++)
                          {
                            int diagcomp = TensorCFInterp::gradIndex(divDir, divDir);
                            irregSigma += a_lambdaCell[dit[mybox]](vof, 0)*a_gradU[dit[mybox]](vof, diagcomp);
                          }
                      }
                    a_sigma[dit[mybox]](vof, gradcomp) = irregSigma;
                  }
              }
          }
      }
  } //end pragma
}
//----------------------------------------
//simple averaging
void
NWOEBViscousTensorOp::
averageToCells(LevelData<EBCellFAB>&      a_cellCoef,
               const LevelData<EBFluxFAB>& a_faceCoef,
               const LevelData<BaseIVFAB<Real> >& a_irregCoef)
{
  CH_TIME("nwoebvto::averageToCells");
  CH_assert(!m_operatorForTimingOnly);
  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
  {
    //    int id=omp_get_thread_num();
    // pout() << "my thread "<< id << endl;
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
      {
        a_cellCoef[dit[mybox]].setVal(0.); //necessary because we are doing increments
        const Box& grid = m_eblg.getDBL()[dit[mybox]];
        //number of faces on a regular cell
        Real nfacesper = 2*SpaceDim;
        EBCellFAB&       cellCo = a_cellCoef[dit[mybox]];
        for (int idir= 0; idir < SpaceDim; idir++)
          {
            const EBFaceFAB& faceCo = a_faceCoef[dit[mybox]][idir];
            FORT_INCRFACETOCELLAVERAGE(CHF_FRA1(      cellCo.getSingleValuedFAB(), 0),
                                       CHF_CONST_FRA1(faceCo.getSingleValuedFAB(), 0),
                                       CHF_INT(idir),
                                       CHF_BOX(grid),
                                       CHF_REAL(nfacesper));
          }
        for (m_vofIterIrreg[dit[mybox]].reset(); m_vofIterIrreg[dit[mybox]].ok(); ++m_vofIterIrreg[dit[mybox]])
          {
            const VolIndex& vof = m_vofIterIrreg[dit[mybox]]();
            Real irregCoef = a_irregCoef[dit[mybox]](vof, 0);
            int nfaces = 1; //1 is the irregular face
            for (int idir = 0; idir < SpaceDim ; idir++)
              {
                for (SideIterator sit; sit.ok();++sit)
                  {
                    Vector<FaceIndex> faces = m_eblg.getEBISL()[dit[mybox]].getFaces(vof, idir, sit());
                    nfaces += faces.size();
                    for (int iface = 0; iface < faces.size(); iface++)
                      {
                        irregCoef += a_faceCoef[dit[mybox]][idir](faces[iface], 0);
                      }
                  }
              }
            irregCoef /= nfaces;   //always >=1 by construction
          }
      }
  }// end pragma
}
void
NWOEBViscousTensorOp::
getCellCenteredCoefficients(LevelData<EBCellFAB>&    a_etaCell,
                            LevelData<EBCellFAB>& a_lambdaCell)
{
  CH_TIME("nwoebvto::getCellCenteredCoefficients");
  CH_assert(!m_operatorForTimingOnly);
  averageToCells(a_etaCell, *m_eta, *m_etaIrreg);
  averageToCells(a_lambdaCell, *m_lambda, *m_lambdaIrreg);
}
//----------------------------------------
void
NWOEBViscousTensorOp::
getShearStressDotGradU(LevelData<EBCellFAB>      & a_shearStressDotGradUU,
                       const LevelData<EBCellFAB> & a_gradU,
                       int a_level)
{

  CH_TIME("nwoebvto::getShearStressDotGradU");
  CH_assert(!m_operatorForTimingOnly);
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB>   lambdaCell(m_eblg.getDBL(), SpaceDim         , IntVect::Zero, fact);
  LevelData<EBCellFAB>      etaCell(m_eblg.getDBL(), SpaceDim         , IntVect::Zero, fact);
  LevelData<EBCellFAB>        sigma(m_eblg.getDBL(), SpaceDim*SpaceDim, IntVect::Zero, fact);
  //get coeffcients centered at cells
  getCellCenteredCoefficients(etaCell, lambdaCell);

  //use these to geta  cell-centered sigma
  getCCSigma(sigma, a_gradU, etaCell, lambdaCell);

  //take the dot product of the gradient and the stress
  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
  {
    //    int id=omp_get_thread_num();
    //pout() << "my thread "<< id << endl;
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
      {
        const Box& grid = m_eblg.getDBL()[dit[mybox]];
        //      const EBISBox& ebisBox = m_eblg.getEBISL()[dit[mybox]];
        a_shearStressDotGradUU[dit[mybox]].setVal(0.);
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            for (int jdir = 0; jdir < SpaceDim; jdir++)
              {
                int gradComp = TensorCFInterp::gradIndex(idir, jdir);
                FORT_INCRPOINTDOTPROD(CHF_FRA1(a_shearStressDotGradUU[dit[mybox]].getSingleValuedFAB(), 0),
                                      CHF_CONST_FRA1(       a_gradU[dit[mybox]].getSingleValuedFAB(), gradComp),
                                      CHF_CONST_FRA1(         sigma[dit[mybox]].getSingleValuedFAB(), gradComp),
                                      CHF_BOX(grid));
              }
          }
        for (m_vofIterIrreg[dit[mybox]].reset(); m_vofIterIrreg[dit[mybox]].ok(); ++m_vofIterIrreg[dit[mybox]])
          {
            const VolIndex& vof = m_vofIterIrreg[dit[mybox]]();
            
            Real dotprod = 0;
            for (int idir = 0; idir < SpaceDim; idir++)
              {
                for (int jdir = 0; jdir < SpaceDim; jdir++)
                  {
                    int gradcomp = TensorCFInterp::gradIndex(idir, jdir);
                    dotprod += sigma[dit[mybox]](vof, gradcomp)*a_gradU[dit[mybox]](vof, gradcomp);
                  }
              }
            a_shearStressDotGradUU[dit[mybox]](vof, 0) = dotprod;
          }
      }
  } // end pragma
}

//--------------------------------------------------------------------
void
NWOEBViscousTensorOp::
getVelDotSigma(LevelData<EBFluxFAB>      & a_velDotSigma,
               const LevelData<EBFluxFAB>& a_vel,
               const LevelData<EBFluxFAB>& a_sigma)
{
  CH_TIME("nwoebvto::getveldotsigma");
  CH_assert(!m_operatorForTimingOnly);
  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
  
  //not thread safe
  //#pragma omp parallel for
  for (int mybox=0;mybox<nbox; mybox++)
    {
      const Box& grid = m_eblg.getDBL()[dit[mybox]];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit[mybox]];
      IntVectSet ivsCell = ebisBox.getIrregIVS(grid);
      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          EBFaceFAB&       dotFace = a_velDotSigma[dit[mybox]][faceDir];
          const EBFaceFAB& velFace =         a_vel[dit[mybox]][faceDir];
          const EBFaceFAB& sigFace =       a_sigma[dit[mybox]][faceDir];
          
          Box faceBox = grid;
          faceBox.surroundingNodes(faceDir);
          FORT_VELDOTSIGMA(CHF_FRA(      dotFace.getSingleValuedFAB()),
                           CHF_CONST_FRA(velFace.getSingleValuedFAB()),
                           CHF_CONST_FRA(sigFace.getSingleValuedFAB()),
                           CHF_BOX(faceBox));
          
          for (FaceIterator faceit(ivsCell, ebisBox.getEBGraph(), faceDir,stopCrit);
               faceit.ok(); ++faceit)
            {
              const FaceIndex& face = faceit();
              Real dotval = 0;
              for (int velDir = 0; velDir < SpaceDim; velDir++)
                {
                  Real sigmaval = a_sigma[dit[mybox]][faceDir](face, velDir);
                  Real   velval =   a_vel[dit[mybox]][faceDir](face, velDir);
                  dotval += sigmaval*velval;
                }
              a_velDotSigma[dit[mybox]][faceDir](faceit(), 0) = dotval;
            }
          stopCrit = FaceStop::AllBoundaryOnly;
          IntVectSet ivsBox(grid);
          //set boundary faces to zero
          for (FaceIterator faceit(ivsCell, ebisBox.getEBGraph(), faceDir,stopCrit);
               faceit.ok(); ++faceit)
            {
              a_velDotSigma[dit[mybox]][faceDir](faceit(), 0) = 0.;
            }
        }
    }
}

//-----------------------------------------------------------------------
Real
NWOEBViscousTensorOp::
getSafety()
{
  Real safety = 0.25;
  return safety;
}
//-----------------------------------------------------------------------
void
NWOEBViscousTensorOp::
setTime(Real a_oldTime, Real a_mu, Real a_dt)
{
}
//-----------------------------------------------------------------------
void
NWOEBViscousTensorOp::
setAlphaAndBeta(const Real& a_alpha,
                const Real& a_beta)
{
  CH_TIME("nwoebvto::setAlphaAndBeta");
  CH_assert(!m_operatorForTimingOnly);
  m_alpha = a_alpha;
  m_beta  = a_beta;
  calculateRelaxationCoefficient();
  calculateAlphaWeight(); //have to do this because a coef has probably been changed under us.
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
NWOEBViscousTensorOp::
calculateAlphaWeight()
{
  CH_TIME("nwoebvto::calculateAlphaWeight");
  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
#pragma omp parallel for
  for (int mybox=0;mybox<nbox; mybox++)
    {
      VoFIterator& vofit = m_vofIterIrreg[dit[mybox]];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();
          Real volFrac = m_eblg.getEBISL()[dit[mybox]].volFrac(VoF);
          Real alphaWeight = (*m_acoef)[dit[mybox]](VoF, 0);
          alphaWeight *= volFrac;

          for (int ivar = 0; ivar < SpaceDim; ivar++)
            {
              m_alphaDiagWeight[dit[mybox]](VoF, ivar) = alphaWeight;
            }
        }
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
NWOEBViscousTensorOp::
calculateRelaxationCoefficient()
{
  CH_TIME("nwoebvto::calcRelCoef");

  //define regular relaxation coefficent
  //define regular relaxation coefficent
  Real safety = getSafety();
  int ncomp = SpaceDim;
  // Not threadsafe yet
  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
#pragma omp parallel for
  for (int mybox=0;mybox<nbox; mybox++)
    {
      const Box& grid = m_eblg.getDBL().get(dit[mybox]);
      const EBCellFAB& acofab = (*m_acoef)[dit[mybox]];
      //initialize lambda = alpha*acoef
      m_relCoef[dit[mybox]].setVal(0.);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_relCoef[dit[mybox]].plus(acofab, 0, idir, 1);
        }
      m_relCoef[dit[mybox]]*= m_alpha;

      BaseFab<Real>& regRel =   m_relCoef[dit[mybox]].getSingleValuedFAB();

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox& regEta = (FArrayBox &)((*m_eta)   [dit[mybox]][idir].getSingleValuedFAB());
          FArrayBox& regLam = (FArrayBox &)((*m_lambda)[dit[mybox]][idir].getSingleValuedFAB());
          FORT_DECRINVRELCOEFVTOP(CHF_FRA(regRel),
                                  CHF_FRA(regEta),
                                  CHF_FRA(regLam),
                                  CHF_CONST_REAL(m_beta),
                                  CHF_BOX(grid),
                                  CHF_REAL(m_dx),
                                  CHF_INT(idir),
                                  CHF_INT(ncomp));
        }

      //now invert so lambda = stable lambda for variable coef lapl
      //(according to phil, this is the correct lambda)
      FORT_INVERTLAMBDAVTOP(CHF_FRA(regRel),
                            CHF_REAL(safety),
                            CHF_BOX(grid),
                            CHF_INT(ncomp));

      VoFIterator& vofit = m_vofIterIrreg[dit[mybox]];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();
          Vector<Real> diagWeight(ncomp, 0);
          Real maxDiag = 0;
          for (int ivar = 0; ivar < ncomp; ivar++)
            {
              Real alphaWeight = m_alphaDiagWeight[dit[mybox]](VoF, ivar);
              Real  betaWeight =  m_betaDiagWeight[dit[mybox]](VoF, ivar);
              alphaWeight *= m_alpha;
              betaWeight  *= m_beta;

              diagWeight[ivar] = alphaWeight + betaWeight;
              maxDiag =Max(Abs(diagWeight[ivar]), maxDiag);
            }
          for (int ivar = 0; ivar < ncomp; ivar++)
            {
              if (Abs(diagWeight[ivar]) > 1.0e-12)
                {
                  m_relCoef[dit[mybox]](VoF, ivar) = safety/diagWeight[ivar];
                }
              else
                {
                  m_relCoef[dit[mybox]](VoF, ivar) = 0.0;
                }
            }
        }
    }
}
//-----------------------------------------------------------------------

/*****/
void
NWOEBViscousTensorOp::
defineStencils()
{
  CH_TIME("nwoebvto::defineStencils");
  CH_assert(!m_operatorForTimingOnly);

  m_divergenceStencil.define(m_eblg.getDBL());
  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
  //not thread safe
  //#pragma omp parallel for
  for (int mybox=0;mybox<nbox; mybox++)
    {

      const DataIndex& datind = dit[mybox];
      BaseIVFAB<Real> dumbiv; //no irregular flux for this bit
      Box smallBox = m_eblg.getDBL()[datind];
      Box grownBox = smallBox;
      grownBox.grow(m_ghostPhi);
      const EBISBox& ebisBox = m_eblg.getEBISL()[datind];
      EBCellFAB dumEBCF(ebisBox, grownBox,1);
      EBFluxFAB dumEBFF(ebisBox, smallBox,1);
      m_divergenceStencil[datind] = RefCountedPtr<DivergenceStencil>( new DivergenceStencil(dumEBCF, dumEBFF, dumbiv, smallBox, ebisBox, m_dx*RealVect::Unit, false));
    }
  m_vofIterIrreg.define(m_eblg.getDBL());
  m_vofIterMulti.define(m_eblg.getDBL());
  m_alphaDiagWeight.define(  m_eblg.getDBL());
  m_betaDiagWeight.define(   m_eblg.getDBL());
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  m_relCoef.define(m_eblg.getDBL(), SpaceDim, IntVect::Zero, ebcellfact);

  Box sideBoxLo[SpaceDim];
  Box sideBoxHi[SpaceDim];
  for (int ivar = 0; ivar < SpaceDim; ivar++)
    {
      int idir = ivar;
      Box domainBox = m_eblg.getDomain().domainBox();
      sideBoxLo[idir] = adjCellLo(domainBox, idir, 1);
      sideBoxLo[idir].shift(idir,  1);
      sideBoxHi[idir] = adjCellHi(domainBox, idir, 1);
      sideBoxHi[idir].shift(idir, -1);
      m_opEBStencil[ivar].define(m_eblg.getDBL());
      m_vofIterDomLo[ivar].define( m_eblg.getDBL()); // vofiterator cache for domain lo
      m_vofIterDomHi[ivar].define( m_eblg.getDBL()); // vofiterator cache for domain hi
    }
  Real fakeBeta = 1;
  m_domainBC->setCoef(m_eblg,   fakeBeta,      m_eta,      m_lambda);
  m_ebBC->setCoef(    m_eblg,   fakeBeta,      m_etaIrreg, m_lambdaIrreg, m_eta, m_lambda);
  m_ebBC->define(    *m_eblg.getCFIVS(), 1.0);
  //first define the iterators
  // DataIterator dit = m_eblg.getDBL().dataIterator();
  // int nbox=dit.size();
  //not thread safe
  //#pragma omp parallel for
  EBCellFactory fact(m_eblg.getEBISL());              
  LevelData<EBCellFAB> phiProxy(m_eblg.getDBL(), SpaceDim, m_ghostPhi, fact);
  LevelData<EBCellFAB> rhsProxy(m_eblg.getDBL(), SpaceDim, m_ghostRHS, fact);
  for (int mybox=0;mybox<nbox; mybox++)
    {
      const DataIndex& datind = dit[mybox];
      const EBISBox& ebisBox = m_eblg.getEBISL()[datind];
      const Box&     grid = m_eblg.getDBL().get(datind);
      //need to grow the irregular set by one near multivalued cells
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);
      IntVectSet ivsMulti = ebisBox.getMultiCells(grid);
      /**/
      //    This nonsense is supposed to be handled in EBIS but is not
      ivsMulti.grow(1);
      ivsMulti &= grid;

      IntVectSet ivsComplement(grid);
      ivsComplement -= ivsMulti;
      ivsComplement -= ivsIrreg;
      //ivscomplement now contains the complement of the cells we need;

      IntVectSet ivsStenc(grid);
      ivsStenc -= ivsComplement;
      /**/

      m_alphaDiagWeight[datind].define(ivsStenc, ebisBox.getEBGraph(), SpaceDim);
      m_betaDiagWeight [datind].define(ivsStenc, ebisBox.getEBGraph(), SpaceDim);
      m_vofIterIrreg   [datind].define(ivsStenc, ebisBox.getEBGraph()   );
      m_vofIterMulti   [datind].define(ivsMulti, ebisBox.getEBGraph()   );

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          IntVectSet loIrreg = ivsStenc;
          IntVectSet hiIrreg = ivsStenc;
          loIrreg &= sideBoxLo[idir];
          hiIrreg &= sideBoxHi[idir];
          m_vofIterDomLo[idir][datind].define(loIrreg,ebisBox.getEBGraph());
          m_vofIterDomHi[idir][datind].define(hiIrreg,ebisBox.getEBGraph());
        }

      //now define the stencils for each variable
      for (int ivar = 0; ivar < SpaceDim; ivar++)
        {
          LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(ivar);
          DataIterator dit = m_eblg.getDBL().dataIterator();

          const EBISBox& ebisBox = m_eblg.getEBISL()[datind];
          VoFIterator &  vofit = m_vofIterIrreg[datind];
          const Vector<VolIndex>& vofvec = vofit.getVector();
          // cast from VolIndex to BaseIndex
          Vector<RefCountedPtr<BaseIndex> >    dstVoF(vofvec.size());
          Vector<RefCountedPtr<BaseStencil> > stencil(vofvec.size());
          for(int ivec = 0; ivec < vofvec.size(); ivec++)
            {
              const VolIndex& vof = vofvec[ivec];
              VoFStencil vofsten;
              getVoFStencil(vofsten, vof, datind, ivar);
              if (fluxStencil != NULL)
                {
                  BaseIVFAB<VoFStencil>& stenFAB = (*fluxStencil)[datind];
                  //only real irregular cells are going to have EB fluxes
                  //but we are working on a bigger set
                  if (stenFAB.getIVS().contains(vof.gridIndex()))
                    {
                      VoFStencil fluxStencilPt = stenFAB(vof, 0);
                      Real bndryArea = ebisBox.bndryArea(vof);
                      //this might need to get multiplied by the -normal[idir]
                      Real factor = bndryArea/m_dx;
                      fluxStencilPt *= factor;
                      vofsten += fluxStencilPt;
                    }
                }
              dstVoF[ivec]  = RefCountedPtr<BaseIndex  >(new  VolIndex(vof));
              stencil[ivec] = RefCountedPtr<BaseStencil>(new VoFStencil(vofsten));

              Real diagWeight = EBArith::getDiagWeight(vofsten, vof, ivar);
              m_betaDiagWeight[datind](vof, ivar) = diagWeight;
            }
          m_opEBStencil[ivar][datind] = RefCountedPtr<VCAggStencil>
            (new VCAggStencil(dstVoF, stencil, phiProxy[datind], rhsProxy[datind], m_relCoef[datind], m_alphaDiagWeight[datind],SpaceDim));

        }
    }
  calculateAlphaWeight();
  calculateRelaxationCoefficient();
}
/*****/
void
NWOEBViscousTensorOp::
defineStencilsForTimingOnly()
{
 CH_TIME("nwoebvto::defineStencilsForTimingONly");
  CH_assert(m_operatorForTimingOnly);

  m_vofIterIrreg.define(m_eblg.getDBL());
  m_vofIterMulti.define(m_eblg.getDBL());
  m_alphaDiagWeight.define(  m_eblg.getDBL());
  m_betaDiagWeight.define(   m_eblg.getDBL());
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  m_relCoef.define(m_eblg.getDBL(), SpaceDim, IntVect::Zero, ebcellfact);

  m_domainBC->setCoef(m_eblg,   m_beta,      m_eta,      m_lambda);
  m_ebBC->setCoef(    m_eblg,   m_beta,      m_etaIrreg, m_lambdaIrreg, m_eta, m_lambda);
  m_ebBC->define(    *m_eblg.getCFIVS(), 1.0);
  for (int ivar = 0; ivar < SpaceDim; ivar++)
    {
      m_opEBStencil[ivar].define(m_eblg.getDBL());
    }
  m_vofIterIrreg.define(m_eblg.getDBL());

  EBCellFactory fact(m_eblg.getEBISL());              
  LevelData<EBCellFAB> phiProxy(m_eblg.getDBL(), SpaceDim, m_ghostPhi, fact);
  LevelData<EBCellFAB> rhsProxy(m_eblg.getDBL(), SpaceDim, m_ghostRHS, fact);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const DataIndex& datind = dit();
      const EBISBox& ebisBox = m_eblg.getEBISL()[datind];
      const Box&     grid = m_eblg.getDBL().get(datind);
      //need to grow the irregular set by one near multivalued cells
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);
      IntVectSet ivsMulti = ebisBox.getMultiCells(grid);
      /**/
      //    This nonsense is supposed to be handled in EBIS but is not
      ivsMulti.grow(1);
      ivsMulti &= grid;

      IntVectSet ivsComplement(grid);
      ivsComplement -= ivsMulti;
      ivsComplement -= ivsIrreg;
      //ivscomplement now contains the complement of the cells we need;

      IntVectSet ivsStenc(grid);
      ivsStenc -= ivsComplement;
      /**/

      m_vofIterIrreg[datind].define(ivsStenc, ebisBox.getEBGraph());
      m_alphaDiagWeight[datind].define(ivsStenc, ebisBox.getEBGraph(), SpaceDim);
      m_betaDiagWeight [datind].define(ivsStenc, ebisBox.getEBGraph(), SpaceDim);

      //now define the stencils for each variable
      for (int ivar = 0; ivar < SpaceDim; ivar++)
        {
          LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(ivar);
          DataIterator dit = m_eblg.getDBL().dataIterator();

          const EBISBox& ebisBox = m_eblg.getEBISL()[datind];
          VoFIterator& vofit  = m_vofIterIrreg[datind];
          const Vector<VolIndex>& vofvec = vofit.getVector();
          // cast from VolIndex to BaseIndex
          Vector<RefCountedPtr<BaseIndex> >    dstVoF(vofvec.size());
          Vector<RefCountedPtr<BaseStencil> > stencil(vofvec.size());
          for(int ivec = 0; ivec < vofvec.size(); ivec++)
            {
              const VolIndex& vof = vofvec[ivec];
              VoFStencil vofsten;
              getVoFStencil(vofsten, vof, datind, ivar);
              if (fluxStencil != NULL)
                {
                  BaseIVFAB<VoFStencil>& stenFAB = (*fluxStencil)[datind];
                  //only real irregular cells are going to have EB fluxes
                  //but we are working on a bigger set
                  if (stenFAB.getIVS().contains(vof.gridIndex()))
                    {
                      VoFStencil fluxStencilPt = stenFAB(vof, 0);
                      Real bndryArea = ebisBox.bndryArea(vof);
                      //this might need to get multiplied by the -normal[idir]
                      Real factor = bndryArea/m_dx;
                      fluxStencilPt *= factor;
                      vofsten += fluxStencilPt;
                    }
                }
              dstVoF[ivec]  = RefCountedPtr<BaseIndex  >(new  VolIndex(vof));
              stencil[ivec] = RefCountedPtr<BaseStencil>(new VoFStencil(vofsten));

              Real diagWeight = EBArith::getDiagWeight(vofsten, vof, ivar);
              m_betaDiagWeight[datind](vof, ivar) = diagWeight;
            }
          m_opEBStencil[ivar][datind] = RefCountedPtr<VCAggStencil>
            (new VCAggStencil(dstVoF, stencil, phiProxy[datind], rhsProxy[datind], m_relCoef[datind], m_alphaDiagWeight[datind],SpaceDim));

        }
    }
  calculateAlphaWeight();
  calculateRelaxationCoefficient();
}
/* generate vof stencil as a divergence of flux stencils */
/***/
void
NWOEBViscousTensorOp::
getVoFStencil(VoFStencil&      a_vofStencil,
              const VolIndex&  a_vof,
              const DataIndex& a_dit,
              int             a_ivar)
{
  CH_TIME("nwoebvto::getVoFStencil");
  const EBISBox& ebisBox = m_eblg.getEBISL()[a_dit];
  a_vofStencil.clear();

  //doing the diagonal thing and everything in line
  Real volFrac = m_eblg.getEBISL()[a_dit].volFrac(a_vof);
  Real alphaWeight =    (*m_acoef)[a_dit](a_vof, 0);
  alphaWeight *= volFrac;
  a_vofStencil.add(a_vof, alphaWeight, a_ivar);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          Vector<FaceIndex> faces = ebisBox.getFaces(a_vof, idir, sit());
          for (int iface = 0; iface < faces.size(); iface++)
            {
              VoFStencil fluxStencil;
              getFluxStencil(fluxStencil, m_eta, m_lambda, m_dx, m_eblg, faces[iface], a_dit, a_ivar);
              Real areaFrac = ebisBox.areaFrac(faces[iface]);
              fluxStencil *= Real(isign)*areaFrac/m_dx;
              a_vofStencil += fluxStencil;
            }
        }
    }
}
/***/
/// stencil for flux computation.   the truly ugly part of this computation
/// beta and eta are multiplied in here
/****/
void
NWOEBViscousTensorOp::
getFluxStencil(VoFStencil                                & a_fluxStencil,
               const RefCountedPtr<LevelData<EBFluxFAB> >& a_eta,
               const RefCountedPtr<LevelData<EBFluxFAB> >& a_lambda,
               const Real                                & a_dx,
               const EBLevelGrid                         & a_eblg,
               const FaceIndex                           & a_face,
               const DataIndex                           & a_dit,
               int a_ivar)
{
  //need to do this by interpolating to centroids
  //so get the stencil at each face center and add with
  //interpolation weights
  FaceStencil interpSten = EBArith::getInterpStencil(a_face,
                                                     (*a_eblg.getCFIVS())[a_dit],
                                                     a_eblg.getEBISL()[a_dit],
                                                     a_eblg.getDomain());

  a_fluxStencil.clear();
  for (int isten = 0; isten < interpSten.size(); isten++)
    {
      const FaceIndex& face = interpSten.face(isten);
      const Real&    weight = interpSten.weight(isten);
      VoFStencil faceCentSten;
      getFaceCenteredFluxStencil(faceCentSten, a_eta, a_lambda, a_dx, a_eblg, face, a_dit, a_ivar);
      faceCentSten *= weight;
      a_fluxStencil += faceCentSten;
    }
  //debug
  //a_fluxStencil.clear();
  //if (a_face.direction() == 1)  getFaceCenteredFluxStencil(a_fluxStencil, a_face, a_dit, a_ivar);
  //end debug

}
void
NWOEBViscousTensorOp::
getFaceCenteredFluxStencil(VoFStencil                                & a_fluxStencil,
                           const RefCountedPtr<LevelData<EBFluxFAB> >& a_eta,
                           const RefCountedPtr<LevelData<EBFluxFAB> >& a_lambda,
                           const Real                                & a_dx,
                           const EBLevelGrid                         & a_eblg,
                           const FaceIndex                           & a_face,
                           const DataIndex                           & a_dit,
                           int a_ivar)
{
  int faceDir= a_face.direction();
  VoFStencil  gBTranSten, gBNormSten;
  //the flux is lambda I diverge(B) + eta(gradB + gradB tran)
  //so the ivar flux for the faceDir direction is
  // lambda*delta(ivar, faceDir)*divergence(B) + eta*(partial(B_ivar, faceDir) + partial(B_faceDir, ivar))

  getGradientStencil(gBNormSten, a_ivar , faceDir, a_face, a_dit, a_dx, a_eblg);
  getGradientStencil(gBTranSten, faceDir, a_ivar , a_face, a_dit, a_dx, a_eblg);
  Real etaFace = (*a_eta)[a_dit][faceDir](a_face, 0);
  gBNormSten *=    etaFace;
  gBTranSten *=    etaFace;

  a_fluxStencil += gBNormSten;
  a_fluxStencil += gBTranSten;
  if (a_ivar == faceDir)
    {
      Real lambdaFace = (*a_lambda)[a_dit][faceDir](a_face, 0);
      VoFStencil divergeSten;
      getDivergenceStencil(divergeSten, a_face, a_dit, a_dx, a_eblg);
      divergeSten   *= lambdaFace;
      a_fluxStencil += divergeSten;
    }
}
/***/
void
NWOEBViscousTensorOp::
getGradientStencil(VoFStencil&  a_gradStencil,
                   int a_ivar,
                   int a_diffDir,
                   const FaceIndex& a_face,
                   const DataIndex& a_dit)
{
  getGradientStencil(a_gradStencil, a_ivar, a_diffDir, a_face, a_dit, m_dx, m_eblg);
}

void
NWOEBViscousTensorOp::
getGradientStencil(VoFStencil&  a_gradStencil,
                   int a_ivar,
                   int a_diffDir,
                   const FaceIndex& a_face,
                   const DataIndex& a_dit,
                   const Real     & a_dx,
                   const EBLevelGrid& a_eblg)
{
  CH_TIME("nwoebvto::getGradientStencil");
  a_gradStencil.clear();
  if ((a_face.direction() == a_diffDir) && (!a_face.isBoundary()))
    {
      a_gradStencil.add(a_face.getVoF(Side::Hi),  1.0/a_dx, a_ivar);
      a_gradStencil.add(a_face.getVoF(Side::Lo), -1.0/a_dx, a_ivar);
    }
  else
    {
      int numGrad = 0;
      for (SideIterator sit; sit.ok(); ++sit)
        {
          const VolIndex& vofSide = a_face.getVoF(sit());
          if (a_eblg.getDomain().contains(vofSide.gridIndex()))
            {
              VoFStencil stenSide;
              IntVectSet& cfivs = (*a_eblg.getCFIVS())[a_dit];
              EBArith::getFirstDerivStencil(stenSide, vofSide,
                                            a_eblg.getEBISL()[a_dit],
                                            a_diffDir, a_dx, &cfivs, a_ivar);
              a_gradStencil += stenSide;
              numGrad++;
            }
        }
      if (numGrad > 1)
        {
          a_gradStencil *= 1.0/Real(numGrad);
        }
    }
}

/***/
void
NWOEBViscousTensorOp::
getDivergenceStencil(VoFStencil&      a_divStencil,
                     const FaceIndex& a_face,
                     const DataIndex& a_dit)
{
  getDivergenceStencil(a_divStencil, a_face, a_dit, m_dx, m_eblg);
}
void
NWOEBViscousTensorOp::
getDivergenceStencil(VoFStencil&      a_divStencil,
                     const FaceIndex& a_face,
                     const DataIndex& a_dit,
                     const Real     & a_dx,
                     const EBLevelGrid& a_eblg)
{
  a_divStencil.clear();
  for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
    {
      VoFStencil diffStencilDir;
      //difference direction and variable are the same
      //because that is the definition of the divergence.
      getGradientStencil(diffStencilDir, diffDir, diffDir, a_face, a_dit, a_dx, a_eblg);
      a_divStencil += diffStencilDir;
    }
}
/***/
NWOEBViscousTensorOp::
~NWOEBViscousTensorOp()
{
}
/***/
void
NWOEBViscousTensorOp::
AMRResidualNC(LevelData<EBCellFAB>&       a_residual,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_rhs,
              bool a_homogeneousBC,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("nwoebvto::amrresNC");
  CH_assert(!m_operatorForTimingOnly);
  //dummy. there is no coarse when this is called
  CH_assert(a_residual.ghostVect() == m_ghostRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostRHS);
  LevelData<EBCellFAB> phiC;
  AMRResidual(a_residual, a_phiFine, a_phi, phiC, a_rhs, a_homogeneousBC, a_finerOp);
}
/***/
void
NWOEBViscousTensorOp::
AMROperatorNC(LevelData<EBCellFAB>&       a_LofPhi,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              bool a_homogeneousBC,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("nwoebvto::amropNC");
  CH_assert(!m_operatorForTimingOnly);
  //dummy. there is no coarse when this is called
  CH_assert(a_LofPhi.ghostVect() == m_ghostRHS);
  LevelData<EBCellFAB> phiC;
  AMROperator(a_LofPhi, a_phiFine, a_phi, phiC,
              a_homogeneousBC, a_finerOp);
}
/*****/
void
NWOEBViscousTensorOp::
residual(LevelData<EBCellFAB>&       a_residual,
         const LevelData<EBCellFAB>& a_phi,
         const LevelData<EBCellFAB>& a_rhs,
         bool                        a_homogeneousBC)
{
  //this is a multigrid operator so only homogeneous CF BC
  //and null coar level
  CH_TIME("nwoebvto::residual");
  CH_assert(!m_operatorForTimingOnly);
  CH_assert(a_residual.ghostVect() == m_ghostRHS);
  CH_assert(a_phi.ghostVect() == m_ghostPhi);
  applyOp(a_residual,a_phi, a_homogeneousBC);
  //  incr(a_residual, a_rhs, -1.0);
  //  scale(a_residual, -1.0);
  axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
}

/*****/
void
NWOEBViscousTensorOp::
preCond(LevelData<EBCellFAB>&       a_phi,
        const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("nwoebvto::precond");
  CH_assert(!m_operatorForTimingOnly);
  relax(a_phi, a_rhs, 1);
}



/*****/
void
NWOEBViscousTensorOp::
create(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("nwoebvto::create");
  CH_assert(!m_operatorForTimingOnly);
  int ncomp = a_rhs.nComp();
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  a_lhs.define(m_eblg.getDBL(), ncomp, a_rhs.ghostVect(), ebcellfact);
  assign(a_lhs, a_rhs);
}

/*****/
void
NWOEBViscousTensorOp::
createCoarsened(LevelData<EBCellFAB>&       a_lhs,
                const LevelData<EBCellFAB>& a_rhs,
                const int&                  a_refRat)
{
  CH_TIME("nwoebvto::createCoar");
  int ncomp = a_rhs.nComp();
  CH_assert(!m_operatorForTimingOnly);
  IntVect ghostVect = a_rhs.ghostVect();

  CH_assert(m_eblg.getDBL().coarsenable(a_refRat));

  //fill ebislayout
  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, m_eblg.getDBL(), a_refRat);

  EBISLayout ebislCoarsenedFine;
  IntVect ghostVec = a_rhs.ghostVect();
  //const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ProblemDomain coarDom = coarsen(m_eblg.getDomain(), a_refRat);
  m_eblg.getEBIS()->fillEBISLayout(ebislCoarsenedFine, dblCoarsenedFine, coarDom , ghostVec[0]);
  if (m_refToCoar > 1)
    {
      ebislCoarsenedFine.setMaxRefinementRatio(m_refToCoar, m_eblg.getEBIS());
    }

  //create coarsened data
  EBCellFactory ebcellfactCoarsenedFine(ebislCoarsenedFine);
  a_lhs.define(dblCoarsenedFine, ncomp,ghostVec, ebcellfactCoarsenedFine);
}

/*****/
Real
NWOEBViscousTensorOp::
AMRNorm(const LevelData<EBCellFAB>& a_coarResid,
        const LevelData<EBCellFAB>& a_fineResid,
        const int& a_refRat,
        const int& a_ord)

{
  MayDay::Error("never called");
  return -1.0;
}

/*****/
void
NWOEBViscousTensorOp::
assign(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("nwoebvto::assign");
  EBLevelDataOps::assign(a_lhs,a_rhs);
}

/*****/
Real
NWOEBViscousTensorOp::
dotProduct(const LevelData<EBCellFAB>& a_1,
           const LevelData<EBCellFAB>& a_2)
{
  CH_TIME("nwoebvto::dotprod");
  ProblemDomain domain;
  Real volume;

  return EBLevelDataOps::kappaDotProduct(volume,a_1,a_2,EBLEVELDATAOPS_ALLVOFS,domain);
}

/*****/
void
NWOEBViscousTensorOp::
incr(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     Real                        a_scale)
{
  CH_TIME("nwoebvto::incr");
  EBLevelDataOps::incr(a_lhs,a_x,a_scale);
}

/*****/
void
NWOEBViscousTensorOp::
axby(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     const LevelData<EBCellFAB>& a_y,
     Real                        a_a,
     Real                        a_b)
{
  CH_TIME("nwoebvto::axby");
  EBLevelDataOps::axby(a_lhs,a_x,a_y,a_a,a_b);
}

/*****/
void
NWOEBViscousTensorOp::
scale(LevelData<EBCellFAB>& a_lhs,
      const Real&           a_scale)
{
  CH_TIME("nwoebvto::scale");
  EBLevelDataOps::scale(a_lhs,a_scale);
}

Real 
NWOEBViscousTensorOp::
norm(const LevelData<EBCellFAB>& a_rhs,
     int                         a_ord)
{
  CH_TIMERS("nwoebvto::norm");
  CH_TIMER("mpi_allreduce",t1);
  CH_assert(!m_operatorForTimingOnly);

  Real maxNorm = 0.0;

  maxNorm = localMaxNorm(a_rhs);

  CH_START(t1);

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&maxNorm, &tmp, 1, MPI_CH_REAL,
                             MPI_MAX, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communcation error on norm");
    }
  maxNorm = tmp;
#endif

  CH_STOP(t1);

  return maxNorm;
}
///
Real 
NWOEBViscousTensorOp::
localMaxNorm(const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("NWOEBViscousTensorOp::localMaxNorm");
  return  NWOEBViscousTensorOp::staticMaxNorm(a_rhs, m_eblg);
}
///
Real 
NWOEBViscousTensorOp::
staticMaxNorm(const LevelData<EBCellFAB>& a_rhs, const EBLevelGrid& a_eblg)
{
  CH_TIME("nwoebvto::staticmaxnorm");

  Real maxNorm = 0.0;

  for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
    {
      int iRegIrregCovered;
      const BaseFab<char>& maskFAB = a_eblg.getEBISL()[dit()].getEBGraph().getMask(iRegIrregCovered);
      for(int icomp = 0; icomp < SpaceDim; icomp++)
        {
          if (iRegIrregCovered != -1)//not all covered
            {
              if (iRegIrregCovered == 0)//has irreg
                {
                  const Box& box = a_eblg.getDBL().get(dit());
                  //             const EBISBox& ebisBox = a_eblg.getEBISL()[dit()];
                  const BaseFab<Real>& rhsFAB = (a_rhs[dit()]).getSingleValuedFAB();
                  FORT_MAXNORMMASK(CHF_REAL(maxNorm),
                                   CHF_CONST_FRA1(rhsFAB,icomp),
                                   CHF_BOX(box),
                                   CHF_CONST_FBA1(maskFAB,0));

                  //CP: this portion between the stars is the new faster code:
                  //****************************
                  int srccomp = icomp;
                  int ncomp   = 1;

                  const MiniIVFAB<Real>& irrBFAB1 = a_rhs[dit()].getMultiValuedFAB();
                  const Real* r = irrBFAB1.dataPtr(srccomp);
                  int nvof    = irrBFAB1.numVoFs();
                  for (int i=0; i<nvof*ncomp; i++)
                    {
                      maxNorm = Max(maxNorm, Abs(r[i]));
                    }
                }
              else//all reg
                {
                  const Box& box = a_eblg.getDBL().get(dit());
                  const BaseFab<Real>& rhsFAB = (a_rhs[dit()]).getSingleValuedFAB();
                  FORT_MAXNORM(CHF_REAL(maxNorm),
                               CHF_CONST_FRA1(rhsFAB,icomp),
                               CHF_BOX(box));
                }
            }
        }
    }
  return maxNorm;
}

/*****/
void
NWOEBViscousTensorOp::
setToZero(LevelData<EBCellFAB>& a_lhs)
{
  CH_TIME("nwoebvto::setToZero");
  EBLevelDataOps::setToZero(a_lhs);
}

/*****/
void
NWOEBViscousTensorOp::
setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value)
{
  CH_TIME("nwoebvto::setVal");
  EBLevelDataOps::setVal(a_lhs, a_value);
}

void
NWOEBViscousTensorOp::
createCoarser(LevelData<EBCellFAB>&       a_coar,
              const LevelData<EBCellFAB>& a_fine,
              bool                        a_ghosted)
{
  CH_TIME("nwoebvto::createCoarser");
  CH_assert(!m_operatorForTimingOnly);
  const DisjointBoxLayout& dbl = m_eblgCoarMG.getDBL();
  ProblemDomain coarDom = coarsen(m_eblg.getDomain(), 2);

  int nghost = a_fine.ghostVect()[0];
  EBISLayout coarEBISL;

  //  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  const EBIndexSpace* const ebisPtr = m_eblg.getEBIS();
  ebisPtr->fillEBISLayout(coarEBISL,
                          dbl, coarDom, nghost);

  EBCellFactory ebcellfact(coarEBISL);
  a_coar.define(dbl, SpaceDim,a_fine.ghostVect(),ebcellfact);
}
/*****/
void
NWOEBViscousTensorOp::
relaxTimed(LevelData<EBCellFAB>&            a_phi,
           const LevelData<EBCellFAB>&      a_rhs,
           TimedDataIterator         &      a_dit,
           int                              a_iterations)
{
  CH_TIME("nwoebvto::relax");
  CH_assert(m_operatorForTimingOnly);
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  // do first red, then black passes
  for (int whichIter =0; whichIter < a_iterations; whichIter++)
    {
      for (int icolor = 0; icolor < m_colors.size(); icolor++)
        {
          if ((icolor == 0) || (!s_doLazyRelax))
            {
              CH_TIME("ghostfill and homogcfinterp");
              fillVelGhostTimed(a_phi, true, a_dit);
              //no cf interp here--this is just for load balancing
            }
          gsrbColorTimed(a_phi, a_rhs, m_colors[icolor], a_dit);
        }
    }

}

/*****/
void
NWOEBViscousTensorOp::
relax(LevelData<EBCellFAB>&       a_phi,
      const LevelData<EBCellFAB>& a_rhs,
      int                         a_iterations)
{
  CH_TIME("nwoebvto::relax");
  CH_assert(!m_operatorForTimingOnly);
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  //change to false to get a lot of printouts
  static bool printedStuff = true;

  if(!printedStuff)
    {
      pout() << "rhs: " << endl;;
      LevelData<EBCellFAB>* rhs = (LevelData<EBCellFAB>*)(&a_rhs);
      printMaxMinLDCell(rhs);
    }
  // do first red, then black passes
  for (int whichIter =0; whichIter < a_iterations; whichIter++)
    {
      for (int icolor = 0; icolor < m_colors.size(); icolor++)
        {
          if ((icolor == 0) || (!s_doLazyRelax))
            {
              CH_TIME("ghostfill and homogcfinterp");
              fillVelGhost(a_phi, true);
              homogeneousCFInterp(a_phi);
            }
          if(!printedStuff)
            {
              pout() << "iter = "<< whichIter << ", icolor = " << icolor << endl;
              pout() << "phi: " ;
              printMaxMinLDCell(&a_phi);
            }
          gsrbColor(a_phi, a_rhs, m_colors[icolor]);
        }
    }

  //MayDay::Abort("leaving after first relax");
  if(!printedStuff)
    {
      pout() << "rhs: " << endl;;
      LevelData<EBCellFAB>* rhs = (LevelData<EBCellFAB>*)(&a_rhs);
      printMaxMinLDCell(rhs);
      pout() << "phi coming out of relax: "  << endl;
      printMaxMinLDCell(&a_phi);
    }
  printedStuff = true;
}

void
NWOEBViscousTensorOp::
homogeneousCFInterp(LevelData<EBCellFAB>&   a_phif)
{
  CH_TIME("nwoebvto::homog_cfinterp");
  if (m_hasCoar)
    {
      m_interpWithCoarser->coarseFineInterpH(a_phif,  0, 0, SpaceDim);
    }
}
/*****/
void
NWOEBViscousTensorOp::
gsrbColor(LevelData<EBCellFAB>&       a_phi,
          const LevelData<EBCellFAB>& a_rhs,
          const IntVect&              a_color)
{
  CH_TIME("nwoebvto::gsrbColor");
  CH_assert(!m_operatorForTimingOnly);
  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
  {
    // int id=omp_get_thread_num();
    // pout() << "my thread "<< id << endl;
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
      {
        //first do the the regular stuff
        Box dblBox  = dbl.get(dit[mybox]);
        const DataIndex datInd = dit[mybox];
        BaseFab<Real>&       regPhi =     a_phi[datInd].getSingleValuedFAB();
        const BaseFab<Real>& regRhs =     a_rhs[datInd].getSingleValuedFAB();
        const BaseFab<Real>& regRel = m_relCoef[datInd].getSingleValuedFAB();


        //assumes ghost cells already filled

        const BaseFab<Real>& acofab  =(*m_acoef)[datInd].getSingleValuedFAB();
        const EBFluxFAB& eta  =     (*m_eta)[datInd];
        const EBFluxFAB& lam  =  (*m_lambda)[datInd];

        Vector<const BaseFab<Real>* > etaside(3, &(eta[0].getSingleValuedFAB()));
        Vector<const BaseFab<Real>* > lamside(3, &(lam[0].getSingleValuedFAB()));
        for(int idir = 0; idir < SpaceDim; idir++)
          {
            etaside[idir] = &(eta[idir].getSingleValuedFAB());
            lamside[idir] = &(lam[idir].getSingleValuedFAB());
          }

        IntVect loIV = dblBox.smallEnd();
        IntVect hiIV = dblBox.bigEnd();
        
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (loIV[idir] % 2 != a_color[idir])
              {
                loIV[idir]++;
              }
          }
        //only has to cached once (though any of the aggstencils could do it)
        m_opEBStencil[0][datInd]->cachePhi(a_phi[datInd]);

        if (loIV <= hiIV)
          {
            Box coloredBox(loIV, hiIV);
            FORT_GSRBNWOEBVTOP(CHF_FRA(regPhi),
                               CHF_CONST_FRA(regRhs),
                               CHF_CONST_FRA(regRel),
                               CHF_CONST_FRA1(acofab,0),
                               CHF_CONST_FRA1((*etaside[0]), 0),
                               CHF_CONST_FRA1((*etaside[1]), 0),
                               CHF_CONST_FRA1((*etaside[2]), 0),
                               CHF_CONST_FRA1((*lamside[0]), 0),
                               CHF_CONST_FRA1((*lamside[1]), 0),
                               CHF_CONST_FRA1((*lamside[2]), 0),
                               CHF_CONST_REAL(m_dx),
                               CHF_CONST_REAL(m_alpha),
                               CHF_CONST_REAL(m_beta),
                               CHF_BOX(coloredBox));
          }


        m_opEBStencil[0][datInd]->uncachePhi(a_phi[datInd]);

        for(int ivar = 0; ivar < SpaceDim; ivar++)
          {
            m_opEBStencil[ivar][datInd]->relax(a_phi[datInd], a_rhs[datInd], 
                                               m_relCoef[datInd], m_alphaDiagWeight[datInd],
                                               m_alpha, m_beta, ivar, a_color);
          }
      }
  }// end pragma
}

/*****/
void
NWOEBViscousTensorOp::
gsrbColorTimed(LevelData<EBCellFAB>       & a_phi,
               const LevelData<EBCellFAB> & a_rhs,
               const IntVect              & a_color,
               TimedDataIterator          & a_dit)
{
  CH_TIME("nwoebvto::gsrbColorTimed");
  CH_assert(m_operatorForTimingOnly);
  DisjointBoxLayout dbl = m_eblg.getDBL();
  for(a_dit.reset(); a_dit.ok(); ++a_dit)
    {
      //first do the the regular stuff
      const DataIndex datInd = a_dit();
      Box dblBox  = dbl.get(datInd);
      BaseFab<Real>&       regPhi =     a_phi[datInd].getSingleValuedFAB();
      const BaseFab<Real>& regRhs =     a_rhs[datInd].getSingleValuedFAB();
      const BaseFab<Real>& regRel = m_relCoef[datInd].getSingleValuedFAB();


      //assumes ghost cells already filled

      const BaseFab<Real>& acofab  =(*m_acoef)[datInd].getSingleValuedFAB();
      const EBFluxFAB& eta  =     (*m_eta)[datInd];
      const EBFluxFAB& lam  =  (*m_lambda)[datInd];

      Vector<const BaseFab<Real>* > etaside(3, &(eta[0].getSingleValuedFAB()));
      Vector<const BaseFab<Real>* > lamside(3, &(lam[0].getSingleValuedFAB()));
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          etaside[idir] = &(eta[idir].getSingleValuedFAB());
          lamside[idir] = &(lam[idir].getSingleValuedFAB());
        }

      IntVect loIV = dblBox.smallEnd();
      IntVect hiIV = dblBox.bigEnd();
        
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (loIV[idir] % 2 != a_color[idir])
            {
              loIV[idir]++;
            }
        }
      //only has to cached once (though any of the aggstencils could do it)
      m_opEBStencil[0][datInd]->cachePhi(a_phi[datInd]);

      if (loIV <= hiIV)
        {
          Box coloredBox(loIV, hiIV);
          FORT_GSRBNWOEBVTOP(CHF_FRA(regPhi),
                             CHF_CONST_FRA(regRhs),
                             CHF_CONST_FRA(regRel),
                             CHF_CONST_FRA1(acofab,0),
                             CHF_CONST_FRA1((*etaside[0]), 0),
                             CHF_CONST_FRA1((*etaside[1]), 0),
                             CHF_CONST_FRA1((*etaside[2]), 0),
                             CHF_CONST_FRA1((*lamside[0]), 0),
                             CHF_CONST_FRA1((*lamside[1]), 0),
                             CHF_CONST_FRA1((*lamside[2]), 0),
                             CHF_CONST_REAL(m_dx),
                             CHF_CONST_REAL(m_alpha),
                             CHF_CONST_REAL(m_beta),
                             CHF_BOX(coloredBox));
        }


      m_opEBStencil[0][datInd]->uncachePhi(a_phi[datInd]);

      for(int ivar = 0; ivar < SpaceDim; ivar++)
        {
          m_opEBStencil[ivar][datInd]->relax(a_phi[datInd], a_rhs[datInd], 
                                             m_relCoef[datInd], m_alphaDiagWeight[datInd],
                                             m_alpha, m_beta, ivar, a_color);
        }
    }
}
/*****/
void
NWOEBViscousTensorOp::
restrictResidual(LevelData<EBCellFAB>&       a_resCoar,
                 LevelData<EBCellFAB>&       a_phi,
                 const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("nwoebvto::restrictRes");
  CH_assert(!m_operatorForTimingOnly);
  LevelData<EBCellFAB> res;
  bool homogeneous = true;
  CH_assert(a_resCoar.ghostVect() == m_ghostRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostRHS);
  create(res, a_rhs);

  // Get the residual on the fine grid
  residual(res,a_phi,a_rhs,homogeneous);

  // now use our nifty averaging operator
  Interval variables(0, SpaceDim-1);
  if (m_hasMGObjects)
    {
      m_ebAverageMG.average(a_resCoar, res, variables);
    }
  else
    {
      m_ebAverage.average(a_resCoar, res, variables);
    }
}

/*****/
void
NWOEBViscousTensorOp::
prolongIncrement(LevelData<EBCellFAB>&       a_phi,
                 const LevelData<EBCellFAB>& a_cor)
{
  CH_TIME("nwoebvto::prolongInc");
  CH_assert(!m_operatorForTimingOnly);
  CH_assert(a_phi.ghostVect() == m_ghostPhi);
  CH_assert(a_cor.ghostVect() == m_ghostPhi);
  Interval vars(0, SpaceDim-1);
  if (m_hasMGObjects)
    {
      m_ebInterpMG.pwlInterp(a_phi, a_cor, vars);
    }
  else
    {
      m_ebInterp.pwlInterp(a_phi, a_cor, vars);
    }
}

/*****/
int
NWOEBViscousTensorOp::
refToCoarser()
{
  return m_refToCoar;
}

/*****/
int
NWOEBViscousTensorOp::
refToFiner()
{
  return m_refToFine;
}

/*****/
void
NWOEBViscousTensorOp::
AMRResidual(LevelData<EBCellFAB>& a_residual,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoarse,
            const LevelData<EBCellFAB>& a_rhs,
            bool a_homogeneousBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("nwoebvto::amrRes");

  CH_assert(!m_operatorForTimingOnly);
  this->cfinterp(a_phi, a_phiCoarse);

  applyOp(a_residual, a_phi, a_homogeneousBC);

  if (a_finerOp != NULL)
    {
      reflux(a_phiFine, a_phi,  a_residual, a_finerOp);
    }
  axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
  //  incr(a_residual, a_rhs, -1.0);
  //  scale(a_residual, -1.0);
}

/*****/
void
NWOEBViscousTensorOp::
AMRResidualNF(LevelData<EBCellFAB>& a_residual,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_phiCoarse,
              const LevelData<EBCellFAB>& a_rhs,
              bool a_homogeneousBC)
{
  CH_TIME("nwoebvto::amrResNF");
  CH_assert(!m_operatorForTimingOnly);
  this->cfinterp(a_phi, a_phiCoarse);
  this->residual(a_residual, a_phi, a_rhs, a_homogeneousBC ); //apply boundary conditions
}

/*****/
void
NWOEBViscousTensorOp::
reflux(const LevelData<EBCellFAB>&        a_phiFine,
       const LevelData<EBCellFAB>&        a_phi,
       LevelData<EBCellFAB>&              a_residual,
       AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("nwoebvto::reflux");
  CH_assert(!m_operatorForTimingOnly);
  Interval interv(0,SpaceDim-1);

  int ncomp = SpaceDim;
  if (!m_fastFR.isDefined())
    {
      m_fastFR.define(m_eblgFine, m_eblg, m_refToFine, ncomp, s_forceNoEBCF);
    }

  m_fastFR.setToZero();

  incrementFRCoar(m_fastFR, a_phiFine, a_phi);

  incrementFRFine(m_fastFR, a_phiFine, a_phi, a_finerOp);


  Real scale = 1.0/m_dx;
  m_fastFR.reflux(a_residual, interv, scale);
}

/****/
void
NWOEBViscousTensorOp::
incrementFRCoar(EBFastFR& a_fluxReg,
                const LevelData<EBCellFAB>& a_phiFine,
                const LevelData<EBCellFAB>& a_phi)
{
  CH_TIME("nwoebvto::incrFRC");
  CH_assert(!m_operatorForTimingOnly);
  int ncomp = SpaceDim;
  Interval interv(0,SpaceDim-1);
  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
  for (int mybox=0;mybox<nbox; mybox++)
    {
      const EBCellFAB& coarfab = a_phi[dit[mybox]];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit[mybox]];
      const Box&  box = m_eblg.getDBL().get(dit[mybox]);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //no boundary faces here.

          Box ghostedBox = box;
          ghostedBox.grow(1);
          ghostedBox.grow(idir,-1);
          ghostedBox &= m_eblg.getDomain();

          EBFaceFAB coarflux(ebisBox, ghostedBox, idir, ncomp);
          if (s_forceNoEBCF)
            {
              Box faceBox = surroundingNodes(box, idir);
              FArrayBox& regFlux = (FArrayBox &)     coarflux.getSingleValuedFAB();
              const FArrayBox& regPhi = (const FArrayBox &) coarfab.getSingleValuedFAB();
              getFlux(regFlux, regPhi,  faceBox, idir, dit[mybox]);
            }
          else
            {
              getFlux(coarflux, coarfab, ghostedBox, box, m_eblg.getDomain(), ebisBox, m_dx, dit[mybox], idir);
            }


          Real scale = 1.0; //beta and bcoef already in flux
          for (SideIterator sit; sit.ok(); ++sit)
            {
              a_fluxReg.incrementCoarseBoth(coarflux, scale, dit[mybox], interv, idir, sit());
            }
        }
    }
}
/***/
void
NWOEBViscousTensorOp::
incrementFRFine(EBFastFR& a_fluxReg,
                const LevelData<EBCellFAB>& a_phiFine,
                const LevelData<EBCellFAB>& a_phi,
                AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("nwoebvto::incrFRF");
  CH_assert(!m_operatorForTimingOnly);
  int ncomp = SpaceDim;
  Interval interv(0,SpaceDim-1);
  NWOEBViscousTensorOp& finerEBAMROp = (NWOEBViscousTensorOp& )(*a_finerOp);
  //ghost cells of phiFine need to be filled
  LevelData<EBCellFAB>& phiFine = (LevelData<EBCellFAB>&) a_phiFine;


  finerEBAMROp.m_interpWithCoarser->coarseFineInterp(phiFine, a_phi, 0, 0, SpaceDim);
  phiFine.exchange(finerEBAMROp.m_exchangeCopier);
  Real dxFine = m_dx/m_refToFine;
  DataIterator ditf = a_phiFine.dataIterator();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const Box&     boxFine = m_eblgFine.getDBL().get(ditf());
      const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[ditf()];
      const EBCellFAB& phiFine = a_phiFine[ditf()];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); sit.next())
            {
              //              Box fabBox = adjCellBox(boxFine, idir,        sit(), 1);
              //              fabBox.shift(idir, -sign(sit()));
              Box fabBox = boxFine;

              Box ghostedBox = fabBox;
              ghostedBox.grow(1);
              ghostedBox.grow(idir,-1);
              ghostedBox &= m_eblgFine.getDomain();

              EBFaceFAB fluxFine(ebisBoxFine, ghostedBox, idir, ncomp);

              if (s_forceNoEBCF)
                {
                  Box faceBox = surroundingNodes(fabBox, idir);
                  FArrayBox& regFlux = (FArrayBox &)     fluxFine.getSingleValuedFAB();
                  const FArrayBox& regPhi = (const FArrayBox &) phiFine.getSingleValuedFAB();
                  finerEBAMROp.getFlux(regFlux, regPhi,  faceBox, idir, ditf());
                }
              else
                {
                  finerEBAMROp.getFlux(fluxFine, phiFine, ghostedBox, fabBox, m_eblgFine.getDomain(),
                                       ebisBoxFine, dxFine, ditf(), idir);
                }

              Real scale = 1.0; //beta and bcoef already in flux

              a_fluxReg.incrementFineBoth(fluxFine, scale, ditf(), interv, idir, sit());
            }
        }
    }
}
void
NWOEBViscousTensorOp::
getFlux(EBFluxFAB&                    a_flux,
        const LevelData<EBCellFAB>&   a_data,
        const Box&                    a_grid,
        const DataIndex&              a_dit,
        Real                          a_scale)
{
  CH_TIME("nwoebvto::getFlux1");
  CH_assert(!m_operatorForTimingOnly);
  a_flux.define(m_eblg.getEBISL()[a_dit], a_grid, SpaceDim);
  a_flux.setVal(0.);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box ghostedBox = a_grid;
      ghostedBox.grow(1);
      ghostedBox.grow(idir,-1);
      ghostedBox &= m_eblg.getDomain();
      fillVelGhost(a_data[a_dit], a_dit, false);
      getFlux(a_flux[idir],
              a_data[a_dit],
              ghostedBox, a_grid,
              m_eblg.getDomain(),
              m_eblg.getEBISL()[a_dit],
              m_dx, a_dit, idir);
    }
}

void
NWOEBViscousTensorOp::
fillVelGhost(const LevelData<EBCellFAB>& a_phi, bool a_homog) const
{
  CH_TIME("nwoebvto::fillVelGhostLD");
  DataIterator dit = m_eblg.getDBL().dataIterator();

  //exchange must happen first to get corners correct
  LevelData<EBCellFAB>& phi = (LevelData<EBCellFAB>&)(a_phi);
  phi.exchange(m_exchangeCopier);

  int nbox=dit.size();
#pragma omp parallel for
  for (int mybox=0;mybox<nbox; mybox++)
    {
      fillVelGhost(a_phi[dit[mybox]], dit[mybox], a_homog);
    }
  
}
void
NWOEBViscousTensorOp::
fillVelGhostTimed(const LevelData<EBCellFAB>& a_phi, bool a_homog, TimedDataIterator& a_dit) const
{
  CH_TIME("nwoebvto::fillVelGhostTimed");

  for (a_dit.reset(); a_dit.ok(); ++a_dit)
    {
      fillVelGhost(a_phi[a_dit()], a_dit(), a_homog);
    }
  //took out exchange because it confuses the issue of timings
}
void
NWOEBViscousTensorOp::
fillVelGhost(const EBCellFAB& a_phi, const DataIndex& a_datInd, bool a_homog) const
{
  CH_TIME("nwoebvto::fillghostfab");
  EBCellFAB& phi = (EBCellFAB&) a_phi;
  ViscousBaseDomainBC* viscbc = dynamic_cast<ViscousBaseDomainBC*>(&(*m_domainBC));
  Box grid = m_eblg.getDBL()[a_datInd];
  Box domBox = m_eblg.getDomain().domainBox();
  if(viscbc == NULL)
    {
      MayDay::Error("dynamic cast failed");
    }
  if (!s_turnOffBCs)
    {
      FArrayBox& fab = phi.getFArrayBox();
      viscbc->fillVelGhost(fab, grid, domBox, m_dx, a_homog);
    }
  else
    {
      Box valid = m_eblg.getDBL()[a_datInd];
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          for(SideIterator sit; sit.ok(); ++sit)
            {
              FArrayBox& fab = phi.getFArrayBox();
              ExtrapolateBC(fab, valid,  m_dx, idir, sit());
            }
          //grow so that we hit corners
          valid.grow(idir, 1);
        }
    }
}
void
NWOEBViscousTensorOp::
getFlux(EBFaceFAB&                    a_fluxCentroid,
        const EBCellFAB&              a_phi,
        const Box&                    a_ghostedBox,
        const Box&                    a_fabBox,
        const ProblemDomain&          a_domain,
        const EBISBox&                a_ebisBox,
        const Real&                   a_dx,
        const DataIndex&              a_datInd,
        const int&                    a_idir)
{
  CH_TIME("nwoebvto::getFlux2.1");
  CH_assert(!m_operatorForTimingOnly);

  const FArrayBox& regPhi  = (const FArrayBox &)            a_phi.getSingleValuedFAB();

  //has some extra cells so...
  a_fluxCentroid.setVal(0.);


  int ncomp = a_phi.nComp();
  CH_assert(ncomp == a_fluxCentroid.nComp());

  Box faceBox = surroundingNodes(a_ghostedBox, a_idir);
  EBFaceFAB fluxCenter(a_ebisBox, a_ghostedBox, a_idir,SpaceDim);
  FArrayBox&       regFlux = (      FArrayBox &)   fluxCenter.getSingleValuedFAB();
  fillVelGhost(a_phi, a_datInd, false);
  getFlux(regFlux, regPhi,  faceBox, a_idir, a_datInd);

  a_fluxCentroid.copy(fluxCenter);

  IntVectSet ivsCell = a_ebisBox.getIrregIVS(a_ghostedBox);
  if (!ivsCell.isEmpty())
    {
      FaceStop::WhichFaces stopCrit = FaceStop::SurroundingNoBoundary;

      for (FaceIterator faceit(ivsCell, a_ebisBox.getEBGraph(), a_idir,stopCrit);
           faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          for (int ivar = 0; ivar < SpaceDim; ivar++)
            {
              VoFStencil fluxStencil;
              getFluxStencil(fluxStencil, m_eta, m_lambda, m_dx, m_eblg, face, a_datInd, ivar);
              //ivar is ignored in this call--vars live in the stencil.
              Real fluxVal =  applyVoFStencil(fluxStencil, a_phi, ivar);
              fluxCenter(faceit(), ivar) = fluxVal;
            }
        }

      //interpolate from face centers to face centroids
      Box cellBox = a_fluxCentroid.getCellRegion();
      EBArith::interpolateFluxToCentroids(a_fluxCentroid,
                                          fluxCenter,
                                          a_fabBox,
                                          a_ebisBox,
                                          a_domain,
                                          a_idir);
    }

  a_fluxCentroid *= m_beta;
}

/*****/
void
NWOEBViscousTensorOp::
AMROperator(LevelData<EBCellFAB>& a_LofPhi,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoarse,
            bool a_homogeneousBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("nwoebvto::amrOp");
  CH_assert(!m_operatorForTimingOnly);
  cfinterp(a_phi, a_phiCoarse);

  applyOp(a_LofPhi, a_phi, a_homogeneousBC);
  if (a_finerOp != NULL)
    {
      reflux(a_phiFine, a_phi,  a_LofPhi, a_finerOp);
    }
}

/*****/
void
NWOEBViscousTensorOp::
AMROperatorNF(LevelData<EBCellFAB>& a_LofPhi,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_phiCoarse,
              bool a_homogeneousBC)
{
  CH_TIME("nwoebvto::amrOpNF");
  CH_assert(!m_operatorForTimingOnly);

  cfinterp(a_phi, a_phiCoarse);

  //apply boundary conditions in applyOp
  this->applyOp(a_LofPhi, a_phi, a_homogeneousBC );
}
/*****/
void
NWOEBViscousTensorOp::
AMRRestrict(LevelData<EBCellFAB>&       a_resCoar,
            const LevelData<EBCellFAB>& a_residual,
            const LevelData<EBCellFAB>& a_correction,
            const LevelData<EBCellFAB>& a_coarseCorr, 
            bool a_skip_res )
{  
  CH_assert(!a_skip_res);
  CH_assert(!m_operatorForTimingOnly);
  CH_TIME("nwoebvto::amrRestrict");
  LevelData<EBCellFAB> res;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_residual.ghostVect();

  res.define(m_eblg.getDBL(), SpaceDim, ghostVec, ebcellfactTL);
  EBLevelDataOps::setVal(res, 0.0);

  cfinterp(a_correction, a_coarseCorr);
  //API says that we must average(a_residual - L(correction, coarCorrection))
  applyOp(res, a_correction,  true);
  incr(res, a_residual, -1.0);
  scale(res,-1.0);

  //use our nifty averaging operator
  Interval variables(0, SpaceDim-1);
  m_ebAverage.average(a_resCoar, res, variables);
}
void
NWOEBViscousTensorOp::
AMRProlong(LevelData<EBCellFAB>&       a_correction,
           const LevelData<EBCellFAB>& a_coarseCorr)
{
  CH_TIME("nwoebvto::amrProlong");
  CH_assert(!m_operatorForTimingOnly);
  Interval variables(0, SpaceDim-1);
  m_ebInterp.pwlInterp(a_correction, a_coarseCorr, variables);
}

/*****/
void
NWOEBViscousTensorOp::
AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                  const LevelData<EBCellFAB>& a_correction,
                  const LevelData<EBCellFAB>& a_coarseCorrection)
{
  CH_assert(!m_operatorForTimingOnly);
  LevelData<EBCellFAB> r;
  this->create(r, a_residual);
  this->assign(r, a_residual);
  this->AMRResidualNF(r, a_correction, a_coarseCorrection, a_residual, true);
  this->assign(a_residual, r);
}
/***/
void
NWOEBViscousTensorOp::
applyOp(LevelData<EBCellFAB>             & a_lhs,
        const LevelData<EBCellFAB>       & a_phi,
        bool                               a_homogeneous)
{
  CH_assert(!m_operatorForTimingOnly);
  CH_TIME("nwoebvto::applyOpFull");
  {
    //this includes exchange.
    CH_TIME("ghostcell fills");
    fillVelGhost(a_phi, a_homogeneous);
  }

  /**
     This is not done anymore because alpha and a are now part of the stencil
     and part of regular apply
     EBLevelDataOps::setToZero(a_lhs);
     incr( a_lhs, a_phi, m_alpha);
  **/
  {
    CH_TIME("applying op without bcs");
    DataIterator dit = m_eblg.getDBL().dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
      {

        applyOpRegular(  a_lhs[dit[mybox]], a_phi[dit[mybox]], a_homogeneous, dit[mybox]);
        applyOpIrregular(a_lhs[dit[mybox]], a_phi[dit[mybox]], a_homogeneous, dit[mybox]);
      }
  }
}

/*****/
void
NWOEBViscousTensorOp::
getFlux(FArrayBox&                    a_flux,
        const FArrayBox&              a_phi,
        const Box&                    a_faceBox,
        const int&                    a_idir,
        const DataIndex&              a_datInd)
{

  CH_TIME("nwoebvto::getFlux3");
  CH_assert(!m_operatorForTimingOnly);
  CH_assert(a_flux.nComp() == SpaceDim);
  CH_assert( a_phi.nComp() == SpaceDim);

  const FArrayBox& lamFace =(const FArrayBox&)((*m_lambda)[a_datInd][a_idir].getSingleValuedFAB());
  const FArrayBox& etaFace =(const FArrayBox&)((*m_eta   )[a_datInd][a_idir].getSingleValuedFAB());

  NWOViscousTensorOp::getFlux(a_flux, a_phi, etaFace, lamFace, a_faceBox, m_eblg.getDomain(), m_dx, m_beta, a_idir, 1);
}
/*****/
void
NWOEBViscousTensorOp::
applyOpRegular(EBCellFAB&             a_lhs,
               const EBCellFAB&       a_phi,
               const bool&            a_homogeneous,
               const DataIndex&       a_datInd)
{
  CH_TIME("nwoebvto::applyopRegularNoBCs");
  CH_assert(!m_operatorForTimingOnly);
  //assumes ghost cells already filled
  const Box& grid = m_eblg.getDBL()[a_datInd];
  BaseFab<Real>      & lphfab  =  a_lhs.getSingleValuedFAB();
  const BaseFab<Real>& phifab  =  a_phi.getSingleValuedFAB();
  const BaseFab<Real>& acofab  =(*m_acoef)[a_datInd].getSingleValuedFAB();
  const EBFluxFAB& eta  =     (*m_eta)[a_datInd];
  const EBFluxFAB& lam  =  (*m_lambda)[a_datInd];

  Vector<const BaseFab<Real>* > etaside(3, &(eta[0].getSingleValuedFAB()));
  Vector<const BaseFab<Real>* > lamside(3, &(lam[0].getSingleValuedFAB()));
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      etaside[idir] = &(eta[idir].getSingleValuedFAB());
      lamside[idir] = &(lam[idir].getSingleValuedFAB());
    }

  FORT_APPLYOPVTOPNOBCS(CHF_FRA(lphfab),
                        CHF_CONST_FRA(phifab),
                        CHF_CONST_FRA1(acofab,0),
                        CHF_CONST_FRA1((*etaside[0]), 0),
                        CHF_CONST_FRA1((*etaside[1]), 0),
                        CHF_CONST_FRA1((*etaside[2]), 0),
                        CHF_CONST_FRA1((*lamside[0]), 0),
                        CHF_CONST_FRA1((*lamside[1]), 0),
                        CHF_CONST_FRA1((*lamside[2]), 0),
                        CHF_CONST_REAL(m_dx),
                        CHF_CONST_REAL(m_alpha),
                        CHF_CONST_REAL(m_beta),
                        CHF_BOX(grid));
}
/*****/
void
NWOEBViscousTensorOp::
applyOpIrregular(EBCellFAB&             a_lhs,
                 const EBCellFAB&       a_phi,
                 const bool&            a_homogeneous,
                 const DataIndex&       a_datInd)
{
  CH_TIME("nwoebvto::applyOpIrregular");
  CH_assert(!m_operatorForTimingOnly);
  RealVect vectDx = m_dx*RealVect::Unit;
  for (int ivar = 0; ivar < SpaceDim; ivar++)
    {
      m_opEBStencil[ivar][a_datInd]->apply(a_lhs, a_phi, m_alphaDiagWeight[a_datInd], m_alpha, m_beta, ivar, false);
    }
  if (!a_homogeneous)
    {
      const Real factor = m_beta/m_dx;
      m_ebBC->applyEBFlux(a_lhs, a_phi, m_vofIterIrreg[a_datInd], (*m_eblg.getCFIVS()),
                          a_datInd, RealVect::Zero, vectDx, factor,
                          a_homogeneous, 0.0);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (int comp = 0; comp < SpaceDim; comp++)
        {
          for (m_vofIterDomLo[idir][a_datInd].reset(); m_vofIterDomLo[idir][a_datInd].ok();  ++m_vofIterDomLo[idir][a_datInd])
            {
              Real flux;
              const VolIndex& vof = m_vofIterDomLo[idir][a_datInd]();
              m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                      RealVect::Zero,vectDx,idir,Side::Lo, a_datInd, 0.0,
                                      a_homogeneous);

              //area gets multiplied in by bc operator
              a_lhs(vof,comp) -= flux*m_beta/m_dx;
            }
          for (m_vofIterDomHi[idir][a_datInd].reset(); m_vofIterDomHi[idir][a_datInd].ok();  ++m_vofIterDomHi[idir][a_datInd])
            {
              Real flux;
              const VolIndex& vof = m_vofIterDomHi[idir][a_datInd]();
              m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                      RealVect::Zero,vectDx,idir,Side::Hi,a_datInd,0.0,
                                      a_homogeneous);

              //area gets multiplied in by bc operator
              a_lhs(vof,comp) += flux*m_beta/m_dx;
            }
        }
    }
}
/**/
void
NWOEBViscousTensorOp::
cfinterp(const LevelData<EBCellFAB>&       a_phi,
         const LevelData<EBCellFAB>&       a_phiCoarse)
{
  CH_TIME("nwoebvto::cfinterp");
  CH_assert(!m_operatorForTimingOnly);
  if (m_hasCoar)
    {
      LevelData<EBCellFAB>& phi = (LevelData<EBCellFAB>&)a_phi;
      if (a_phiCoarse.isDefined())
        {
          m_interpWithCoarser->coarseFineInterp(phi, a_phiCoarse, 0, 0, SpaceDim);
        }
    }
}
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
