#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "nwoebquadcfinterp.H"

#include "LayoutIterator.H"
#include "DataIterator.H"
#include "EBAlias.H"
#include "EBArith.H"
#include "EBCellFactory.H"
#include "EBLevelDataOps.H"
#include "NamespaceHeader.H"

/***********************/
nwoebquadcfinterp::
nwoebquadcfinterp(const DisjointBoxLayout&       a_gridsFine,
                  const DisjointBoxLayout&       a_gridsCoar,
                  const EBISLayout&              a_ebislFine,
                  const EBISLayout&              a_ebislCoar,
                  const ProblemDomain&           a_domainCoar,
                  const int&                     a_nref,
                  const int&                     a_nvar,
                  const Real&                    a_dxFine,
                  const IntVect &                a_ghost,
                  const LayoutData<IntVectSet>&  a_cfivs,
		  const EBIndexSpace* const a_ebisPtr)             
{
  m_gridsFine  =  a_gridsFine ;          
  m_gridsCoar  =  a_gridsCoar ;     
  m_ebislFine  =  a_ebislFine ;      
  m_ebislCoar  =  a_ebislCoar ;      
  m_domainCoar =  a_domainCoar;
  m_nref       =  a_nref      ;      
  m_nvar       =  a_nvar      ;      
  m_dxFine     =  a_dxFine    ;      
  m_ghost      =  a_ghost;
  defineInternals(a_cfivs, a_ebisPtr);
}
void 
nwoebquadcfinterp::
defineInternals(const LayoutData<IntVectSet>&  a_cfivs, const EBIndexSpace* const a_ebisPtr)             
{
  m_gridsCoFi= DisjointBoxLayout();
  coarsen(m_gridsCoFi, m_gridsFine, m_nref);

  a_ebisPtr->fillEBISLayout(m_ebislCoFi,
                          m_gridsCoFi,
                          m_domainCoar, 4);

  EBCellFactory fact(m_ebislCoFi);
  m_bufferCoFi.define(m_gridsCoFi, m_nvar, 4*IntVect::Unit, fact);

  defineStencils(a_cfivs);
}
void 
nwoebquadcfinterp::
defineStencils(const LayoutData<IntVectSet>&  a_cfivs)             
{
  EBCellFactory fact(m_ebislFine);
  LevelData<EBCellFAB> proxyLevel(m_gridsFine, 1, m_ghost, fact);
  m_stencil.define(m_gridsFine);
  for(DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisCoFi = m_ebislCoFi[dit()];
      const EBISBox& ebisFine = m_ebislFine[dit()];
      Vector< RefCountedPtr<BaseIndex  > > baseDstVoFs;
      Vector< RefCountedPtr<BaseStencil> > baseSten;
      getStencils(baseSten, baseDstVoFs, a_cfivs[dit()], ebisFine, ebisCoFi,  dit());

      EBCellFAB& coarProxy = m_bufferCoFi[dit()];
      EBCellFAB& fineProxy = proxyLevel[dit()];
      m_stencil[dit()] = RefCountedPtr<AggStencil <EBCellFAB, EBCellFAB>  >
        (new AggStencil<EBCellFAB, EBCellFAB >(baseDstVoFs, baseSten, coarProxy, fineProxy));
    }
}
/***********************/
void
nwoebquadcfinterp::
getStencils(Vector<RefCountedPtr< BaseStencil> >  & a_stencils, 
            Vector<RefCountedPtr< BaseIndex  > >  & a_baseDstVoFs,
            const IntVectSet                      & a_cfivs,
            const EBISBox                         & a_ebisFine,
            const EBISBox                         & a_ebisCoFi,
            const DataIndex                       & a_dit)
{
  VoFIterator vofit(a_cfivs, a_ebisFine.getEBGraph());
  const Vector<VolIndex>& volvec = vofit.getVector();
  a_baseDstVoFs.resize(volvec.size());
  a_stencils.resize(   volvec.size());
  for(int ivec = 0; ivec < volvec.size(); ivec++)
    {
      VoFStencil sten;
      getStencil(sten,  volvec[ivec], a_ebisFine, a_ebisCoFi, a_dit);
      a_stencils    [ivec]  = RefCountedPtr<BaseStencil>(new VoFStencil(sten));
      a_baseDstVoFs [ivec]  = RefCountedPtr<BaseIndex  >(new VolIndex(volvec[ivec]));
    }
}
/***********************/
void
nwoebquadcfinterp::
getStencil(VoFStencil           & a_stencil,
           const VolIndex       & a_vofFine,
           const EBISBox        & a_ebisFine,
           const EBISBox        & a_ebisCoFi,
           const DataIndex      & a_dit)
{
  VolIndex fineVoF = a_vofFine;
  Real dxFine = m_dxFine;  Real dxCoar = m_nref*m_dxFine;
  a_stencil.clear();
  VolIndex coarVoF = a_ebisFine.coarsen(a_vofFine);
  RealVect coarLoc = EBArith::getVofLocation(coarVoF, dxCoar*RealVect::Unit, RealVect::Zero);
  RealVect fineLoc = EBArith::getVofLocation(fineVoF, dxFine*RealVect::Unit, RealVect::Zero);
  RealVect dist = fineLoc - coarLoc;
  EBArith::getExtrapolationStencil(a_stencil, dist, dxCoar*RealVect::Unit, coarVoF, a_ebisCoFi);
}  
/***********************/
void 
nwoebquadcfinterp::
coarseFineInterp(LevelData<EBCellFAB>&       a_phif,
                 const LevelData<EBCellFAB>& a_phic,
                 const int                   a_isrc,
		 const int                   a_idst,
		 const int                   a_inco)
{
  const Interval interv(0, a_inco-1);
  a_phic.copyTo(interv, m_bufferCoFi, interv);
  DataIterator dit = m_gridsFine.dataIterator();
  int nbox = dit.size();
#pragma omp parallel for
  for(int ibox = 0; ibox < nbox; ibox++)
    {
      //false is for increment only
      m_stencil[dit[ibox]]->apply(a_phif[dit[ibox]],
                                  m_bufferCoFi[dit[ibox]],
                                  a_isrc, a_idst, a_inco, false);
    }
}

/***********************/
void 
nwoebquadcfinterp::
coarseFineInterpH(LevelData<EBCellFAB>& a_phif,
                  const int             a_isrc,
		  const int             a_idst,
		  const int             a_inco)
{
  EBLevelDataOps::setVal(m_bufferCoFi, 0.0);
  DataIterator dit = m_gridsFine.dataIterator();
  const int nbox = dit.size();
#pragma omp parallel for
  for(int ibox = 0; ibox < nbox; ibox++)
    {
      //false is for increment only
      m_stencil[dit[ibox]]->apply(a_phif[dit[ibox]],
                                  m_bufferCoFi[dit[ibox]],
                                  a_isrc, a_idst, a_inco, false);
    }

}

#include "NamespaceFooter.H"
