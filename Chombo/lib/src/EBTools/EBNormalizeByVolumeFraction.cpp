#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBNormalizeByVolumeFraction.H"
#include "EBArith.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "CH_OpenMP.H"
#include "NamespaceHeader.H"
using namespace std;
//----------------------------------------------------------------------------
void 
EBNormalizeByVolumeFraction::
getLocalStencil(VoFStencil      & a_stencil, 
                const VolIndex  & a_vof, 
                const DataIndex & a_dit)
{
  Vector<VolIndex> neighbors;
  EBArith::getAllVoFsWithinRadius(neighbors, a_vof, m_eblg.getEBISL()[a_dit], m_radius);
  Real sumkap = 0;
  a_stencil.clear();
  for(int ivof = 0; ivof < neighbors.size(); ivof++)
    {
      Real kappa = m_eblg.getEBISL()[a_dit].volFrac(neighbors[ivof]);
      //data is assumed to already be kappa weighted
      a_stencil.add(neighbors[ivof], 1.0);
      sumkap += kappa;
    }
  if(sumkap > 1.0e-9)
    {
      a_stencil *= (1.0/sumkap);
    }
}
//----------------------------------------------------------------------------
void 
EBNormalizeByVolumeFraction::
define(const LevelData<EBCellFAB> & a_Q)
{
  CH_TIME("EBNormalizeByVolFrac::define");
  m_isDefined = true;
  m_stencil.define(m_eblg.getDBL());

  DataIterator dit = m_eblg.getDBL().dataIterator(); 
  int nbox=dit.size();
#pragma omp parallel for
  for (int mybox=0;mybox<nbox; mybox++)
    {
      CH_TIME("vof stencil definition");
      const     Box& grid = m_eblg.getDBL()  [dit[mybox]];
      const EBISBox& ebis = m_eblg.getEBISL()[dit[mybox]];
 
      IntVectSet ivs = ebis.getIrregIVS(grid);
      VoFIterator vofit(ivs, ebis.getEBGraph());
      const Vector<VolIndex>& vofvec = vofit.getVector();
      // cast from VolIndex to BaseIndex
      Vector<RefCountedPtr<BaseIndex> >    dstVoF(vofvec.size());
      Vector<RefCountedPtr<BaseStencil> > stencil(vofvec.size());
      // fill stencils for the vofs
      for(int ivec = 0; ivec < vofvec.size(); ivec++)
        {
          VoFStencil localStencil;
          getLocalStencil(localStencil, vofvec[ivec], dit[mybox]);

          // another cast from VolIndex to BaseIndex
          dstVoF[ivec]  = RefCountedPtr<BaseIndex  >(new  VolIndex(vofvec[ivec]));
          stencil[ivec] = RefCountedPtr<BaseStencil>(new VoFStencil(localStencil));
        }
      m_stencil[dit[mybox]] = RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB > >
        (new AggStencil<EBCellFAB, EBCellFAB >(dstVoF, stencil, a_Q[dit[mybox]], a_Q[dit[mybox]]));

    }//dit
}
//----------------------------------------------------------------------------
void
EBNormalizeByVolumeFraction::
operator()(LevelData<EBCellFAB>& a_Q,
           const Interval& a_compInterval)
{
  CH_TIME("EBNormalizer::operator()");
  if(!m_isDefined)
    {
      define(a_Q);
    }

  a_Q.exchange();
  EBCellFactory fact(m_eblg.getEBISL());

  //for this to be correct, need the source and destination to be different
//  LevelData<EBCellFAB> save(m_eblg.getDBL(), a_Q.nComp(), a_Q.ghostVect(), fact);
//  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
//    {
//      save[dit()].copy(a_Q[dit()]);
//    }
  // Endpoints of the given interval.
  int begin  = a_compInterval.begin();
  int length = a_compInterval.size();

  //save.exchange();
   
  DataIterator dit = m_eblg.getDBL().dataIterator(); 
  int nbox=dit.size();
#pragma omp parallel for
  for (int mybox=0;mybox<nbox; mybox++)
    {
      EBCellFAB save(fact.ebisBox(dit[mybox]), grow(m_eblg.getDBL()[dit[mybox]],a_Q.ghostVect()),a_Q.nComp()) ;
      save.copy(a_Q[dit[mybox]]);
      m_stencil[dit[mybox]]->apply(a_Q[dit[mybox]], save, begin, begin, length, false);
    }
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
void
EBNormalizeByVolumeFraction::
operator()(LevelData<EBCellFAB>& a_Q)
{
  return (*this)(a_Q, Interval(0, a_Q.nComp()-1));
}
//----------------------------------------------------------------------------


#include "NamespaceFooter.H"
