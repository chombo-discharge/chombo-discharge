#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif
#include <cmath>
#include "MFAMRDataOps.H"
#include "MFArith.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

/**********/
void MFAMRDataOps::setCoveredAMRVal(Vector<LevelData<MFCellFAB>* >&      a_data,
                                    const Vector< Vector<EBISLayout> >&  a_ebislv,
                                    const Vector<int>&                   a_refRat,
                                    const Real&                          a_value)
{
  CH_TIME("MFAMRDataOps::setCoveredAMRVal");
  int numLevels = a_data.size();

  for (int ilev = 0; ilev < numLevels-1; ilev++)
    {
      int nComp = a_data[ilev]->nComp();
      int coarLev = ilev;
      int fineLev = ilev+1;
      Interval interv(0, nComp-1);
      const DisjointBoxLayout& coarGrids = a_data[coarLev]->disjointBoxLayout();
      const DisjointBoxLayout& fineGrids = a_data[fineLev]->disjointBoxLayout();

      //zero out under finer grids
      for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
        {
          MFCellFAB& coarFAB = (*a_data[coarLev])[dit()];

          LayoutIterator litFine = fineGrids.layoutIterator();
          for (litFine.reset(); litFine.ok(); ++litFine)
            {
              Box overlayBox = coarFAB.box();
              Box coarsenedGrid = coarsen(fineGrids[litFine()], a_refRat[coarLev]);

              overlayBox &= coarsenedGrid;
              if (!overlayBox.isEmpty())
                {
                  IntVectSet ivsZero(overlayBox);
                  for (int iphase=0; iphase<a_ebislv.size(); iphase++)
                    {
                      EBCellFAB& coarEBFAB = coarFAB.getPhase(iphase);
                      const EBISBox& ebisBox = a_ebislv[iphase][coarLev][dit()];
                      for (VoFIterator vofit(ivsZero, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                        {
                          for (int ivar =0; ivar < nComp; ivar++)
                            {
                              coarEBFAB(vofit(), ivar) = a_value;
                            }
                        }
                    }
                }
            }
        }
    }
}
/**********/
void MFAMRDataOps::setCoveredAMRVal(Vector<LevelData<MFCellFAB>* >&       a_data,
                                    const Vector< Vector<EBLevelGrid> >&  a_eblgv,
                                    const Vector<int>&                    a_refRat,
                                    const Real&                           a_value)
{
  CH_TIME("MFAMRDataOps::setCoveredAMRVal(eblg)");
  int numLevels = a_data.size();

  for (int ilev = 0; ilev < numLevels-1; ilev++)
    {
      int nComp = a_data[ilev]->nComp();
      int coarLev = ilev;
      int fineLev = ilev+1;
      Interval interv(0, nComp-1);
      const DisjointBoxLayout& coarGrids = a_data[coarLev]->disjointBoxLayout();
      const DisjointBoxLayout& fineGrids = a_data[fineLev]->disjointBoxLayout();

      //zero out under finer grids
      for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
        {
          MFCellFAB& coarFAB = (*a_data[coarLev])[dit()];

          LayoutIterator litFine = fineGrids.layoutIterator();
          for (litFine.reset(); litFine.ok(); ++litFine)
            {
              Box overlayBox = coarFAB.box();
              Box coarsenedGrid = coarsen(fineGrids[litFine()], a_refRat[coarLev]);

              overlayBox &= coarsenedGrid;
              if (!overlayBox.isEmpty())
                {
                  IntVectSet ivsZero(overlayBox);
                  for (int iphase=0; iphase<a_eblgv.size(); iphase++)
                    {
                      EBCellFAB& coarEBFAB = coarFAB.getPhase(iphase);
                      const EBISBox& ebisBox = a_eblgv[iphase][coarLev].getEBISL()[dit()];
                      for (VoFIterator vofit(ivsZero, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                        {
                          for (int ivar =0; ivar < nComp; ivar++)
                            {
                              coarEBFAB(vofit(), ivar) = a_value;
                            }
                        }
                    }
                }
            }
        }
    }
}
/**********/
void MFAMRDataOps::setCoveredVal(Vector<LevelData<MFCellFAB>* >& a_data,
                                 const Real&                     a_value)
{
  CH_TIME("MFAMRDataOps::setCoveredVal");
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::setCoveredVal(*a_data[ilev],
                                    a_value);
    }

}
/**********/
void MFAMRDataOps::setToZero(Vector<LevelData<MFCellFAB>* >& a_result)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::setToZero(*a_result[ilev]);
    }
}
/**********/
void MFAMRDataOps::setVal(Vector<LevelData<MFCellFAB>* >& a_result,
                          const Real&                     a_value)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::setVal(*a_result[ilev],
                             a_value);
    }
}
/**********/
void MFAMRDataOps::setVal(Vector<LevelData<MFCellFAB>* >& a_result,
                          const Real&                     a_value,
                          const int&                      a_comp)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::setVal(*a_result[ilev],
                             a_value,
                             a_comp);
    }
}
/**********/
void MFAMRDataOps::axby(Vector<LevelData<MFCellFAB>* >&          a_lhs,
                        const Vector<LevelData<MFCellFAB>* >&    a_x,
                        const Vector<LevelData<MFCellFAB>* >&    a_y,
                        const Real&                              a_a,
                        const Real&                              a_b)
{
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::axby(*a_lhs[ilev],
                           *a_x[ilev],
                           *a_y[ilev],
                           a_a,
                           a_b);
    }
}
/**********/
void MFAMRDataOps::axby(Vector<LevelData<MFCellFAB>* >&       a_lhs,
                        const Vector<LevelData<MFCellFAB>* >& a_x,
                        const Vector<LevelData<MFCellFAB>* >& a_y,
                        const Real&                           a_a,
                        const Real&                           a_b,
                        const int&                            a_lhsComp,
                        const int&                            a_xComp,
                        const int&                            a_yComp)
{
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::axby(*a_lhs[ilev],
                           *a_x[ilev],
                           *a_y[ilev],
                           a_a,
                           a_b,
                           a_lhsComp,
                           a_xComp,
                           a_yComp);
    }
}
/**********/
void MFAMRDataOps::assign(Vector<LevelData<MFCellFAB>* >&       a_to,
                          const Vector<LevelData<MFCellFAB>* >& a_from,
                          const Interval&                       a_toInterval,
                          const Interval&                       a_fromInterval)
{
  int numLevels = a_to.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::assign(*a_to[ilev],
                             *a_from[ilev],
                             a_toInterval,
                             a_fromInterval);
    }
}
/**********/
void MFAMRDataOps::assign(Vector<RefCountedPtr<LevelData<MFCellFAB> > >&  a_to,
                          const Vector<LevelData<MFCellFAB>* >&           a_from,
                          const Interval&                                 a_toInterval,
                          const Interval&                                 a_fromInterval)
{
  CH_TIME("MFAMRDataOps::assign(cell_to,cell_from,int_to,int_from)");
  int numLevels = a_to.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::assign(*a_to[ilev],
                             *a_from[ilev],
                             a_toInterval,
                             a_fromInterval);
    }
}
/**********/
void MFAMRDataOps::assign(Vector<LevelData<MFCellFAB>* >&       a_lhs,
                          const Vector<LevelData<MFCellFAB>* >& a_rhs)
{
  CH_TIME("MFAMRDataOps::assign(cell_lhs,cell_rhs)");
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::assign(*a_lhs[ilev],
                             *a_rhs[ilev]);
    }
}
/**********/
void MFAMRDataOps::assign(Vector<RefCountedPtr<LevelData<MFCellFAB> > >&  a_lhs,
                          const Vector<LevelData<MFCellFAB>* >&           a_rhs)
{
  CH_TIME("MFAMRDataOps::assign(cellRef_lhs,cellRef_rhs)");
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::assign(*a_lhs[ilev],
                             *a_rhs[ilev]);
    }
}
/**********/
void MFAMRDataOps::incr(Vector<LevelData<MFCellFAB>* >&       a_lhs,
                        const Vector<LevelData<MFCellFAB>* >& a_rhs,
                        const Real&                           a_scale)
{
  CH_TIME("MFAMRDataOps::incr(lhs,rhs,scale)");
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::incr(*a_lhs[ilev],
                           *a_rhs[ilev],
                           a_scale);
    }
}
/**********/
void MFAMRDataOps::incr(Vector<LevelData<MFCellFAB>* >&  a_lhs,
                        const Real&                      a_scale)
{
  CH_TIME("MFAMRDataOps::incr(lhs,scale)");
  int numLevels = a_lhs.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::incr(*a_lhs[ilev],
                           a_scale);
    }
}
/**********/
void MFAMRDataOps::scale(Vector<LevelData<MFCellFAB>* >& a_result,
                         const Real&                     a_value)
{
  CH_TIME("MFAMRDataOps::scale(cell_result,scale)");
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::scale(*a_result[ilev],
                            a_value);
    }
}
/**********/
void MFAMRDataOps::scale(Vector<LevelData<MFCellFAB>* >& a_result,
                         const Real&                     a_value,
                         const int&                      a_comp)
{
  CH_TIME("MFAMRDataOps::scale(flux_result,scale,comp)");
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::scale(*a_result[ilev],
                            a_value,
                            a_comp);
    }
}
/**********/
void MFAMRDataOps::sum(Vector<LevelData<MFCellFAB>* >&       a_result,
                       const Vector<LevelData<MFCellFAB>* >& a_in1,
                       const Vector<LevelData<MFCellFAB>* >& a_in2)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::sum(*a_result[ilev],
                          *a_in1[ilev],
                          *a_in2[ilev]);
    }
}
/**********/
void MFAMRDataOps::addConstant(Vector<LevelData<MFCellFAB>* >& a_data,
                               const Real&                     a_constant)
{
  for (int ilev = 0; ilev < a_data.size(); ilev++)
    {
      MFLevelDataOps::addConstant(*a_data[ilev], a_constant);
    }
}
/**********/
void MFAMRDataOps::product(Vector<LevelData<MFCellFAB>* >&       a_result,
                           const Vector<LevelData<MFCellFAB>* >& a_in1,
                           const Vector<LevelData<MFCellFAB>* >& a_in2)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::product(*a_result[ilev],
                              *a_in1[ilev],
                              *a_in2[ilev]);
    }
}
/**********/
void MFAMRDataOps::product(Vector<LevelData<MFCellFAB>* >&       a_result,
                           const Vector<LevelData<MFCellFAB>* >& a_in1,
                           const Vector<LevelData<MFCellFAB>* >& a_in2,
                           const int&                            a_rComp,
                           const int&                            a_1Comp,
                           const int&                            a_2Comp)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::product(*a_result[ilev],
                              *a_in1[ilev],
                              *a_in2[ilev],
                              a_rComp,
                              a_1Comp,
                              a_2Comp);
    }
}
/**********/
void MFAMRDataOps::divideVectorByScalar(Vector<LevelData<MFCellFAB>* >&       a_vectorOut,
                                        const Vector<LevelData<MFCellFAB>* >& a_vectorIn,
                                        const Vector<LevelData<MFCellFAB>* >& a_scalar)
{
  //a_vector /= a_scalar
  int numLevels = a_vectorOut.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::divideVectorByScalar(*a_vectorOut[ilev],
                                           *a_vectorIn[ ilev],
                                           *a_scalar[   ilev]);
    }
}
/**********/
void MFAMRDataOps::divide(Vector<LevelData<MFCellFAB>* >&       a_result,
                          const Vector<LevelData<MFCellFAB>* >& a_in1,
                          const Vector<LevelData<MFCellFAB>* >& a_in2)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::divide(*a_result[ilev],
                              *a_in1[ilev],
                              *a_in2[ilev]);
    }
}
/**********/
void MFAMRDataOps::divide(Vector<LevelData<MFCellFAB>* >&       a_result,
                          const Vector<LevelData<MFCellFAB>* >& a_in1,
                          const Vector<LevelData<MFCellFAB>* >& a_in2,
                          const int&                            a_rComp,
                          const int&                            a_1Comp,
                          const int&                            a_2Comp)
{
  int numLevels = a_result.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::divide(*a_result[ilev],
                             *a_in1[ilev],
                             *a_in2[ilev],
                             a_rComp,
                             a_1Comp,
                             a_2Comp);
    }
}
/**********/
void MFAMRDataOps::kappaWeight(Vector<LevelData<MFCellFAB>* >& a_data)
{
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::kappaWeight(*a_data[ilev]);
    }
}
/**********/
void MFAMRDataOps::kappaScale(Vector<LevelData<MFCellFAB>* >& a_data,
                              const Real&                     a_value)
{
  int numLevels = a_data.size();
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      MFLevelDataOps::kappaScale(*a_data[ilev],
                                 a_value);
    }
}
/**********/
Real MFAMRDataOps::subtractOffMean(Vector<LevelData<MFCellFAB>* >&     a_data,
                                   const Vector<DisjointBoxLayout>&    a_grids,
                                   const Vector< Vector<EBISLayout> >& a_ebislv,
                                   const Vector<int>&                  a_refRat)
{
  int pmode = -2;  //norm = (1/v)(sum(a_data dv)) ---no absolute values and multiply kappa as you go
  const Real mean = MFArith::norm(a_data, a_grids, a_ebislv, a_refRat, 0, pmode);

  //subtract off mean of data
  incr(a_data, -mean);
  return mean;
}
/**********/
Real MFAMRDataOps::subtractOffMean(Vector< LevelData<MFCellFAB>* >&      a_data,
                                   const Vector< Vector<EBLevelGrid> >&  a_eblgv,
                                   const Vector<int>&                    a_refRat)
{
  int numLevels = a_data.size();
  int numPhases = a_eblgv.size();
  Vector< Vector<EBISLayout> > ebislv(numPhases);
  ebislv[0].resize(numLevels);
  ebislv[1].resize(numLevels);
  Vector<DisjointBoxLayout>    grids(numPhases);
  Vector<ProblemDomain>        domain(numPhases);
  for (int ilev=0; ilev<numLevels; ilev++)
    {
      ebislv[0][ilev]  = a_eblgv[0][ilev].getEBISL();
      ebislv[1][ilev]  = a_eblgv[1][ilev].getEBISL();
      grids[ilev]  = a_eblgv[0][ilev].getDBL();
      domain[ilev] = a_eblgv[0][ilev].getDomain();
    }

  int pmode = -2;  //norm = (1/v)(sum(a_data dv)) ---no absolute values and multiply kappa as you go
  const Real mean = MFArith::norm(a_data, grids, ebislv, a_refRat, 0, pmode);

  //subtract off mean of data
  incr(a_data, -mean);
  return mean;
}

#include "NamespaceFooter.H"
