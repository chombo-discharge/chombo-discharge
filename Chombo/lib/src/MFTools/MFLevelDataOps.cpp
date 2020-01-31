#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "SPMD.H"

#include "MFLevelDataOps.H"
#include "EBLevelDataOpsF_F.H"
#include "FaceIterator.H"
#include "VoFIterator.H"
#include "BoxIterator.H"
#include "CH_Timer.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "EBIndexSpace.H"
#include "EBLoadBalance.H"
#include "CornerCopier.H"
#include "NamespaceHeader.H"

/*****/
void MFLevelDataOps::setCoveredVal(LevelData<MFCellFAB>& a_data,
                                   const Real&           a_value)
{
  CH_TIME("MFLevelDataOps::setCoveredVal(cell)");
  int ncomp = a_data.nComp();
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      MFCellFAB& dataMF = a_data[dit()];
      for (int iphase=0; iphase<dataMF.numPhases(); iphase++)
        {
          EBCellFAB& dataEB = dataMF.getPhase(iphase);
          for (int icomp = 0; icomp < ncomp;++icomp)
            {
              dataEB.setCoveredCellVal(a_value,icomp);
            }
        }
    }
}
/*****/
void MFLevelDataOps::setIrregVal(LevelData<MFCellFAB>&     a_data,
                                 const DisjointBoxLayout&  a_dbl,
                                 const Vector<EBISLayout>& a_ebislv,
                                 const Real&               a_value)
{
  int ncomp = a_data.nComp();
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      for (int icomp=0; icomp<ncomp;++icomp)
        {
          const Box& grid = a_dbl.get(dit());
          MFCellFAB& dataMF = a_data[dit()];
          for (int iphase=0; iphase<dataMF.numPhases(); iphase++)
            {
              EBCellFAB& dataEB = dataMF.getPhase(iphase);
              const EBISBox& ebisBox = a_ebislv[iphase][dit()];
              const IntVectSet& ivsIrreg = ebisBox.getIrregIVS(grid);
              for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                {
                  const VolIndex vof = vofit();
                  dataEB(vof,icomp) = a_value;
                }
            }
        }
    }
}
/*****/
void MFLevelDataOps::setCoveredVal(LevelData<MFCellFAB>& a_data,
                                   const int&            a_comp,
                                   const Real&           a_value)
{
  CH_TIME("MFLevelDataOps::setCoveredVal(cell,comp)");
  for (DataIterator dit = a_data.dataIterator();dit.ok();++dit)
    {
      MFCellFAB& dataMF = a_data[dit()];
      for (int iphase=0; iphase<dataMF.numPhases(); iphase++)
        {
          EBCellFAB& dataEB = dataMF.getPhase(iphase);
          dataEB.setCoveredCellVal(a_value, a_comp);
        }
    }
}
/*****/
void MFLevelDataOps::averageMultiVofsToRegFAB(LevelData<MFCellFAB>&     a_data,
                                              const DisjointBoxLayout&  a_dbl,
                                              const Vector<EBISLayout>& a_ebislv)
{
  const IntVect& ghostVect = a_data.ghostVect();
  const int&         nComp = a_data.nComp();
  const int&     numPhases = a_ebislv.size();
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box& dblBox = a_dbl.get(dit());
      MFCellFAB& dataMF = a_data[dit()];
      for (int iphase=0; iphase<numPhases; iphase++)
        {
          EBCellFAB&         data     = dataMF.getPhase(iphase);
          BaseFab<Real>&     regData  = data.getSingleValuedFAB();
          const EBISBox&     ebisBox  = a_ebislv[iphase][dit()];

          Box bigBox = dblBox;
          bigBox.grow(ghostVect);
          const IntVectSet&  multiIVS = ebisBox.getMultiCells(bigBox);
          for (IVSIterator ivsit(multiIVS); ivsit.ok(); ++ivsit)
            {
              const IntVect&            iv = ivsit();
              const Vector<VolIndex>& vofs = ebisBox.getVoFs(iv);

              for (int comp = 0; comp < nComp; comp++)
                {
                  Real sumData = 0.0;
                  Real sumFrac = 0.0;
                  for (int i = 0; i < vofs.size(); i++)
                    {
                      const VolIndex& thisVof  = vofs[i];
                      const Real&     kappa    = ebisBox.volFrac(thisVof);
                      sumFrac += kappa;
                      sumData += kappa*data(thisVof,comp);
                    }
                  sumData /= sumFrac;
                  regData(iv, comp) = sumData;
                }
            }
        }
    }
}
/*****/
void MFLevelDataOps::copyToMultiVofsFromRegFAB(LevelData<MFCellFAB>&     a_data,
                                               const DisjointBoxLayout&  a_dbl,
                                               const Vector<EBISLayout>& a_ebislv)
{
  const IntVect& ghostVect = a_data.ghostVect();
  const int&         nComp = a_data.nComp();
  const int&     numPhases = a_ebislv.size();
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box& dblBox = a_dbl.get(dit());
      MFCellFAB& dataMF = a_data[dit()];
      for (int iphase=0; iphase<numPhases; iphase++)
        {
          EBCellFAB&         data     = dataMF.getPhase(iphase);
          BaseFab<Real>&     regData  = data.getSingleValuedFAB();
          const EBISBox&     ebisBox  = a_ebislv[iphase][dit()];

          Box bigBox = dblBox;
          bigBox.grow(ghostVect);
          const IntVectSet&  multiIVS = ebisBox.getMultiCells(bigBox);
          for (IVSIterator ivsit(multiIVS); ivsit.ok(); ++ivsit)
            {
              const IntVect&            iv = ivsit();
              const Vector<VolIndex>& vofs = ebisBox.getVoFs(iv);

              for (int comp = 0; comp < nComp; comp++)
                {
                  const Real& rData = regData(iv,comp);
                  for (int i = 0; i < vofs.size(); i++)
                    {
                      const VolIndex& thisVof  = vofs[i];
                      data(thisVof,comp) = rData;
                    }
                }
            }
        }
    }
}
/*****/
void MFLevelDataOps::defineLevelData(LevelData<MFCellFAB>&     a_levelData,
                                     const Vector<EBISLayout>& a_ebislv,
                                     const DisjointBoxLayout&  a_dbl,
                                     const IntVect&            a_ghosts,
                                     const int&                a_nComp)
{
  Vector<int> nComp(a_ebislv.size(), a_nComp);
  Vector<EBISLayout> ebislv = a_ebislv;
  MFCellFactory mfcellfact(ebislv, nComp);
  a_levelData.define(a_dbl, a_nComp, a_ghosts, mfcellfact);
}
/*****/
void MFLevelDataOps::setToZero(LevelData<MFCellFAB>& a_result)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      a_result[dit()].setVal(0.0);
    }
}
/*****/
void MFLevelDataOps::setVal(LevelData<MFCellFAB>& a_result,
                            const Real&           a_value)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      a_result[dit()].setVal(a_value);
    }
}
/*****/
void MFLevelDataOps::setVal(LevelData<MFCellFAB>& a_result,
                            const Real&           a_value,
                            const int&            a_comp)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      MFCellFAB& result = a_result[dit()];
      for (int iphase=0; iphase<result.numPhases(); iphase++)
        {
          EBCellFAB& data = result.getPhase(iphase);
          data.setVal(a_comp,a_value);
        }
    }
}
/*****/
void MFLevelDataOps::axby( LevelData<MFCellFAB>&       a_lhs,
                           const LevelData<MFCellFAB>& a_x,
                           const LevelData<MFCellFAB>& a_y,
                           const Real&                 a,
                           const Real&                 b)
{
  //  CH_assert(a_lhs.disjointBoxLayout() == a_x.disjointBoxLayout());

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      MFCellFAB& lhs = a_lhs[dit()];
      const MFCellFAB& x = a_x[dit()];
      const MFCellFAB& y = a_y[dit()];
      for (int iphase=0; iphase<lhs.numPhases(); iphase++)
        {
          EBCellFAB& dataLHS = lhs.getPhase(iphase);
          const EBCellFAB& dataX = x.getPhase(iphase);
          const EBCellFAB& dataY = y.getPhase(iphase);
          dataLHS.axby(dataX, dataY, a, b);
        }
    }
}
/*****/
/*
void MFLevelDataOps::axby( LevelData<MFCellFAB>&       a_lhs,
                           const LevelData<MFCellFAB>& a_x,
                           const LevelData<MFCellFAB>& a_y,
                           const Real&                 a,
                           const Real&                 b,
                           const int&                  a_lhsComp,
                           const int&                  a_xComp,
                           const int&                  a_yComp)
{
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      MFCellFAB& lhs = a_lhs[dit()];
      const MFCellFAB& x = a_x[dit()];
      const MFCellFAB& y = a_y[dit()];
      for (int iphase=0; iphase<lhs.numPhases(); iphase++)
        {
          EBCellFAB& dataLHS = lhs.getPhase(iphase);
          const EBCellFAB& dataX = x.getPhase(iphase);
          const EBCellFAB& dataY = y.getPhase(iphase);
          dataLHS.axby(dataX, dataY, a, b, a_lhsComp, a_xComp, a_yComp);
        }
    }
}
*/
/*****/
void MFLevelDataOps::assign(LevelData<MFCellFAB>&       a_to,
                            const LevelData<MFCellFAB>& a_from,
                            const Interval&             a_toInterval,
                            const Interval&             a_fromInterval)
{
  CH_TIME("MFLevelDataOps::assign(to,from,toInterval,fromInterval)");
  a_from.copyTo(a_fromInterval, a_to, a_toInterval);
}
/*****/
void MFLevelDataOps::assign(LevelData<MFCellFAB>&       a_lhs,
                            const LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFLevelDataOps::assign(to,from)");
  Interval interv(0, a_rhs.nComp()-1);
  a_rhs.copyTo(interv, a_lhs, interv);
}
/*****/
void MFLevelDataOps::clone(LevelData<MFCellFAB>&       a_lhs,
                           const LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFLevelDataOps::clone");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      MFCellFAB& lhs = a_lhs[dit()];
      const MFCellFAB& rhs = a_rhs[dit()];
      for (int iphase=0; iphase<lhs.numPhases(); iphase++)
        {
          EBCellFAB& dataLHS = lhs.getPhase(iphase);
          const EBCellFAB& dataRHS = rhs.getPhase(iphase);
          dataLHS.copy(dataRHS);
        }
    }
}
/*****/
void MFLevelDataOps::incr(LevelData<MFCellFAB>&       a_lhs,
                          const LevelData<MFCellFAB>& a_rhs,
                          const Real&                 a_scale)
{
  //  CH_assert(a_lhs.disjointBoxLayout() == a_rhs.disjointBoxLayout());
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      a_lhs[d].plus(a_rhs[d], a_scale);
    }
}
/*****/
void MFLevelDataOps::incr( LevelData<MFCellFAB>& a_lhs,
                           const Real&           a_scale)
{
  //  CH_assert(a_lhs.disjointBoxLayout() == a_rhs.disjointBoxLayout());
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      a_lhs[dit()] += a_scale;
    }
}
/*****/
void MFLevelDataOps::scale(LevelData<MFCellFAB>& a_result,
                           const Real&           a_value)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      a_result[dit()].mult(a_value);
    }
}
/*****/
void MFLevelDataOps::scale(LevelData<MFCellFAB>& a_result,
                           const Real&           a_value,
                           const int&            a_comp)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      MFCellFAB& result = a_result[d];
      for (int iphase=0; iphase<result.numPhases(); iphase++)
        {
          EBCellFAB& dataEB = result.getPhase(iphase);
          dataEB.mult(a_value, a_comp, 1);
                   
        }
    }
}
/*****/
void MFLevelDataOps::sum(LevelData<MFCellFAB>&       a_result,
                         const LevelData<MFCellFAB>& a_in1,
                         const LevelData<MFCellFAB>& a_in2)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      MFCellFAB& result = a_result[d];
      const MFCellFAB& in1 = a_in1[d];
      const MFCellFAB& in2 = a_in2[d];

      result.copy(in1);
      result += in2;
    }
}
/*****/
void MFLevelDataOps::addConstant(LevelData<MFCellFAB>& a_data,
                                 const Real&           a_constant)
{
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      MFCellFAB& data = a_data[d];
      data += a_constant;
    }
}
/*****/
void MFLevelDataOps::power(LevelData<MFCellFAB>& a_data,
                           const Real&           a_exponent,
                           const int&            a_comp)
{
  CH_TIME("MFLevelDataOps::power");
  CH_assert(a_comp<a_data.nComp());

  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      MFCellFAB& dataMF = a_data[dit()];
      for (int iphase=0; iphase<dataMF.numPhases(); iphase++)
        {
          EBCellFAB& data = dataMF.getPhase(iphase);
          const Box& region = data.getRegion();
          const EBISBox& dataEBISBox = data.getEBISBox();
          const EBGraph& dataEBGraph = dataEBISBox.getEBGraph();

          const IntVectSet dataIVS = IntVectSet(region);
          for (VoFIterator vofit(dataIVS,dataEBGraph); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              const Real& dataIn = data(vof,a_comp);
              data(vof,a_comp) = pow(dataIn,a_exponent);
            }
        }
    }
}
/*****/
void MFLevelDataOps::product(LevelData<MFCellFAB>&       a_result,
                             const LevelData<MFCellFAB>& a_in1,
                             const LevelData<MFCellFAB>& a_in2)
{

  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      MFCellFAB& result = a_result[d];
      const MFCellFAB& in1 = a_in1[d];
      const MFCellFAB& in2 = a_in2[d];
      result.copy(in1);
      result *= in2;
    }
}
/*****/
void MFLevelDataOps::product(LevelData<MFCellFAB>&       a_result,
                             const LevelData<MFCellFAB>& a_in1,
                             const LevelData<MFCellFAB>& a_in2,
                             const int&                  a_rComp,
                             const int&                  a_1Comp,
                             const int&                  a_2Comp)
{
  const DisjointBoxLayout& dbl = a_result.disjointBoxLayout();
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      const Box& dblBox = dbl.get(dit());

      MFCellFAB& result = a_result[d];
      const MFCellFAB& in1 = a_in1[d];
      const MFCellFAB& in2 = a_in2[d];

      const Box& regionr = dblBox;
      const Box& region1 = dblBox;

      const Interval intR(a_rComp,a_rComp);
      const Interval int1(a_1Comp,a_1Comp);
      // const Interval int2(a_2Comp,a_2Comp);

      result.copy(region1,intR,regionr,in1,int1);
      result.mult(in2,a_2Comp,a_rComp,1);
    }
}
/*****/
void MFLevelDataOps::divideVectorByScalar(LevelData<MFCellFAB>&       a_vectorOut,
                                          const LevelData<MFCellFAB>& a_vectorIn,
                                          const LevelData<MFCellFAB>& a_scalar)
{
  //a_vector /= a_scalar
  int nCompVectorOut = a_vectorOut.nComp();
  int nCompVectorIn  = a_vectorIn.nComp();
  int nCompScalar = a_scalar.nComp();
  CH_assert(nCompScalar==1);
  CH_assert(nCompVectorOut==nCompVectorIn);

  for (int icomp = 0;icomp<nCompVectorOut;++icomp)
    {
      divide(a_vectorOut,a_vectorIn,a_scalar,icomp,icomp,0);
    }
}
/*****/
void MFLevelDataOps::vectorMagnitude(LevelData<MFCellFAB>&       a_scalarOut,
                                     const LevelData<MFCellFAB>& a_vectorIn,
                                     const int&                  a_pval)
{
  int ncomp  = a_vectorIn.nComp();

  setVal(a_scalarOut, 0.0);

  for (DataIterator dit = a_scalarOut.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      MFCellFAB& scalarOut = a_scalarOut[d];
      const MFCellFAB& vectorIn = a_vectorIn[d];

      for (int iphase=0; iphase<scalarOut.numPhases(); iphase++)
        {
          const EBCellFAB& src = vectorIn.getPhase(iphase);
          EBCellFAB& dst = scalarOut.getPhase(iphase);
          const Box& region = dst.getRegion();
          const EBISBox& ebisbox = dst.getEBISBox();
          const EBGraph& ebgraph = ebisbox.getEBGraph();

          const IntVectSet ivs = IntVectSet(region);
          for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              if (a_pval==0)
                {
                  for (int icomp=0; icomp<ncomp; ++icomp)
                    {
                      dst(vof,0) = Max(dst(vof,0),
                                       Abs(src(vof, icomp)));
                    }
                }
              else if (a_pval==1)
                {
                  for (int icomp=0; icomp<ncomp; ++icomp)
                    {
                      dst(vof,0) += Abs(src(vof, icomp));
                    }
                }
              else
                {
                  Real invp = 1./a_pval;
                  for (int icomp=0; icomp<ncomp; ++icomp)
                    {
                      dst(vof,0) += pow(Abs(src(vof, icomp)),a_pval);
                    }
                  dst(vof,0) = pow(dst(vof,0), invp);
                }
            }
        }
    }
}
/*****/
void MFLevelDataOps::divide(LevelData<MFCellFAB>&       a_result,
                            const LevelData<MFCellFAB>& a_in1,
                            const LevelData<MFCellFAB>& a_in2)
{

  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      MFCellFAB& result = a_result[d];
      const MFCellFAB& in1 = a_in1[d];
      const MFCellFAB& in2 = a_in2[d];

      result.copy(in1);
      result /= in2;
    }
}
/*****/
void MFLevelDataOps::divide(LevelData<MFCellFAB>&       a_result,
                            const LevelData<MFCellFAB>& a_in1,
                            const LevelData<MFCellFAB>& a_in2,
                            const int&                  a_rComp,
                            const int&                  a_1Comp,
                            const int&                  a_2Comp)
{
  const DisjointBoxLayout& dbl = a_result.disjointBoxLayout();
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      const Box& dblBox = dbl.get(dit());

      MFCellFAB& result = a_result[d];
      const MFCellFAB& in1 = a_in1[d];
      const MFCellFAB& in2 = a_in2[d];

      const Box& regionr = dblBox;
      const Box& region1 = dblBox;

      const Interval intR(a_rComp,a_rComp);
      const Interval int1(a_1Comp,a_1Comp);
      // const Interval int2(a_2Comp,a_2Comp);

      result.copy(region1,intR,regionr,in1,int1);
      result.divide(in2,a_2Comp,a_rComp,1);
    }
}
/*****/
void MFLevelDataOps::invert(LevelData<MFCellFAB>&       a_result,
                            const LevelData<MFCellFAB>& a_in1)
{
  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();

      MFCellFAB& result = a_result[d];
      const MFCellFAB& in1 = a_in1[d];
      result.setVal(1.0);
      result /= in1;
    }
}
/*****/
void MFLevelDataOps::kappaWeight(LevelData<MFCellFAB>& a_data)
{
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      MFCellFAB& data = a_data[d];

      kappaWeight(data);
    }
}
/*****/
void MFLevelDataOps::kappaWeight(MFCellFAB& a_data)
{
  for (int iphase=0; iphase<a_data.numPhases(); iphase++)
    {
      int nComp = a_data.nComp(iphase);
      EBCellFAB& data = a_data.getPhase(iphase);
      const Box& dataBox = data.box();
      const EBISBox& dataEBISBox = data.getEBISBox();
      const EBGraph& dataEBGraph = dataEBISBox.getEBGraph();

      const IntVectSet& dataIrreg = dataEBISBox.getIrregIVS(dataBox);

      for (VoFIterator vofit(dataIrreg,dataEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          for (int comp = 0; comp < nComp; comp++)
            {
              data(vof,comp) *= dataEBISBox.volFrac(vof);
            }
        }
    }
}
/*****/
void MFLevelDataOps::kappaScale(LevelData<MFCellFAB>& a_data,
                                const Real&           a_scale)
{
  kappaWeight(a_data);
  scale(a_data,a_scale);
}

#include "NamespaceFooter.H"
