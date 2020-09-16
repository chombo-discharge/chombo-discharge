#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EulerRK4.H"
#include "DataIterator.H"
#include "ProtoInterface.H"

using   Proto::Point;
typedef Proto::Box Bx;
using ProtoCh::getProtoBox;
using ProtoCh::getPoint;
using ProtoCh::aliasBoxData;

/****/
EulerState::
EulerState(shared_ptr<LevelData<FArrayBox> > a_U)
{
  m_U     = a_U;
  m_grids = m_U->disjointBoxLayout();
}
/****/
void 
EulerState::
increment(const EulerDX & a_DX)
{
  CH_TIME("EulerState::increment");
  LevelData<FArrayBox> & data  = *m_U;
  LevelData<FArrayBox> & delta = *(a_DX.m_DU);
  
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    data[dit[ibox]] += delta[dit[ibox]];
  }
}
/****/
void 
EulerDX::
increment(double        & a_weight,
          const EulerDX & a_DX)
{
  CH_TIME("EulerDX::increment");
  LevelData<FArrayBox>& data  = *m_DU;
  LevelData<FArrayBox>& delta = *(a_DX.m_DU);
  
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    data[dit[ibox]].plus(delta[dit[ibox]], a_weight);
  }
}

/****/
void 
EulerDX::
init(const EulerState& a_State)
{
  CH_TIME("EulerDX::init");
  m_grids = a_State.m_grids;
  m_DU = shared_ptr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grids, a_State.m_U->nComp(), a_State.m_U->ghostVect()));

  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    (*m_DU)[dit[ibox]].setVal(0.);
  }
}

/****/
void 
EulerDX::
operator*=(const double& a_weight)
{
  CH_TIME("EulerDX::operator*=");
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    (*m_DU)[dit[ibox]] *= a_weight;
  }
}

/****/
void 
EulerRK4Op::
operator()(EulerDX& a_DX,
           double a_time,
           double a_dt,
           EulerState& a_State) const
{
  CH_TIMERS("EulerRKOp::operator()");
  CH_TIMER("defining leveldatas",tdef);
  CH_TIMER("copying to temporary",tcop);
  CH_TIMER("RK arithmetic_no_comm",  trk);

  CH_START(tdef);
  int ncomp  =  a_State.m_U->nComp();
  IntVect gv =  a_State.m_U->ghostVect();
  DisjointBoxLayout grids = a_State.m_grids;
  LevelData<FArrayBox> U_ave(grids, ncomp, gv);
  LevelData<FArrayBox>& delta = *(a_DX.m_DU);
  CH_STOP(tdef);

  CH_START(tcop);
  Interval interv(0, ncomp-1);
  a_State.m_U->copyTo(interv, U_ave, interv);
  CH_STOP(tcop);


  CH_START(trk);
  DataIterator dit = grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    U_ave[dit[ibox]] += delta[dit[ibox]];
  }

  double velmax = EulerOp::step(*a_DX.m_DU, U_ave);
  a_State.m_velSave = std::max(a_State.m_velSave,velmax);

#pragma omp parallel
  for(int ibox = 0; ibox <  dit.size(); ibox++)
  {
    delta[dit[ibox]] *= a_dt;
  }
  CH_STOP(trk);
}


double
EulerRK4Op::
maxWave(EulerState& a_State)
{
  CH_TIME("EulerRKOp::maxwave_init");
  EulerDX DX;
  DX.init(a_State);
  double velmax = EulerOp::step(*(DX.m_DU),*(a_State.m_U));
  return velmax;
}
