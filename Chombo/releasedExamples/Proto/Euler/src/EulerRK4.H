#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EULERRK4_H_
#define _EULERRK4_H_

#include "EulerOp.H"
#include "Proto.H"
#include "LevelData.H"
#include "FArrayBox.H"


class EulerDX;
using std::shared_ptr;
///
class EulerState
{
public:
  EulerState(){};
  ~EulerState(){};
  EulerState(shared_ptr<LevelData<FArrayBox> > a_U);
  void increment(const EulerDX& a_DX);
  //for max wave speeds
  double m_velSave;
  shared_ptr<LevelData<FArrayBox> > m_U;
  DisjointBoxLayout m_grids;
  Copier m_exchangeCopier;
};

///
class EulerDX
{
public:
  EulerDX(){};
  ~EulerDX(){};
  shared_ptr<LevelData<FArrayBox> > m_DU;
  DisjointBoxLayout m_grids;

  void increment(double& a_weight,const EulerDX& a_incr);
  void init(const EulerState& a_State);
  void operator*=(const double& a_weight);

};

///
class EulerRK4Op
{

public:
  void operator()(EulerDX& a_DX,
                  double a_time,
                  double a_dt,
                  EulerState& a_State) const;
  
  static double maxWave(EulerState& a_State);
};

#endif
