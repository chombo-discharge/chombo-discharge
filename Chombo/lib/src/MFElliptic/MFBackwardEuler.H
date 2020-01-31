#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFBACKWARDEULER_H_
#define _MFBACKWARDEULER_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "AMRMultiGrid.H"
#include "LevelData.H"
#include "MFCellFAB.H"
#include "AMRTGA.H"
#include "NamespaceHeader.H"

///
/**
   Solves dphi/dt = L phi + rho
   using backward euler.
   Uses TGAHelmOp functionaltity.
 **/
class MFBackwardEuler
{
public:


  ///
  ~MFBackwardEuler();

  ///
  /**
   **/
  MFBackwardEuler(const RefCountedPtr<AMRMultiGrid< LevelData<MFCellFAB> > > & a_solver,
                  const AMRLevelOpFactory<LevelData<MFCellFAB> >&              a_factory,
                  const ProblemDomain&                                         a_level0Domain,
                  const Vector<int>&                                           a_refRat,
                  int a_numLevels = -1, int a_verbosity = 3);


  ///
  /**
     This advances a parabolic pde from a_phiOld to a_phiNew using TGA on a non-moving domain with source term a_source
  **/
  void oneStep(Vector<LevelData<MFCellFAB>*>&             a_phiNew,
               Vector<LevelData<MFCellFAB>*>&             a_phiOld,
               Vector<LevelData<MFCellFAB>*>&             a_source,
               const Real&             a_dt,
               int                     a_lbase,
               int                     a_lmax,
               Real                    a_told = 0,
               bool                    a_zeroPhi = true);

  ///
  void residual(Vector<LevelData<MFCellFAB>* >&       a_error,
                Vector<LevelData<MFCellFAB>* >&       a_phiNew,
                Vector<LevelData<MFCellFAB>* >&       a_phiOld,
                Vector<LevelData<MFCellFAB>* >&       a_source,
                const Real&                           a_dt,
                int                                   a_lbase,
                int                                   a_lmax,
                Real                                  a_told);

  ///
  void resetAlphaAndBeta(const Real& a_alpha,
                         const Real& a_beta);

  ///
  void setTime(Real a_time);

protected:
  void solveHelm(Vector<LevelData<MFCellFAB>* >&       a_ans,
                 Vector<LevelData<MFCellFAB>* >&       a_rhs,
                 int               a_lbase,
                 int               a_lmax,
                 Real              a_dt,
                 bool              a_zeroPhi);

  void applyHelm(Vector<LevelData<MFCellFAB>*>&       a_ans,
                 Vector<LevelData<MFCellFAB>*>&       a_phi,
                 int                                  a_lbase,
                 int                                  a_lmax,
                 Real                                 a_mu,
                 Real                                 a_dt,
                 bool                                 a_homogeneousBC);

  //fills a_ans = dt*kappa*a_source + dt*phiOld
  void createEulerRHS(Vector<LevelData<MFCellFAB>* >&       a_ans,
                      Vector<LevelData<MFCellFAB>* >&       a_source,
                      Vector<LevelData<MFCellFAB>* >&       a_phiOld,
                      int               a_lbase,
                      int               a_lmax,
                      Real              a_dt);
  void createData(Vector<LevelData<MFCellFAB>* >&       a_source,
                  int                a_lbase,
                  int                a_lmax);

  TGAHelmOp<LevelData<MFCellFAB> >*
  newOp(const ProblemDomain&                             a_indexSpace,
        const AMRLevelOpFactory<LevelData<MFCellFAB> >&  a_opFact);

private:
  Vector<RefCountedPtr<TGAHelmOp<LevelData<MFCellFAB> > > >  m_ops;
  Vector< LevelData<MFCellFAB>* >                            m_rhst;
  ProblemDomain                                              m_level0Domain;
  Vector<int>                                                m_refRat;
  RefCountedPtr<AMRMultiGrid< LevelData<MFCellFAB> > >       m_solver;
  int m_verbosity, m_numLevels;
  bool                                           m_dataCreated;

  //copy constructor and operator= disallowed for all the usual reasons
  MFBackwardEuler(const MFBackwardEuler& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const MFBackwardEuler& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  /// weak construction is bad.   Ref Counted pointers are your friends.
  MFBackwardEuler()
  {
    MayDay::Error("invalid operator");
  }


};

#include "NamespaceFooter.H"
#endif
