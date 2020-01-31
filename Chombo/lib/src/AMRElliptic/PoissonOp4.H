#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POISSONOP4_H_
#define _POISSONOP4_H_

#include "AMRMultiGrid.H"
#include "REAL.H"
#include "Box.H"
#include "LevelDataOps.H"
#include "BCFunc.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "NamespaceHeader.H"

class PoissonOp4 : public MGLevelOp<LevelData<FArrayBox> >
{
public:
  PoissonOp4()
  {
  }

  virtual ~PoissonOp4()
  {
  }

  static const int m_nGhost = 2;

  /**
     \name LinearOp functions */
  /*@{*/

  void define(const RealVect& a_dx,
              const ProblemDomain& a_domain,
              BCFunc a_bc);

  virtual void createCoarsened(LevelData<FArrayBox>& a_lhs,
			       const LevelData<FArrayBox>& a_rhs,
			       const int& a_refRat);

  virtual void residual(LevelData<FArrayBox>& a_lhs,
			const LevelData<FArrayBox>& a_phi,
			const LevelData<FArrayBox>& a_rhs,
			bool a_homogeneous = false);

  virtual void preCond(LevelData<FArrayBox>& a_correction,
		       const LevelData<FArrayBox>& a_residual);

  virtual void applyOp(LevelData<FArrayBox>& a_lhs,
		       const LevelData<FArrayBox>& a_phi,
		       bool a_homogeneous = false);

  virtual void create(LevelData<FArrayBox>& a_lhs,
		      const LevelData<FArrayBox>& a_rhs);

  virtual void assign(LevelData<FArrayBox>& a_lhs,
		      const LevelData<FArrayBox>& a_rhs) ;

  virtual Real dotProduct(const LevelData<FArrayBox>& a_1,
                          const LevelData<FArrayBox>& a_2) ;

  virtual void incr(LevelData<FArrayBox>& a_lhs,
		    const LevelData<FArrayBox>& a_x,
		    Real a_scale) ;

  virtual void axby(LevelData<FArrayBox>& a_lhs, 
		    const LevelData<FArrayBox>& a_x,
		    const LevelData<FArrayBox>& a_y,
		    Real a, Real b) ;

  virtual void scale(LevelData<FArrayBox>& a_lhs, 
		     const Real& a_scale) ;

  virtual Real norm(const LevelData<FArrayBox>& a_x, int a_ord);

  virtual void setToZero(LevelData<FArrayBox>& a_x);

  /*@}*/

  /**
     \name MGLevelOp functions */
  /*@{*/

  virtual void relax(LevelData<FArrayBox>& a_e,
                     const LevelData<FArrayBox>& a_residual,
                     int iterations);

  virtual void createCoarser(LevelData<FArrayBox>& a_coarse,
                             const LevelData<FArrayBox>& a_fine,
                             bool ghost);
  /**
     calculate restricted residual
     a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  */
  virtual void restrictResidual(LevelData<FArrayBox>& a_resCoarse,
                                LevelData<FArrayBox>& a_phiFine,
                                const LevelData<FArrayBox>& a_rhsFine);

  /**
     correct the fine solution based on coarse correction
     a_phiThisLevel += I[2h->h](a_correctCoarse)
  */
  virtual void prolongIncrement(LevelData<FArrayBox>& a_phiThisLevel,
                                const LevelData<FArrayBox>& a_correctCoarse);

  /*@}*/


protected:

  RealVect                m_dx;
  RealVect                m_dxCrse;
  ProblemDomain           m_domain;
  BCFunc                  m_bc;

  LevelDataOps<FArrayBox> m_levelOps;
  DisjointBoxLayout       m_coarsenedMGrids;
  
  void levelGSRB(LevelData<FArrayBox>& a_e,  
		 const LevelData<FArrayBox>& a_residual);

};

class PoissonOp4Factory : public MGLevelOpFactory<LevelData<FArrayBox> >
{
public:

  PoissonOp4Factory();

  PoissonOp4Factory(RealVect& a_dx, BCFunc a_bc);

  void define(const RealVect& m_dx, BCFunc m_bc);

  virtual ~PoissonOp4Factory()
  {
  }

  virtual PoissonOp4* MGnewOp(const ProblemDomain& a_FineindexSpace,
				     int   a_depth,
				     bool  a_homoOnly = true);

  virtual void MGreclaim(PoissonOp4* a_reclaim);

  RealVect                m_dx;
  BCFunc                  m_bc;

};

#include "NamespaceFooter.H"
#endif
