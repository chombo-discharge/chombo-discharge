#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// BVS, June 18, 2003

// We can assume that template class T has null construction.

#ifndef _AMRFASMULTIGRID_H_
#define _AMRFASMULTIGRID_H_

#include "AMRMultiGrid.H"
#include "FASMultiGrid.H"

#include "NamespaceHeader.H"


enum OLD_FASMG_type {FULL=0,VCYCLE=1,FCYCLE=2};

///
/**
   Class to solve elliptic equations using the FAS multigrid
 */
template <class T>
class AMRFASMultiGrid : public AMRMultiGrid<T>
{
public:
  AMRFASMultiGrid();
  virtual ~AMRFASMultiGrid();

  virtual void define(const ProblemDomain& a_coarseDomain,
                      AMRLevelOpFactory<T>& a_factory,
                      LinearSolver<T>* a_bottomSolver,
                      int a_numLevels);

  virtual void solveNoInit(Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                           int l_max, int l_base, bool a_zeroPhi=true,
                           bool forceHomogeneous = false);

  void setCycleType( OLD_FASMG_type a_type )
  {
    m_type = a_type;
  }
  void setAvoidNorms( bool b = true )
  {
    m_avoid_norms = b;
  }
  void setNumVcycles( int n )
  {
    m_numVcycles = n;
  }

  // Number of multigrid iterations taken
  int m_iter;

private:
  virtual void FMG(Vector<T*>& a_phi,
                   const Vector<T*>& a_rhs,
                   int l, int l_max, int l_base);
  
  virtual void VCycle(Vector<T*>& a_phi,
                      const Vector<T*>& a_rhs,
                      int l, int l_max, int l_base);

  virtual void init(const Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                    int l_max, int l_base);

  void clear_private();

  // data
  OLD_FASMG_type m_type;
  bool m_avoid_norms;  // flag to avoid norms and residuals (for convergence checking)
  int m_numVcycles;
  Vector<Copier> m_HOCopier;
  Vector<Copier> m_revHOCopier;

  Vector<Copier> m_HOCornerCopier;
  ProblemDomain m_coarseDomain; // need to cache this
};

//*******************************************************
// AMRFASMultigrid Implementation
//*******************************************************

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::AMRFASMultiGrid
////////////////////////////////////////////////////////////////////////
template <class T>
AMRFASMultiGrid<T>::AMRFASMultiGrid() 
  : AMRMultiGrid<T>()
{
  m_numVcycles = 1;
  m_type = VCYCLE;
  m_avoid_norms = false;
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::~AMRFASMultiGrid
////////////////////////////////////////////////////////////////////////
template <class T>
AMRFASMultiGrid<T>::~AMRFASMultiGrid()
{
  CH_TIME("~AMRFASMultiGrid");
  clear_private();
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::define
////////////////////////////////////////////////////////////////////////
template <class T>
void AMRFASMultiGrid<T>::define( const ProblemDomain& a_coarseDomain,
                                 AMRLevelOpFactory<T>& a_factory,
                                 LinearSolver<T>* a_bottomSolver,
                                 int a_maxAMRLevels )
{
  CH_TIME("AMRFASMultiGrid::define");

  AMRMultiGrid<T>::define( a_coarseDomain,a_factory,a_bottomSolver,a_maxAMRLevels);

  m_HOCopier.resize( a_maxAMRLevels );
  m_revHOCopier.resize( a_maxAMRLevels );
  m_HOCornerCopier.resize( a_maxAMRLevels );
  m_coarseDomain = a_coarseDomain;

  // We actually want the MG solver to be FASMultigrid, not the default MultiGrid
  // This is so we can solve down to the very coarsest grid in a nonlinear way
  ProblemDomain current = a_coarseDomain;
  for (int i = 0; i < a_maxAMRLevels; i++)
    {

      // first delete initially defined standard MG object
      if (this->m_mg[i] != NULL)
        {
          delete this->m_mg[i];
          this->m_mg[i] = NULL;
        }

      // Create FASMultigrid object instead
      this->m_mg[i]= new FASMultiGrid<T>();
      this->m_mg[i]->define(a_factory, &this->m_nosolve, current, this->m_maxDepth, this->m_op[i]);

      // Only do this if it will be used (avoiding a reference to invalid
      // and/or unavailable refinement ratios)
      if (i < a_maxAMRLevels-1)
        {
          current.refine(a_factory.refToFiner(current));
        }
    }
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::init
////////////////////////////////////////////////////////////////////////
template <class T>
void AMRFASMultiGrid<T>::init( const Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                               int l_max, int l_base )
{
  CH_TIME("AMRFASMultiGrid::init");

  AMRMultiGrid<T>::init( a_phi, a_rhs, l_max,l_base );

  // set new copiers for high order interp
  // why were we doing this?

  ProblemDomain dom = m_coarseDomain;
  for (int i = l_base+1; i <= l_max; i++)
    {
      AMRLevelOp<T>& op = *(this->m_op[i]);
      int r = op.refToCoarser();
      m_HOCopier[i].define( a_phi[i-1]->disjointBoxLayout(), this->m_resC[i]->disjointBoxLayout(),
                            a_phi[i-1]->ghostVect() );

      m_revHOCopier[i] = m_HOCopier[i];
      m_revHOCopier[i].reverse();


      m_HOCornerCopier[i].define( this->m_resC[i]->disjointBoxLayout(), this->m_resC[i]->disjointBoxLayout(), dom, this->m_resC[i]->ghostVect(), true ); //-- smart corner copier -- needed for AMR
      //m_HOCornerCopier[i].define( this->m_resC[i]->disjointBoxLayout(), this->m_resC[i]->disjointBoxLayout(), this->m_resC[i]->ghostVect() ); -- dumb copier - not needed w/o AMR

      dom = dom.refine(r);
    }
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::clear_private
////////////////////////////////////////////////////////////////////////
template <class T>
void AMRFASMultiGrid<T>::clear_private()
{
  CH_TIME("AMRFASMultiGrid::clear_private");

  for (int i = 0; i < this->m_op.size(); i++)
    {
    }
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::solveNoInit
////////////////////////////////////////////////////////////////////////
template<class T>
void AMRFASMultiGrid<T>::solveNoInit( Vector<T*>& a_phi,
                                      const Vector<T*>& a_rhs,
                                      int l_max, int l_base, bool a_zeroPhi,
                                      bool a_forceHomogeneous)
{
  CH_TIMERS("AMRFASMultiGrid::solveNoInit");
  CH_TIMER("AMRFASMultiGrid::cycle", vtimer);
//  CH_assert(!a_forceHomogeneous); // this is now allowed!
  bool outputIntermediates = false;

  this->setBottomSolver( l_max, l_base );

  CH_assert(l_base <= l_max);
  CH_assert(a_rhs.size() == a_phi.size());

  if (a_zeroPhi)
    for (int ilev = l_base; ilev <=l_max; ++ilev)
      {
        this->m_op[ilev]->setToZero(*(a_phi[ilev]));
      }

  Real initial_rnorm = 0;
  if ( !m_avoid_norms )
    {
      CH_TIME("Initial AMR Residual");
      initial_rnorm = this->computeAMRResidual( a_phi, a_rhs, l_max, l_base );
    }

  if (this->m_convergenceMetric != 0.)
    {
      initial_rnorm = this->m_convergenceMetric;
    }

  // set bottom solver convergence norm and solver tolerance
  this->m_bottomSolver->setConvergenceMetrics(initial_rnorm, this->m_bottomSolverEpsCushion*this->m_eps);

  Real rnorm = initial_rnorm;
  Real norm_last = 2*initial_rnorm;

  int iter=0;
  if (this->m_verbosity >= 2 && !m_avoid_norms )
    {
    // indentation in the output is to ensure that the :: lines up with the output from AMRMultigrid
      pout() << " AMRFASMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm << std::endl;
    }

  bool goNorm = rnorm > this->m_normThresh || m_avoid_norms;       //iterate if norm is not small enough
  bool goRedu = rnorm > this->m_eps*initial_rnorm || m_avoid_norms;//iterate if initial norm is not reduced enough
  bool goIter = iter < this->m_iterMax;                            //iterate if iter < max iteration count
  bool goHang = iter < this->m_imin || rnorm <(1-this->m_hang)*norm_last; //iterate if we didn't hang
  bool goMin = iter < this->m_iterMin ; // iterate if iter < min
  while (goMin || (goIter && goRedu && goHang && goNorm))
    {
      norm_last = rnorm;

      //this generates a correction from the current residual
      CH_START(vtimer);
      if ( m_type == FULL )
        {
          FMG( a_phi, a_rhs, l_max, l_max, l_base);
        }
      else if ( m_type == VCYCLE )
        {
          // Set m_residual[l_max] = a_rhs[l_max]
          this->m_op[l_max]->assignLocal( *(this->m_residual[l_max]), *(a_rhs[l_max]) );
          VCycle( a_phi, a_rhs, l_max, l_max, l_base);
        }
      else 
        {
          MayDay::Error("unknown FAS type");
        }
      CH_STOP(vtimer);
      //increment phi by correction and reset correction to zero
      for (int ilev = l_base; ilev <= l_max; ilev++)
        {
          if (outputIntermediates)
            {
              char strcorname[100];
              sprintf(strcorname, "amrFASmg.phi.iter.%03d", iter);
              string namecor(strcorname);
              this->outputAMR(a_phi, namecor, l_max, l_base);
            }
        }

      // For solvers with accuracy higher than 2nd order
      //  consistency between levels has to be explicitly enforced.
      if (this->m_op[0]->orderOfAccuracy()>2)
        {
          for (int ilev=l_max; ilev>l_base; ilev--)
            {
              this->m_op[ilev]->enforceCFConsistency(*a_phi[ilev-1], *a_phi[ilev]);
            }
        }
      
      // recompute residual
      iter++;
      if ( !m_avoid_norms )
        {
          rnorm = this->computeAMRResidual( a_phi, a_rhs, l_max, l_base );

          if (this->m_verbosity >= 2)
            {
              pout() << " AMRFASMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm;
              if (rnorm > 0.0)
                {
                  pout() << ", rate = " << norm_last/rnorm;
                }
              pout() << std::endl;
            }

          goNorm = rnorm > this->m_normThresh;                        //keep iterating if norm is not small enough
          goRedu = rnorm > this->m_eps*initial_rnorm;                 //keep iterating if initial norm is not reduced enough
          goHang = iter < this->m_imin || rnorm <(1-this->m_hang)*norm_last;//keep iterating if we didn't hang
        }
      goIter = iter < this->m_iterMax;                            //keep iterating if iter < max iteration count
      goMin = iter < this->m_iterMin ; // keep iterating if iter < min
    }

  this->m_exitStatus = int(!goRedu) + int(!goIter)*2 + int(!goHang)*4 + int(!goNorm)*8;
  if (this->m_verbosity >= 2 && !m_avoid_norms)
    {
      pout() << "    AMRFASMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm << std::endl;
    }
  if (this->m_verbosity > 1)
    {
      if (!goIter && goRedu && goNorm) // goRedu=T, goIter=F, goHang=?, goNorm=T
        { // m_exitStatus == 0 + 2 + 0|4 + 0 = 2|6
          pout() << "    AMRFASMultiGrid:: WARNING: Exit because max iteration count exceeded" << std::endl;
        }
      if (!goHang && goRedu && goNorm) // goRedu=T, goIter=?, goHang=F, goNorm=T
        { // m_exitStatus == 0 + 0|2 + 4 + 0 = 4|6
          pout() << "    AMRFASMultiGrid:: WARNING: Exit because of solver hang" << std::endl;
        }
      if (this->m_verbosity > 4)
        {
          pout() << "    AMRFASMultiGrid:: exitStatus = " << this->m_exitStatus << std::endl;
        }
    }

  // Store number or iterations taken so we can access this information
  // from other objects
  this->m_iter = iter;
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::FMG
////////////////////////////////////////////////////////////////////////
template<class T>
void AMRFASMultiGrid<T>::FMG( Vector<T*>& a_phi,
                              const Vector<T*>& a_rhs,
                              int i, int l_max, int l_base)
{
  // coarse grid
  this->m_op[l_base]->assignLocal( *(this->m_residual[l_base]), *(a_rhs[l_base]) );
  this->m_mg[l_base]->oneCycle( *(a_phi[l_base]), *(this->m_residual[l_base]) );

  for ( int ilev = l_base+1 ; ilev <= l_max ; ilev++ )
    {
      // FMG interpolation -- high order
      this->m_op[ilev]->AMRProlongS_2( *(a_phi[ilev]), *(a_phi[ilev-1]),
                                       *(this->m_resC[ilev]), this->m_revHOCopier[ilev],
                                       this->m_HOCornerCopier[ilev], this->m_op[ilev-1] );
      // V-cycle
      this->m_op[ilev]->assignLocal( *(this->m_residual[ilev]), *(a_rhs[ilev]) );
      for (int i=0;i<m_numVcycles;i++)
        {
          VCycle( a_phi, a_rhs, ilev, l_max, l_base );
        }
    }
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::VCycle
//   'm_residual' is really the full RHS ala FAS.
//   'm_correction' is just used as a temp
////////////////////////////////////////////////////////////////////////
//#define FAS_PRINT
template<class T>
void AMRFASMultiGrid<T>::VCycle( Vector<T*>& a_phi,
                                 const Vector<T*>& a_rhs_dummy,
                                 int ilev, int l_max, int l_base )
{
  if (ilev == l_base)
    {
      CH_TIME("AMRFASMultiGrid::VCycle:coarse_grid_solver");
      // exact solver

#ifdef FAS_PRINT
      this->m_op[ilev]->write(a_phi[ilev],"zh_coarsest_phi_0.hdf5");
      this->m_op[ilev]->write(this->m_residual[ilev],"zh_coarsest_rhs.hdf5");
#endif

      // Keep doing FAS solves on the coarser levels
      // Remember that m_residual is just the full rhs with FAS
      // Cast as FASMultiGrid so we can use a different onyCycle
      T* phiCoarsePtr = NULL;
      if (l_base > 0)
        {
          phiCoarsePtr = a_phi[l_base-1];
        }
      FASMultiGrid<T>* fasMg = dynamic_cast<FASMultiGrid<T>*> (this->m_mg[l_base]);

      // This gives us back the full solution (NOT a correction)
      // m_residual[l_base] is the full FAS right hand side (with tau correction)
      fasMg->oneCycle( *(a_phi[l_base]), *(this->m_residual[l_base]), phiCoarsePtr);

#ifdef FAS_PRINT
this->m_op[ilev]->write(a_phi[ilev],"zh_coarsest_phi_1.hdf5");
#endif
    }
  else
    {
      //Note: this else clause is currently untested
      //============= Downsweep ========================

    T* rhs = NULL;

    // on the finest level, the rhs is the full rhs
    // on coarser levels, it is whatever we calculate and put in m_residual
    // need to be careful about this, as the fine level m_residual get's overwritten
    // whenever we compute the AMR residual (e.g. at the end of each V-cycle)
    if (ilev == l_max)
    {
      rhs = a_rhs_dummy[ilev];
    }
    else
    {
      rhs = this->m_residual[ilev];
    }


#ifdef FAS_PRINT
this->m_op[ilev]->write(this->m_residual[ilev],"za_fine_res_0.hdf5");
this->m_op[ilev]->write(a_phi[ilev],"zb_fine_phi_0.hdf5");
#endif
                // Replacing relax with relaxNF
 //     this->m_op[ilev]->relax( *(a_phi[ilev]), *(this->m_residual[ilev]), this->m_pre );
                this->m_op[ilev]->relaxNF( *(a_phi[ilev]), (a_phi[ilev-1]), *rhs, this->m_pre );

#ifdef FAS_PRINT
this->m_op[ilev]->write(this->m_residual[ilev],"za_fine_res_1.hdf5");
this->m_op[ilev]->write(a_phi[ilev],"zb_fine_phi_1.hdf5");
this->m_op[ilev]->residual( *(this->m_correction[ilev]),  *a_phi[ilev], *(this->m_residual[ilev]), false);
this->m_op[ilev]->write(this->m_correction[ilev],"zz_fine_res_0.hdf5");
#endif


     // Compute residual on next coarser level
     //  for the valid region NOT covered by this level.
     // includes reflux corrections
    // do this before we fill a_phi[ilev-1] with junk
     this->computeAMRResidualLevel(this->m_residual,
                                   a_phi,
                                   a_rhs_dummy,
                                   l_max, l_base, ilev-1,
                                   false);


      // Compute the restriction of the solution to the coarser level phiC -- R(u_f) 
     // Need this so that later we can compute L_c(R(u_f))
      // This just restricts a_phi[ilev] to m_resC[ilev] (due to the skip residuals flag)
      this->m_op[ilev]->AMRRestrictS(*(this->m_resC[ilev]),      // output
                                     *(a_phi[ilev]),             // input
                                     *(a_phi[ilev]),             // dummy
                                     *(a_phi[ilev-1]),           // dummy
                                     *(this->m_correction[ilev]),// scratch
                                     true );                     // skip residuals

#ifdef FAS_PRINT
this->m_op[ilev]->write(a_phi[ilev],"zb_fine_phi_2.hdf5");
this->m_op[ilev-1]->write(this->m_resC[ilev],"zb_coarse_phi_1.hdf5");
#endif
      // Overwrite R(u_f) on the valid region of the next coarser level a_phi[ilev-1]
      //this->m_op[ilev-1]->assignCopier( *(a_phi[ilev-1]), *(m_phiC[ilev]), m_phiCopier[ilev] );
      this->m_resC[ilev]->copyTo(this->m_resC[ilev]->interval(), *(a_phi[ilev-1]), a_phi[ilev-1]->interval(), this->m_resCopier[ilev]);


      // Compute the restriction of the _residual_ to the coarser level
      //  for the valid region covered by this level
      // resC -- R(f - L_f(u_f))
      this->m_op[ilev]->AMRRestrictS(*(this->m_resC[ilev]),      // output
                                     *rhs,  // const input
                                     *(a_phi[ilev]),             // const but C-F interp modified
                                     *(a_phi[ilev-1]),           // coarse phi, for C-F interp.
                                     *(this->m_correction[ilev])); // scratch


#ifdef FAS_PRINT
this->m_op[ilev]->write(this->m_resC[ilev],"zd_resC.hdf5");
this->m_op[ilev]->write(a_phi[ilev-1],"ze_coarse_phi_2.hdf5");
this->m_op[ilev]->write(this->m_residual[ilev],"za_fine_res_2.hdf5");
#endif

      // Overwrite residual on the valid region of the next coarser level
      //  with coarsened residual from this level
      //this->m_op[ilev-1]->assignCopier( *(this->m_residual[ilev-1]), *(this->m_resC[ilev]), this->m_resCopier[ilev]);
      this->m_resC[ilev]->copyTo(this->m_resC[ilev]->interval(), *(this->m_residual[ilev-1]), this->m_residual[ilev-1]->interval(), this->m_resCopier[ilev] );

#ifdef FAS_PRINT
this->m_op[ilev-1]->write(this->m_residual[ilev-1],"zf_coarse_res_1.hdf5");
#endif

      // Compute the tau correction on the covered region, L_c (R(u_f)), and put it in m_correction[ilev-1]
       // Do not pass in fine level here - we've already done refluxing
//      this->m_op[ilev-1]->AMROperatorNC( *(this->m_correction[ilev-1]), *(a_phi[ilev]), *(a_phi[ilev-1]), true, this->m_op[ilev]);
      T undefined;
      this->m_op[ilev-1]->AMROperatorNC( *(this->m_correction[ilev-1]), undefined, *(a_phi[ilev-1]), true, this->m_op[ilev]);

#ifdef FAS_PRINT
this->m_op[ilev-1]->write(this->m_correction[ilev-1],"zf_coarse_Lu.hdf5");
this->m_op[ilev]->write(a_phi[ilev],"zb_fine_phi_3.hdf5");
#endif

      // compute the full RHS, R( f - L_f(u_f) ) + L_c( R(u_f) )
      this->m_op[ilev-1]->axby( *(this->m_residual[ilev-1]), *(this->m_residual[ilev-1]), *(this->m_correction[ilev-1]), 1.0, 1.0);

#ifdef FAS_PRINT
this->m_op[ilev-1]->write(this->m_residual[ilev-1],"zf_coarse_res_3.hdf5");
#endif

      {
        // store correction in R_u_f -- R(u_f)
        T R_u_f;
        this->m_op[ilev-1]->create( R_u_f, *a_phi[ilev-1]);     
        this->m_op[ilev-1]->assignLocal( R_u_f, *(a_phi[ilev-1]) ); 
        //============finish Compute residual for the next coarser level======
        for (int img = 0; img < this->m_numMG; img++)
          {
            VCycle( a_phi, a_rhs_dummy, ilev-1, l_max, l_base );
          }

#ifdef FAS_PRINT
this->m_op[ilev-1]->write(a_phi[ilev-1],"zi_coarse_phi_2.hdf5");
#endif

        // subtract off initial solution to get an increment (m_correction[ilev-1] is an increment)
        this->m_op[ilev-1]->axby( *this->m_correction[ilev-1], *a_phi[ilev-1], R_u_f, 1.0, -1.0 );
        this->m_op[ilev-1]->clear( R_u_f );
      }

#ifdef FAS_PRINT
this->m_op[ilev-1]->residual( *(this->m_correction[ilev-1]),  *a_phi[ilev-1], *(this->m_residual[ilev-1]), true);
this->m_op[ilev-1]->write(this->m_correction[ilev-1],"zz_coarse_res.hdf5");
#endif

      //================= Upsweep ======================
      // increment the correction with coarser version

// There's a bug with the copiers and AMRProlongS_2, so doing this the old way for now
//      this->m_op[ilev]->AMRProlongS_2( *(a_phi[ilev]), *(a_phi[ilev-1]),
//                                       *(this->m_resC[ilev]), this->m_HOCopier[ilev],
//                                       this->m_HOCornerCopier[ilev], this->m_op[ilev-1] );


      this->m_op[ilev]->AMRProlongS( *(a_phi[ilev]), *(this->m_correction[ilev-1]),
                                            *(this->m_resC[ilev]),
                                            this->m_reverseCopier[ilev]  );


#ifdef FAS_PRINT 
this->m_op[ilev]->write(a_phi[ilev],"zl_fine_phi_1.hdf5");
this->m_op[ilev]->residual( *(this->m_correction[ilev]),  *a_phi[ilev], *(this->m_residual[ilev]), false);
this->m_op[ilev]->write(this->m_correction[ilev],"zz_fine_res_1.hdf5");
#endif
                // Replacing relax with relaxNF
//      this->m_op[ilev]->relax( *(a_phi[ilev]), *(this->m_residual[ilev]), this->m_post );
                this->m_op[ilev]->relaxNF( *(a_phi[ilev]),  (a_phi[ilev-1]), *rhs, this->m_post );

#ifdef FAS_PRINT
this->m_op[ilev]->write(a_phi[ilev],"zl_fine_phi_2.hdf5");
this->m_op[ilev]->write(this->m_residual[ilev],"za_fine_res_3.hdf5");
this->m_op[ilev]->write(a_rhs_dummy[ilev],"za_rhs.hdf5");
this->m_op[ilev]->residual( *(this->m_correction[ilev]),  *a_phi[ilev], *(this->m_residual[ilev]), false);
this->m_op[ilev]->write(this->m_correction[ilev],"zz_fine_res_2.hdf5");
#endif
    }
}

#include "NamespaceFooter.H"

#endif
