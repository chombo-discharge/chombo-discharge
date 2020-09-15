#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _FASMULTIGRID_H_
#define _FASMULTIGRID_H_

#include "AMRMultiGrid.H" // So we can cast MGLevelOp's as AMRLevelOp's
#include "MultiGrid.H"


#include "NamespaceHeader.H"

template <class T>
class FASMultiGrid : public MultiGrid<T>
{
public:
  FASMultiGrid();
  virtual ~FASMultiGrid();


  /// Function to define a FASMultiGrid object.
  /* Author: Jamie Parkinson, February 2017.
   An extension of AMRFASMultiGrid to even coarser levels*/
  /**
      a_factory is the factory for generating operators.
      a_bottomSolver is called at the bottom of v-cycle.
      a_domain is the problem domain at the top of the vcycle.
      maxDepth defines the location of the bottom of the v-cycle.
      The vycle will terminate (hit bottom) when the factory returns NULL for a paticular
      depth if maxdepth = -1.  Otherwise the vycle terminates at maxdepth.
   */
  virtual void define(MGLevelOpFactory<T>& a_factory,
                      LinearSolver<T>* a_bottomSolver,
                      const ProblemDomain& a_domain,
                      int  a_maxDepth = -1,
                      MGLevelOp<T>* a_finestLevelOp = NULL);

  ///
  /**
      solve L(a_phi) = a_rhs.   Tolerance is how much you want the norm of the error reduced.
      verbosity is how chatty you want the function to be.   maxIterations is the maximum number
      of v-cycles.   This does the whole residual correction switcharoo and calls oneCycle up to
      maxIterations times, evaluating the residual as it goes.
   */
  virtual void solve(T& a_phi, const T& a_rhs, Real a_tolerance, int a_maxIterations, int verbosity= 0);

  /// Execute ONE v-cycle of multigrid.
  /**
       If you want the solution to converge, you need to iterate this.
       See solve() or AMRFASMultiGrid::solve for a more automatic solve() function.
   */
  virtual void oneCycle(T& a_e, const T& a_res);
  void oneCycle(T& a_e, const T& a_res, T* a_phiCoarse);

  //internal function
  virtual void cycle(int a_depth, T& a_correction, const T& a_residual);

protected:

  Vector< T* >   m_phiCoarse;



};


//Constructor
template <class T>
FASMultiGrid<T>::FASMultiGrid()
: MultiGrid<T>()
  {
  }


// Destructor
template <class T>
FASMultiGrid<T>::~FASMultiGrid()
{
  CH_TIME("~FASMultiGrid");
}

template <class T>
void FASMultiGrid<T>::define(MGLevelOpFactory<T>& a_factory,
                             LinearSolver<T>*     a_bottomSolver,
                             const ProblemDomain& a_domain,
                             int                  a_maxDepth,
                             MGLevelOp<T>*        a_finestOp)
                             {
  MultiGrid<T>::define(a_factory,a_bottomSolver,a_domain, a_maxDepth, a_finestOp);

  // Should never be homogeneous in the FAS case, as we're solving the full problem
  // so need domain BCs
  this->m_homogeneous = false;

}


template <class T>
void FASMultiGrid<T>::solve(T& a_phi, const T& a_rhs, Real a_tolerance, int a_maxIterations, int verbosity)
{
  CH_TIME("FASMultiGrid::solve");
  this->init(a_phi, a_rhs);

  T correction, residual;
  this->m_op[0]->create(correction, a_phi);
  this->m_op[0]->create(residual, a_rhs);
  this->m_op[0]->setToZero(a_phi);
  this->m_op[0]->residual(residual, a_phi, a_rhs, false);

  Real errorno = this->m_op[0]->norm(residual, 0);
  if (verbosity > 2)
  {
    pout() << "FASmultigrid::solve initial residual = " << errorno << std::endl;
  }
  Real compval = a_tolerance*errorno;
  Real epsilon = 1.0e-16;
  compval = Max(compval, epsilon);
  Real error = errorno;
  int iter = 0;
  while ((error > compval) && (error > a_tolerance*errorno) && (iter < a_maxIterations))
  {
    this->m_op[0]->setToZero(correction);
    this->m_op[0]->residual(residual, a_phi, a_rhs, false);
    error = this->m_op[0]->norm(residual, 0);
    if (verbosity > 3)
    {
      pout() << "FASMultigrid::solve iter = " << iter <<  ",  residual = " << error << std::endl;
    }

    this->cycle(0, correction, residual);
    this->m_op[0]->incr(a_phi, correction, 1.0);

    iter++;
  }
  if (verbosity > 2)
  {
    pout() << "FASMultigrid::solve final residual = " << error << std::endl;
  }

  this->m_op[0]->clear(correction);
  this->m_op[0]->clear(residual);
}

// Do one cycle
// in FAS we solve the full problem, so instead of having the
// error and residual, we have phi and the rhs

template <class T>
void FASMultiGrid<T>::oneCycle(T& a_phi, const T& a_rhs)
{
  // We should never do this - always need a coarse phi
  CH_assert(0);
}

template <class T>
void FASMultiGrid<T>::oneCycle(T& a_phi, const T& a_rhs, T* a_phiCoarse)
{
  CH_TIME("FASMultigrid::oneCycle");

  // First need to create coarser versions of phiCoarse at all MG depths
  m_phiCoarse.resize(this->m_depth);
  m_phiCoarse[0] = a_phiCoarse;

  for (int depth = 1; depth < this->m_depth; depth++)
  {
    if (m_phiCoarse[depth-1] == NULL)
    {
      m_phiCoarse[depth] = NULL;
    }
    else
    {
      const T& phiCoarseGridFiner = *m_phiCoarse[depth-1];

      m_phiCoarse[depth] = new T();

      T scratchCrse;
      this->m_op[depth]->create(scratchCrse,  phiCoarseGridFiner);
      this->m_op[depth]->setToZero(scratchCrse);

      this->m_op[depth]->createCoarser(*(m_phiCoarse[depth]), phiCoarseGridFiner, true);
      this->m_op[depth]->restrictR(*(m_phiCoarse[depth]), phiCoarseGridFiner);

      // If the coarse level is so coarse that it only has one cell then we can't do
      // CF interpolation with it, so delete it
      const DisjointBoxLayout& dbl =  m_phiCoarse[depth]->disjointBoxLayout();
      const Box& coarseDomBox = dbl.physDomain().domainBox();
      if (coarseDomBox.smallEnd() == coarseDomBox.bigEnd())
      {
        delete m_phiCoarse[depth];
        m_phiCoarse[depth] = NULL;
      }
    }
  }

  // We pass in a_phi and a_rhs here
  this->cycle(0, a_phi, a_rhs);

  //Cleanup
  // Note that m_phiCoarse[0] = a_phiCoarse, so we start the loop from i=1
  for (int i=1; i<m_phiCoarse.size(); i++)
  {
    if (m_phiCoarse[i]!=NULL)
    {
      delete m_phiCoarse[i];
      m_phiCoarse[i] = NULL;
    }
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <class T>
void FASMultiGrid<T>::cycle(int depth, T& a_phi, const T& a_rhs)
{
  CH_TIME("FASMultigrid::cycle");

  // currently I can't drop a timer in this function because it is recursive.  doh
  if (depth == this->m_depth-1)
  {
    CH_TIME("FASMultigrid::cycle:bottom-solve");

    // Store the original phi so we can compute the 'correction' later
    this->m_op[0]->setToZero(*this->m_correction[depth]);
    this->m_op[0]->incr(*this->m_correction[depth], a_phi, -1.0);

    // On the coarsest level, just do loads of relaxing as the bottom solve
    this->m_op[depth  ]->relaxNF(a_phi, m_phiCoarse[depth], a_rhs, this->m_bottom);

    // Subtract off new phi to get actual correction from this level
    this->m_op[0]->incr(*this->m_correction[depth], a_phi, 1.0);


  }
  else
  {
    int cycles = this->m_cycle;
    if ( cycles < 0 )
    {
      MayDay::Error("FASMultigrid::cycle - non-V cycle options not implemented");
    }
    else
    {
      // Create some temporary data structures
      T phi_coarse, op_coarse;

      this->m_op[depth+1]->create(phi_coarse,  *(this->m_correction[depth+1]));
      this->m_op[depth+1]->setToZero(phi_coarse);


      this->m_op[depth+1]->create(op_coarse,  *(this->m_correction[depth+1]));
      this->m_op[depth+1]->setToZero(op_coarse);

      // Store initial guess at phi on this level
      // Note this is negative, as the eventual correction = u_new - u_old (and this is u_old)
      this->m_op[depth]->setToZero( *(this->m_correction[depth]));
      this->m_op[depth]->incr(*(this->m_correction[depth]), a_phi, -1.0);

      // Relax phi towards actual phi
      this->m_op[depth  ]->relaxNF( a_phi, m_phiCoarse[depth], a_rhs, this->m_pre );


      // restrict this new phi to a coarser grid and store it in 'm_correction'

      T scratch;
      this->m_op[depth]->create(scratch,  a_phi);
      this->m_op[depth]->setToZero(scratch);

      this->m_op[depth]->createCoarser(phi_coarse, a_phi, true); // create coarser grids for phi

      this->m_op[depth  ]->restrictR(phi_coarse, a_phi);

      // Construct modified RHS

      // Calculate Op(Restrict(phi_fine))
      this->m_op[depth+1]->applyOpMg(op_coarse, phi_coarse, m_phiCoarse[depth+1], false); // false - not homogeneous

      // This sets rhs = Restrict(residual fine)
      this->m_op[depth  ]->restrictResidual(*(this->m_residual[depth+1]), a_phi, m_phiCoarse[depth], a_rhs, false); // false - not homogeneous

      // This sets rhs = Restrict(residual fine) + Op(Restrict(phi fine))
      this->m_op[depth+1]->incr(*(this->m_residual[depth+1]), op_coarse, 1.0);

      for (int img = 0; img < cycles; img++)
      {
        cycle(depth+1, phi_coarse, *(this->m_residual[depth+1]));
      }

      // Interpolate the correction from the coarser level and add to phi on this level
      this->m_op[depth  ]->prolongIncrement(a_phi, *(this->m_correction[depth+1]));

      // Do some more smoothing of phi on this level
      this->m_op[depth  ]->relaxNF(a_phi, m_phiCoarse[depth], a_rhs, this->m_post);

      // Finally, compute the 'correction' (initial phi on this level - current phi on this level)
      // so we can interpolate it when we get back to  finer levels
      this->m_op[depth]->incr(*(this->m_correction[depth]), a_phi, 1.0);

      // Clear temporary data objects
      this->m_op[depth+1]->clear(phi_coarse);
      this->m_op[depth+1]->clear(op_coarse);
    }
  }
}


#include "NamespaceFooter.H"
#endif /*_FASMULTIGRID_H_*/
