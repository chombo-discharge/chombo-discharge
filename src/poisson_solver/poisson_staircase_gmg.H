/*!
  @file poisson_staircase_gmg.H
  @brief Geometric multigrid poisson solver
  @details This solver uses a staircasing approximation for dielectrics
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _POISSON_STAIRCASE_GMG_
#define _POISSON_STAIRCASE_GMG_

#include "poisson_solver.H"

#include <AMRMultiGrid.H>


/*!
  @brief Abstract Poisson solver class. Just an interface, so pretty lightweight stuff.
*/
class poisson_staircase_gmg : public poisson_solver {
public:

  /*!
    @brief Constructor
  */
  poisson_staircase_gmg();

  /*!
    @brief Constructor
  */
  virtual ~poisson_staircase_gmg();

  /*!
    @brief Override this one for testing purposes. 
  */
  virtual void solve();

  /*!
    @brief Solve Poisson onto state 
  */
  virtual void solve(MFAMRCellData& a_state, const MFAMRCellData& a_source);

  /*!
    @brief 2 ghost cells. That's what we need. 
  */
  virtual int query_ghost() const;
  
protected:
  
  /*!
    @brief Bottom solver
  */
  RefCountedPtr<LinearSolver<LevelData<EBCellFAB> > > m_bottom_solver;

  /*!
    @brief a coefficient
  */
  Vector<RefCountedPtr<LevelData<EBCellFAB> > > m_aco;

  /*!
    @brief b-coefficient on face centers
  */
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > > m_bco;

  /*!
    @brief b-coefficient on eb centroids
  */
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >m_bco_irreg;

  /*!
    @brief Geometric multigrid solver
  */
  AMRMultiGrid<LevelData<EBCellFAB> > m_gmg_solver;

  /*!
    @brief Needs setup
  */
  bool m_needs_setup;

  /*!
    @brief Set up the geometric multigrid solver 
  */
  virtual void setup_gmg();

};

#endif
