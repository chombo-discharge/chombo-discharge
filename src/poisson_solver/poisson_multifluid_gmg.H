/*!
  @file poisson_multifluid_gmg.H
  @brief Geometric multigrid multifluid poisson solver
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _POISSON_MULTIFLUID_GMG_
#define _POISSON_MULTIFLUID_GMG_

#include "poisson_solver.H"
#include "mf_helmholtz_opfactory.H"

/*!
  @brief Abstract Poisson solver class. Just an interface, so pretty lightweight stuff.
*/
class poisson_multifluid_gmg : public poisson_solver {
public:

  /*!
    @brief Constructor
  */
  poisson_multifluid_gmg();

  /*!
    @brief Constructor
  */
  virtual ~poisson_multifluid_gmg();

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
    @brief Needs setup
  */
  bool m_needs_setup;

  /*!
    @brief Operator factory
  */
  RefCountedPtr<mf_helmholtz_opfactory> m_opfact;

  /*!
    @brief Set up the geometric multigrid solver 
  */
  virtual void setup_gmg();
};
#endif
