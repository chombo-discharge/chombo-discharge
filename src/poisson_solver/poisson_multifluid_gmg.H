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
    @brief b-coefficient
  */
  MFAMRCellData m_aco;

  /*!
    @brief b-coefficient
  */
  MFAMRFluxData m_bco;

  /*!
    @brief b-coefficient
  */
  MFAMRIVData m_bco_irreg;

  /*!
    @brief Set up the geometric multigrid solver 
  */
  virtual void setup_gmg();

  /*!
    @brief Set coefficients
  */
  virtual void set_coefficients();

  /*!
    @brief Set the permittivities
  */
  virtual void set_permittivities(const Vector<dielectric>& a_dielectrics);

  /*!
    @brief Set face-centered permittivities
  */
  virtual void set_face_perm(EBFluxFAB&                a_perm,
			     const RealVect&           a_origin,
			     const Real&               a_dx,
			     const Vector<dielectric>& a_dielectrics);

  /*!
    @brief Set EB-centered permittivities
  */
  virtual void set_eb_perm(BaseIVFAB<Real>&          a_perm,
			   const RealVect&           a_origin,
			   const Real&               a_dx,
			   const Vector<dielectric>& a_dielectrics);


  // Putting this in a function of its owns
  virtual void base_tests();
};
#endif
