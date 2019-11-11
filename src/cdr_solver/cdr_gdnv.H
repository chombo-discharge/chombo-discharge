/*!
  @file cdr_gdnv.H
  @brief Godunov solver for CDR equations
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _CDR_GDNV_
#define _CDR_GDNV_

#include "cdr_tga.H"

#include <EBAdvectLevelIntegrator.H>

/*!
  @brief Godunov implementation for advection. 
*/
class cdr_gdnv : public cdr_tga {
public:

  /*!
    @brief Constructor
  */
  cdr_gdnv();

  /*!
    @brief Destructor
  */
  virtual ~cdr_gdnv();

  /*!
    @brief Parse options
  */
  virtual void parse_options();

  /*!
    @brief Query number of ghost cells
  */
  virtual int query_ghost() const;

  /*!
    @brief Allocate internal stuff
  */
  virtual void allocate_internals();

  /*!
    @brief Regrid internals
  */
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Average velocities to faces
  */
  virtual void average_velo_to_faces();

  /*!
    @brief Compute div(nv - D*grad(n)) explicitly
    @note Overriden because cdr_gdnv has other ways of computing the nonconservative divergence
  */
  virtual void compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt);

  /*!
    @brief Compute divF
    @note This function is implemented because cdr_gdnv has several ways of computing the nonconservative divergence
  */
  virtual void compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt);
  
protected:

  /*!
    @brief Simple switch for nonconservative average computation
  */
  enum nonConsDiv {
    conservative_average = 0,
    covered_face = 1,
  };

  /*!
    @brief Allocate covered stuff
  */
  virtual void allocate_covered();

  /*!
    @brief Call parent function: Then fix no influx
  */
  virtual void average_velo_to_faces(EBAMRFluxData& a_velo_face, const EBAMRCellData& a_velo_cell);

  /*!
    @brief Delete covered storage
  */
  virtual void delete_covered();

  /*!
    @brief Extrapolate velocity to covered faces
  */
  virtual void extrapolate_vel_to_covered_faces();

  /*!
    @brief Godunov face extrapolation method for advection
  */
  virtual void advect_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state, const Real a_extrap_dt);

  /*!
    @brief Non-conservative divergence. Overwriting because this class has two ways of doing this. 
  */
  virtual void nonconservative_divergence(EBAMRIVData& a_div, const EBAMRCellData& a_divF, const EBAMRFluxData& a_face_state);

  /*!
    @brief Covered faces on side Lo
  */
  Vector<RefCountedPtr<LayoutData<Vector<Vector<VolIndex> > > > > m_covered_face_lo;

  /*!
    @brief Covered faces on side Hi
  */
  Vector<RefCountedPtr<LayoutData<Vector<Vector<VolIndex> > > > > m_covered_face_hi;

  /*!
    @brief Holder for covered data. This is ugly, but the final Vector<IntVectSet> represents a SpaceDim vector holding
    the grid indices of the VolIndexes which has a covered face on the Side::Lo side in direction dir=0,SpaceDim-1
  */
  Vector<RefCountedPtr<LayoutData<Vector<IntVectSet> > > > m_covered_sets_lo;

  /*!
    @brief Holder for covered data. This is ugly, but the final Vector<IntVectSet> represents a SpaceDim vector holding
    the grid indices of the VolIndexes which has a covered face on the Side::Hi side in direction dir=0,SpaceDim-1
  */
  Vector<RefCountedPtr<LayoutData<Vector<IntVectSet> > > > m_covered_sets_hi;

  /*!
    @brief State on covered faces on side Lo 
  */
  Vector<RefCountedPtr<LayoutData<Vector<BaseIVFAB<Real>* > > > > m_covered_phi_lo;

  /*!
    @brief State on covered faces on side Hi
  */
  Vector<RefCountedPtr<LayoutData<Vector<BaseIVFAB<Real>* > > > > m_covered_phi_hi;

  /*!
    @brief Advective velocity on covered faces on side Lo 
  */
  Vector<RefCountedPtr<LayoutData<Vector<BaseIVFAB<Real>* > > > > m_covered_velo_lo;

  /*!
    @brief Advective velocity on covered faces on side Hi
  */
  Vector<RefCountedPtr<LayoutData<Vector<BaseIVFAB<Real>* > > > > m_covered_velo_hi;

  /*!
    @brief Level advect integrator
  */
  Vector<RefCountedPtr<EBAdvectLevelIntegrator> > m_level_advect;

  /*!
    @brief Switch between nonconservative divergence methods
  */
  nonConsDiv m_nonConsDiv;

  /*!
    @brief Turn on/off slope limiting
  */
  bool m_slopelim;

  /*!
    @brief Turn on/off source terms when time-extrapolating
  */
  bool m_extrap_source;

  /*!
    @brief Parse options for setting the hybrid divergence
  */
  virtual void parse_hybrid_div();   

  /*!
    @brief Parses slope limiter options
  */
  virtual void parse_slopelim();     // Parses slope limiter settings

  /*!
    ®brief Parse source term stuff
  */
  virtual void parse_extrap_source();
};

#endif
