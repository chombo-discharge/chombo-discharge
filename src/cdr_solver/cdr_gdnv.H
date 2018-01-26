/*!
  @file cdr_gdnv.H
  @brief Godunov solver for CDR equations
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _CDR_GDNV_
#define _CDR_GDNV_

#include "cdr_solver.H"

/*!
  @brief Abstract class for cdr-equations. This is (almost) a whole solver. 
*/
class cdr_gdnv : public cdr_solver {
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
    @brief Query number of ghost cells
  */
  virtual int query_ghost() const;

  /*!
    @brief Allocate internal stuff
  */
  virtual void allocate_internals();

  /*!
    @brief Set the nonconservative divergence method. Use = 0 for div(F)_nc = sum(divF(c))/norm, 1 for covered face extrapolation.
  */
  virtual void set_divF_nc(const int a_which_divFnc);
  
protected:

  /*!
    @brief Allocate covered stuff
  */
  virtual void allocate_covered();

  /*!
    @brief Extrapolate velocity to covered faces
  */
  virtual void extrapolate_vel_to_covered_faces();

  /*!
    @brief Godunov face extrapolation method
  */
  virtual void extrapolate_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state, const Real a_extrap_dt);

  /*!
    @brief Non-conservative divergence
  */
  virtual void nonconservative_divergence(EBAMRIVData& a_div_nc, const EBAMRCellData& a_divF, const EBAMRFluxData& a_face_state);

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
    @brief Switch between Trebotich and Colella nonconservative divergence
  */
  int m_which_divFnc;
};

#endif
