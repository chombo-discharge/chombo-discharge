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
  @brief Godunov implementation for convection-diffusion-reaction equations. 

  @details Options that are available for this class
  
  Class options
  -------------

      cdr_gdnv.divF_nc      = covered_face    # Valid options are "covered_face" and "conservative_average"
      cdr_gdnv.limit_slopes = false           # Use van Leer slope limiter or not (default is false)

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
    @brief Advective advance
  */
  virtual void advance_advect(EBAMRCellData& a_state, const Real a_dt);

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

  /*!
    @brief Regrid internals
  */
  virtual void regrid(const int a_old_finest_level, const int a_new_finest_level);
  
protected:

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
    @brief Compute div(J)
  */
  virtual void compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt);

  /*!
    @brief Non-conservative divergence. Overwriting because this class has two ways of doing this. 
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
    @brief Level advect integrator
  */
  Vector<RefCountedPtr<EBAdvectLevelIntegrator> > m_level_advect;

  /*!
    @brief Switch between Trebotich and Colella nonconservative divergence
  */
  int m_which_divFnc;

  /*!
    @brief Turn on/off slope limiting
  */
  bool m_slopelim;
};

#endif
