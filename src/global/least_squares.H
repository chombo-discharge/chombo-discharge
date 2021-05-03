/*!
  @file LeastSquares.hpp
  @brief Agglomeration of various routines that use least squares computations
  @author Robert Marskar
  @date May 2016
*/

#ifndef _LEASTSQUARES_
#define _LEASTSQUARES_

#include "MultiIndex.hpp"

#include "Stencils.H"
#include "EBISBox.H"
#include "IntVectSet.H"
#include "RealVect.H"
#include "ProblemDomain.H"
#include "VolIndex.H"
#include "PolyGeom.H"
#include "EBISLayout.H"

/*!
  @brief Agglomeration of various useful stencil computations
*/  
class LeastSquares {

public:

  //
  LeastSquares(){
  }

  //
  ~LeastSquares(){
  }

  /*!
    @brief Tolerance for stencil computations. 
  */
  static Real s_tol;
  
  /*!
    @brief Compute an interpolation stencil that interpolates from cell-centered data to centroid-centered data
    @param[out] a_stencil Stencil
    @param[in]  a_vof Input vof
    @param[in]  a_ebisbox ebisbox
    @param[in]  a_dx Resolution
    @param[in]  a_radius order Interpolation order
    @param[in]  a_radius Radius for least squares problem. Typically a_radius >= a_order. If you put radius = -1, we 
    will first try a_radius = a_order, and then increase until the system has more neighbors than unknowns. 
    @return True if stencil was found. False otherwise. 
  */
  static bool getCenterToCentroidInterpStencil(VoFStencil&     a_stencil,
					       const VolIndex& a_vof,
					       const EBISBox&  a_ebisbox,
					       const Real&     a_dx,
					       const int       a_order = 1,
					       const int       a_radius = 1);

  /*!
    @brief Compute an interpolation stencil that interpolates from cell-centered data to EB-centered data
    @param[out] a_stencil Stencil
    @param[in]  a_vof Input vof
    @param[in]  a_ebisbox ebisbox
    @param[in]  a_dx Resolution
    @param[in]  a_radius order Interpolation order
    @param[in]  a_radius Radius for least squares problem. Typically a_radius >= a_order. If you put radius = -1, we 
    will first try a_radius = a_order, and then increase until the system has more neighbors than unknowns. 
    @return True if stencil was found. False otherwise. 
  */
  static bool getCenterToEbCentroidInterpStencil(VoFStencil&     a_stencil,
						 const VolIndex& a_vof,
						 const EBISBox&  a_ebisbox,
						 const Real&     a_dx,
						 const int       a_order = 1,
						 const int       a_radius = 1);
  /*!
    @brief Compute an interpolation stencil that interpolates from centroid-centered data to cell-centered data
    @param[out] a_stencil Stencil
    @param[in]  a_vof Input vof
    @param[in]  a_ebisbox ebisbox
    @param[in]  a_dx Resolution
    @param[in]  a_order Interpolation order
    @param[in]  a_radius Radius for least squares problem. Typically a_radius >= a_order. If you put radius = -1, we 
    will first try a_radius = a_order, and then increase until the system has more neighbors than unknowns. 
    @return True if stencil was found. False otherwise. 
  */
  static bool getCentroidToCenterInterpStencil(VoFStencil&     a_stencil,
					       const VolIndex& a_vof,
					       const EBISBox&  a_ebisbox,
					       const Real&     a_dx,
					       const int       a_order = 1,
					       const int       a_radius = 1);

    /*!
    @brief Compute a least squares gradient stencil on the EB centroid
    @param[out] a_stencil Stencil
    @param[in]  a_vof Input vof
    @param[in]  a_ebisbox ebisbox
    @param[in]  a_dx Resolution
    @param[in]  a_order Interpolation order
    @param[in]  a_radius Radius for least squares problem. Typically a_radius >= a_order. If you put radius = -1, we 
    will first try a_radius = a_order, and then increase until the system has more neighbors than unknowns. 
    @return True if stencil was found. False otherwise. 
  */
  static bool getEbCentroidGradientStencil(VoFStencil&     a_stencil,
					   const VolIndex& a_vof,
					   const EBISBox&  a_ebisbox,
					   const Real&     a_dx,
					   const int       a_order = 1,
					   const int       a_radius = 1);
  
  /*!
    @brief Build the vector of distances that goes into LeastSquares::getCenterToCentroidInterpStencil
    @param[out] a_displ The displacement vector \f$x_i - x_0\f$ where \f$x_i\f$ are cell centers and \f$x_0\f$ is the cell centroid
    @param[in]  a_vof Input vof
    @param[in]  a_monoVoFs VoFs that can be reached. This should include a_vof
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_dx Resolution
  */
  static void getCenterToCentroidDisplacements(Vector<RealVect>&       a_displ,
					       const VolIndex&         a_vof,
					       const Vector<VolIndex>& a_monoVoFs,
					       const EBISBox&          a_ebisbox,
					       const Real&             a_dx);

  /*!
    @brief Build the vector of distances that goes into LeastSquares::getCenterToEbCentroidInterpStencil
    @param[out] a_displ The displacement vector \f$x_i - x_0\f$ where \f$x_i\f$ are cell centers and \f$x_0\f$ is the EB centroid
    @param[in]  a_vof Input vof
    @param[in]  a_monoVoFs VoFs that can be reached. This should include a_vof
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_dx Resolution
  */
  static void getCenterToEbCentroidDisplacements(Vector<RealVect>&       a_displ,
						 const VolIndex&         a_vof,
						 const Vector<VolIndex>& a_monoVoFs,
						 const EBISBox&          a_ebisbox,
						 const Real&             a_dx);

  /*!
    @brief Build the vector of distances that goes into LeastSquares::getCentroidToCenterInterpStencil
    @param[out] a_displ The displacement vector \f$x_i - x_0\f$ where \f$x_i\f$ are cell centroid and \f$x_0\f$ is the cell center.
    @param[in]  a_vof Inputs vof
    @param[in]  a_monoVoFs VoFs that can be reached. This should include a_vof
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_dx Resolution
  */
  static void getCentroidToCenterDisplacements(Vector<RealVect>&       a_displ,
					       const VolIndex&         a_vof,
					       const Vector<VolIndex>& a_monoVoFs,
					       const EBISBox&          a_ebisbox,
					       const Real&             a_dx);
  
  /*!
    @brief Get the least squares interpolation stencil
    @param[out] a_stencil Stencil
    @param[in]  a_monoVoFs VolIndexes
    @param[in]  a_displ Displacements
    @param[in]  a_Q Interpolation order. 
  */
  static bool getLSqInterpStencil(VoFStencil&             a_stencil,
				  const Vector<VolIndex>& a_monoVoFs,
				  const Vector<RealVect>& a_displ,
				  const int               a_Q);

  /*!
    @brief Get the least squares gradient stencil
    @param[out] a_stencil Gradient stencil (has SpaceDim components)
    @param[in]  a_monoVoFs VolIndexes
    @param[in]  a_displ Displacements
    @param[in]  a_Q Interpolation order. 
  */
  static bool getLSqGradientStencil(VoFStencil&             a_stencil,
				    const Vector<VolIndex>& a_monoVoFs,
				    const Vector<RealVect>& a_displ,
				    const int               a_Q);

  /*!
    @brief Get neighbouring VolIndex. Return false if there are fewer neighbours than a_minVoFs. If a_radius = -1, adjust radius
    until we get enough VoFs.
    @param[out] a_allVofs All VolIndex in a given radius (including a_vof)
    @param[in]  a_vof Input vof
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_minVoFs Minimum amoumt of vofs. Return false if a_allVoFsRad < a_minVoFs
    @param[in]  a_order Taylor expansion order
    @param[in]  a_radius Radius to search. If radius = -1, adjust radius until we get enough vofs

  */
  static bool getVoFsRadius(Vector<VolIndex>& a_allVoFsRad,
			    const VolIndex&   a_vof,
			    const EBISBox&    a_ebisbox,
			    const int         a_minVoFs,
			    const int         a_order, 
			    const int         a_radius = -1);
  
  /*!
    @brief Generate a sequence of multi-indices |q| <= Q in lexicographical order
    @param[in] a_Q Largest multi-index
  */
  static Vector<MultiIndex> getMultiIndicesLexiOrder(const int a_Q);

  /*!
    @brief Get the size of a Taylor expansion for a given order
  */
  static int getTaylorExpansionSize(const int a_Q);
};

#endif
