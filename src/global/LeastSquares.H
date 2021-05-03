/*!
  @file   LeastSquares.H
  @brief  Agglomeration of various routines that use least squares computations
  @author Robert Marskar
  @date   May 2016
*/

#ifndef _LEASTSQUARES_H_
#define _LEASTSQUARES_H_

#include "MultiIndex.H"

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
namespace LeastSquares {

  /*!
    @brief Tolerance for stencil computations. 
  */
  static Real s_tol = 1.E-8;

  /*!
    @brief Create a list of weights. This routine returns a list of diagonal weights for a least squares system. The weights
    are given as 1/||x1-x0||^a where a is the power. 
    @param[in] a_displacements Displacement vectors, i.e. (x1-x0,x2-x0, x3-x0,...,(xN-x0)^T. 
    @param[in] a_power Power factor for weights.
  */
  Vector<Real> makeDiagWeights(const Vector<RealVect>& a_displacements, const int a_pow);

  /*!
    @brief Build the vector of distances that goes into LeastSquares::getCenterToCentroidInterpStencil
    @param[in]  a_vof Input vof
    @param[in]  a_monoVoFs VoFs that can be reached. This should include a_vof
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_dx Resolution
  */
  Vector<RealVect> getCenterToCentroidDisplacements(const VolIndex&         a_curVoF,
						    const Vector<VolIndex>& a_allVoFs,
						    const EBISBox&          a_ebisbox,
						    const Real&             a_dx);

  /*!
    @brief Build the vector of distances that goes into LeastSquares::getCenterToEBInterpStencil
    @param[in]  a_vof Input vof
    @param[in]  a_monoVoFs VoFs that can be reached. This should include a_vof
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_dx Resolution
  */
  Vector<RealVect> getCenterToEBDisplacements(const VolIndex&         a_vof,
					      const Vector<VolIndex>& a_monoVoFs,
					      const EBISBox&          a_ebisbox,
					      const Real&             a_dx);

  /*!
    @brief Build the vector of distances that goes into LeastSquares::getCentroidToCenterInterpStencil
    @param[in]  a_vof Inputs vof
    @param[in]  a_monoVoFs VoFs that can be reached. This should include a_vof
    @param[in]  a_ebisbox EBISBox
    @param[in]  a_dx Resolution
  */
  Vector<RealVect> getCentroidToCenterDisplacements(const VolIndex&         a_vof,
						    const Vector<VolIndex>& a_monoVoFs,
						    const EBISBox&          a_ebisbox,
						    const Real&             a_dx);

  /*!
    @brief Get the least squares gradient stencil
    @param[out] a_stencil Gradient stencil (has SpaceDim components)
    @param[in]  a_monoVoFs VolIndexes
    @param[in]  a_displ Displacements
    @param[in]  a_Q Interpolation order. 
    @param[in]  a_pow Weighting "power", i.e. w = 1/|x1-x0|^p where p is power. 
  */
  bool getLSqGradientStencil(VoFStencil&             a_stencil,
			     const Vector<VolIndex>& a_monoVoFs,
			     const Vector<RealVect>& a_displ,
			     const int               a_Q,
			     const int               a_pow);

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
  bool getEbCentroidGradientStencil(VoFStencil&     a_stencil,
				    const VolIndex& a_vof,
				    const EBISBox&  a_ebisbox,
				    const Real&     a_dx,
				    const int       a_order = 1,
				    const int       a_radius = 1,
				    const int       a_weight = 1);
  
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
  bool getCenterToCentroidInterpStencil(VoFStencil&     a_stencil,
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
  bool getCenterToEBInterpStencil(VoFStencil&     a_stencil,
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
  bool getCentroidToCenterInterpStencil(VoFStencil&     a_stencil,
					const VolIndex& a_vof,
					const EBISBox&  a_ebisbox,
					const Real&     a_dx,
					const int       a_order = 1,
					const int       a_radius = 1);


  /*!
    @brief Get the least squares interpolation stencil
    @param[out] a_stencil Stencil
    @param[in]  a_monoVoFs VolIndexes
    @param[in]  a_displ Displacements
    @param[in]  a_Q Interpolation order. 
  */
  bool getLSqInterpStencil(VoFStencil&             a_stencil,
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
  bool getVoFsRadius(Vector<VolIndex>& a_allVoFsRad,
		     const VolIndex&   a_vof,
		     const EBISBox&    a_ebisbox,
		     const int         a_minVoFs,
		     const int         a_order, 
		     const int         a_radius = -1);


  
  /*!
    @brief Generate a sequence of multi-indices |q| <= Q in lexicographical order
    @param[in] a_Q Largest multi-index
  */
  Vector<MultiIndex> getMultiIndicesLexiOrder(const int a_Q);

  /*!
    @brief Get the size of a Taylor expansion for a given order
  */
  int getTaylorExpansionSize(const int a_Q);
};

#include "LeastSquaresI.H"

#endif
