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

  enum class CellPosition {
    Center,
    Centroid,
    Boundary
  };

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
    @brief Returns position of a VoF. This ignores the "origin" of the coordinate system 
    since this is only used together with LeastSquares::displacement(...)
  */
  RealVect position(const CellPosition a_position,
		    const VolIndex&    a_vof,
		    const EBISBox&     a_ebisbox,
		    const Real&        a_dx);
  
  /*!
    @brief Distance between two VoFs
  */
  RealVect displacement(const CellPosition      a_from,
			const CellPosition      a_to,
			const VolIndex&         a_fromVoF,
			const VolIndex&         a_toVoF,
			const EBISBox&          a_ebisbox,
			const Real&             a_dx);

  /*!
    @brief Get displacement vectors from a cell position to another cell position for a list of VoFs, i.e. x = (xFrom - xTo).
    @param[in] a_from    From this position
    @param[in] a_to      To this position. 
    @param[in] a_curVoF  Origin vof
    @param[in] a_toVoFs VoFs to compute the distance to/from
    @param[in] a_ebisbox EBISBox
    @param[in] a_dx      Grid resolution
  */
  Vector<RealVect> getDisplacements(const CellPosition      a_from,
				    const CellPosition      a_to,
				    const VolIndex&         a_fromVoF,
				    const Vector<VolIndex>& a_toVoFs,
				    const EBISBox&          a_ebisbox,
				    const Real&             a_dx);

  /*!
    @brief Create a list of weights. This routine returns a list of diagonal weights for a least squares system. The weights
    are given as 1/||x1-x0||^a where a is the power. 
    @param[in] a_displacements Displacement vectors, i.e. (x1-x0,x2-x0, x3-x0,...,(xN-x0)^T. 
    @param[in] a_power Power factor for weights.
  */
  Vector<Real> makeDiagWeights(const Vector<RealVect>& a_displacements, const int a_pow);
				    

  /*!
    @brief Get a least squares gradient stencil by computing the Moore-Penrose pseudoinverse.
  */
  VoFStencil getGradStenOrderOne(const Vector<VolIndex>& a_allVoFs,
				 const Vector<RealVect>& a_displacements,
				 const int&              a_p);

  /*!
    @brief Compute a least squares gradient stencil on the EB centroid
    @param[in]  a_vof Input vof
    @param[in]  a_ebisbox ebisbox
    @param[in]  a_dx Resolution
    @param[in]  a_p Scaling
  */
  VoFStencil getBndryGradStenOrderOne(const VolIndex& a_vof,
				      const EBISBox&  a_ebisbox,
				      const Real&     a_dx,
				      const int       a_p);
  /*!
    @brief Get the least squares gradient stencil
    @param[out] a_stencil Gradient stencil (has SpaceDim components)
    @param[in]  a_monoVoFs VolIndexes
    @param[in]  a_displ Displacements
    @param[in]  a_Q Interpolation order. 
    @param[in]  a_pow Weighting "power", i.e. w = 1/|x1-x0|^p where p is power. 
  */
  bool getGradSten(VoFStencil&             a_stencil,
		   const Vector<VolIndex>& a_monoVoFs,
		   const Vector<RealVect>& a_displ,
		   const int               a_Q,
		   const int               a_weight);

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
  bool getBndryGradStencil(VoFStencil&     a_stencil,
			   const VolIndex& a_vof,
			   const EBISBox&  a_ebisbox,
			   const Real&     a_dx,
			   const int       a_order = 1,
			   const int       a_radius = 1,
			   const int       a_weight = 1);

  /*!
    @brief Assuming that a_stencil is a gradient stencil, project it along a direction. 
    @param[in] a_stencil Stencil
    @param[in] a_projection Vector to project stencil along. 
  */
  VoFStencil projectGradSten(const VoFStencil& a_stencil, const RealVect& a_projection);

  /*!
    @brief Return the sum of weights in the stencil, for a specific variable. 
    @param[in] a_stencil Stencil
    @param[in] a_variable Which variable to sum. 
  */
  Real sumWeights(const VoFStencil& a_stencil, const int a_variable = 0);

  /*!
    @brief Return the sum of weights in the stencil. 
    @param[in] a_stencil Stencil
  */
  Real sumAllWeights(const VoFStencil& a_stencil);
  
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
