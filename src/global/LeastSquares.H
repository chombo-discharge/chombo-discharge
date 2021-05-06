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

  enum class CellPosition {
    Center,
    Centroid,
    Boundary
  };

  /*!
    @brief Get All VoFs in quadrant
  */
  Vector<VolIndex> getAllVoFsInQuadrant(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const RealVect a_normal);
  
  /*!
    @brief Get all VoFs within radius a_radius. Add the start vof if you want it. 
  */
  Vector<VolIndex> getAllVoFsInRadius(const VolIndex& a_startVoF, const EBISBox& a_ebisbox, const int a_radius, const bool a_addStartVoF);

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
  RealVect displacement(const CellPosition a_from,
			const CellPosition a_to,
			const VolIndex&    a_fromVoF,
			const VolIndex&    a_toVoF,
			const EBISBox&     a_ebisbox,
			const Real&        a_dx);

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
    @brief Removes equations from LSQ system if displacement length is shorter than a_tolerance.
  */
  void removeEquations(Vector<VolIndex>& a_allVoFs, Vector<RealVect>& a_displacements, const Real a_tolerance);

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
				 const Vector<Real>&     a_weights);
  
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
    @brief Get the size of a Taylor expansion for a given order
  */
  int getTaylorExpansionSize(const int a_Q);

  /*!
    @brief Generate a sequence of multi-indices |q| <= Q in lexicographical order
    @param[in] a_Q Largest multi-index
  */
  Vector<MultiIndex> getMultiIndicesLexiOrder(const int a_Q);
};

#endif
