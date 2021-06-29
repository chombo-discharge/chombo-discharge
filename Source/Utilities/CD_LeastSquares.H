/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LeastSquares.H
  @brief  Agglomeration of various routines that use least squares computations
  @author Robert Marskar
*/

#ifndef CD_LeastSquares_H
#define CD_LeastSquares_H

// Chombo includes
#include <Stencils.H>
#include <EBISBox.H>
#include <IntVectSet.H>
#include <RealVect.H>
#include <ProblemDomain.H>
#include <VolIndex.H>
#include <PolyGeom.H>
#include <EBISLayout.H>

// Our includes
#include <CD_VofUtils.H>
#include <CD_NamespaceHeader.H>

class LeastSquares {
public:

  /*!
    @brief Static class, no construction here.
  */
  LeastSquares() = delete;

  /*!
    @brief Static class, no destruction here. 
  */
  ~LeastSquares() = delete;
  
  enum class CellPosition {
    Center,
    Centroid,
    Boundary
  };

  using Connectivity = VofUtils::Connectivity;
  using Neighborhood = VofUtils::Neighborhood;

  /*!
    @brief Get an interpolation stencil to a point in a cell with a specified order and radius. 
    @param[in] a_cellPos        Position in the cell of the starting Vof.
    @param[in] a_otherCellsPos  Positions in the other cells.
    @param[in] a_startVof       The vof to compute for.
    @param[in] a_ebisbox        EBISBox.
    @param[in] a_dx             Grid resolution.
    @param[in] a_p              Power factor for weights.
    @param[in] a_radius         Stencil radius.
    @param[in] a_order          Interpolation order. 
    @param[in] a_addStartingVof Add starting vof or not (can't be used with a_p > 0) as it would lead to ill-conditioned systems. 
  */
  static VoFStencil getInterpolationStencil(const CellPosition a_cellPos,
					    const CellPosition a_otherCellsPos,
					    const Connectivity a_connectivity,					    
					    const VolIndex&    a_startVof,
					    const EBISBox&     a_ebisbox,
					    const Real         a_dx,
					    const int          a_p,
					    const int          a_radius,
					    const int          a_order,
					    const bool         a_addStartingVof);

  /*!
    @brief Compute a least squares gradient stencil on the EB centroid with specified order. 
    @details This routine eliminates the value on the boundary and the user will have to fetch the weight for this point from the returned stencil. 
    @param[in]  a_vof Input vof
    @param[in]  a_ebisbox ebisbox
    @param[in]  a_dx Resolution
    @param[in]  a_p Scaling
    @param[in]  a_order
    @return     Returns a least squares stencil for evaluating the a first order approximation to the normal derivative on the boundary centroid. 
    @details    This will first try to get all Vofs defined in the quadrant defined by the normal vector, i.e. 3 cells in 2D and up to 7 cells in 3D. If that
    does not work the routine will try to fetch a symmetric stencil. If that does not work it will get all Vofs within a certain radius. 
  */
  static VoFStencil getBndryGradSten(const VolIndex&    a_vof,
				     const Neighborhood a_neighborhood,
				     const CellPosition a_cellPositions,
				     const EBISBox&     a_ebisbox,
				     const Real         a_dx,
				     const int          a_radius,
				     const int          a_p,
				     const int          a_order);

  /*!
    @brief Return the sum of weights in the stencil, for a specific variable. 
    @param[in] a_stencil Stencil
    @param[in] a_variable Which variable to sum. 
    @return    Returns sum of weights for a specific variable in the stencil. 
  */
  static Real sumWeights(const VoFStencil& a_stencil, const int a_variable);

  /*!
    @brief Compute the sum of all weights in a stencil. 
    @param[in] a_stencil Stencil
    @return    Returns sum of weights in a_stencil (over all variables).
  */
  static Real sumAllWeights(const VoFStencil& a_stencil);

  /*!
    @brief Assuming that a_stencil is a gradient stencil, project it along a direction. 
    @param[in] a_stencil Stencil describing a gradient. 
    @param[in] a_projection Vector to project stencil along. 
    @return  Returns gradient stencil projected along a_project. 
    @note Users need to make sure ||a_projection|| = 1 for this routine to make sense. 
  */
  static VoFStencil projectGradSten(const VoFStencil& a_stencil, const RealVect& a_projection);
  
  /*!
    @brief Compute the position (ignoring the "origin) of a Vof. 
    @param[in] a_position  Position within the grid cell (e.g. center, centroid, boundary).
    @param[in] a_vof       VolIndex to compute the position of
    @param[in] a_ebisbox   EBISBox
    @param[in] a_dx        Grid resolution
    @return Tosition of a Vof, ignoring the "origin" of the coordinate system since this is only used together with LeastSquares::displacement(...)
  */
  static RealVect position(const CellPosition a_position,
			   const VolIndex&    a_vof,
			   const EBISBox&     a_ebisbox,
			   const Real&        a_dx);
  
  /*!
    @brief Computes the distance between two Vofs that are defined on the same grid level. 
    @param[in] a_from     Specifying position in the cell defined by a_fromVof
    @param[in] a_to       Specifying position in the cell defined by a_toVof
    @param[in] a_fromVof  Vof to compute the distance from
    @param[in] a_toVof    Vof to compute the distance to
    @param[in] a_ebisbox  EBISBox
    @param[in] a_dx       Grid resolution
    @returns The vector pos(a_toVof) - pos(a_fromVof). 
  */
  static RealVect displacement(const CellPosition a_from,
			       const CellPosition a_to,
			       const VolIndex&    a_fromVof,
			       const VolIndex&    a_toVof,
			       const EBISBox&     a_ebisbox,
			       const Real&        a_dx);

  /*!
    @brief Computes the distance between two Vofs that are defined on two different grid levels. 
    @param[in] a_from        Specifying position in the cell defined by a_fromVof
    @param[in] a_to          Specifying position in the cell defined by a_toVof
    @param[in] a_fromVof     Vof to compute the distance from
    @param[in] a_toVof       Vof to compute the distance to
    @param[in] a_ebisboxFrom EBISBox which can reach a_fromVof
    @param[in] a_ebisboxTo   EBISBox which can reach a_toVof
    @param[in] a_dxFrom      Grid resolution for the a_fromVof vof
    @param[in] a_dxTo        Grid resolution for the a_toVof vof
    @returns The vector pos(a_toVof) - pos(a_fromVof). 
  */
  static RealVect displacement(const CellPosition a_from,
			       const CellPosition a_to,
			       const VolIndex&    a_fromVof,
			       const VolIndex&    a_toVof,
			       const EBISBox&     a_ebisboxFrom,
			       const EBISBox&     a_ebisboxTo,
			       const Real&        a_dxFrom,
			       const Real&        a_dxTo);

  /*!
    @brief Get displacement vectors from a cell position to another cell position for a list of Vofs, i.e. x = (xFrom - xTo).
    @param[in] a_from    From this position
    @param[in] a_to      To this position. 
    @param[in] a_curVof  Origin vof
    @param[in] a_toVofs Vofs to compute the distance to/from
    @param[in] a_ebisbox EBISBox
    @param[in] a_dx      Grid resolution
  */
  static Vector<RealVect> getDisplacements(const CellPosition      a_from,
					   const CellPosition      a_to,
					   const VolIndex&         a_fromVof,
					   const Vector<VolIndex>& a_toVofs,
					   const EBISBox&          a_ebisbox,
					   const Real&             a_dx);

  /*!
    @brief Create a list of weights. This routine returns a list of diagonal weights for a least squares system. The weights
    are given as 1/||x1-x0||^a where a is the power. 
    @param[in] a_displacements Displacement vectors, i.e. (x1-x0,x2-x0, x3-x0,...,(xN-x0)^T. 
    @param[in] a_power Power factor for weights.
  */
  static Vector<Real> makeDiagWeights(const Vector<RealVect>& a_displacements, const int a_pow);

  /*!
    @brief Get a least squares gradient stencil solution for computing the gradient at a point. This routine eliminates the value in the point where
    the gradient is computed (which is assumed to be known.)
    @param[in] a_allVofs       Vofs for which the displacement vectors were computed. 
    @param[in] a_displacements Displacement vectors from Vofs to the interpolation point
    @param[in] a_p             Weighting exponent for system
    @param[in] a_order         Desired order. 
    @note This computes weights using makeDiagWeights(a_displacements, a_p) and then calls the other version. 
    @return Stencil for obtaining the gradient in a point. Stencils for each direction are put in the VoFStencil.variable(i). 
  */
  static VoFStencil computeGradSten(const Vector<VolIndex>& a_allVofs,
				    const Vector<RealVect>& a_displacements,
				    const int               a_p,
				    const int               a_order);

  /*!
    @brief Get a least squares gradient stencil solution for a system of first-order Taylor extrapolation for computing the gradient at a point.
    @details This routine eliminates the value in the point where the gradient is computed (which is assumed to be known.)
    @param[in] a_allVofs       Vofs for which the displacement vectors were computed. 
    @param[in] a_displacements Displacement vectors from Vofs to the interpolation point
    @param[in] a_weights       Weights for each system of equations. 
    @param[in] a_order         Interpolation order
    @return Stencil for obtaining the gradient in a point. Stencils for each direction are put in the VoFStencil.variable(i). 
    @note This routine will cause a run-time error if you have less than SpaceDim equations, which is the minimum number of equations needed for obtaining a gradient. 
    However, the routine does *not* check if the inputs are degenerate. 
  */
  static VoFStencil computeGradSten(const Vector<VolIndex>& a_allVofs,
				    const Vector<RealVect>& a_displacements,
				    const Vector<Real>&     a_weights,
				    const int               a_order);
  
  /*!
    @brief Get the size of a Taylor expansion for a given order
    @param[in] a_order Desired order in Taylor expansion
    @return Returns number of terms in Taylor expansion. 
  */
  static int getTaylorExpansionSize(const int a_order);

  /*!
    @brief Compute an interpolation stencil to specified order by solving a least squares system
    @param[in] a_allVofs Vofs to include in the interpolation method. 
    @param[in] a_displacements Displacement vectors from the Vofs to the interpolation point
    @param[in] a_order         Order of the interpolation. 
    @param[in] a_weightP       Weighting order for the least squares system. Must be > 0 to have an effect. 
    @return Stencil for the least squares solution to the interpolation problem. 
    @note This computes weights using makeDiagWeights(a_displacements, a_pow) and then calls the other version. 
  */
  static VoFStencil computeInterpolationStencil(const Vector<VolIndex>& a_allVofs,
						const Vector<RealVect>& a_displacements,
						const int               a_pow,
						const int               a_order);
  
  /*!
    @brief Compute an interpolation stencil to specified order by solving a least squares system
    @param[in] a_allVofs Vofs to include in the interpolation method. 
    @param[in] a_displacements Displacement vectors from the Vofs to the interpolation point
    @param[in] a_order         Order of the interpolation. 
    @param[in] a_weightP       Weighting order for the least squares system. Must be > 0 to have an effect. 
    @return Stencil for the least squares solution to the interpolation problem. 
    @note This calls the more general version with a_derivs = IntVectSet(IntVect(0,0)). 
    @note This will throw an error if you don't have enough equations for obtaining the specified order. 
  */
  static VoFStencil computeInterpolationStencil(const Vector<VolIndex>& a_allVofs,
						const Vector<RealVect>& a_displacements,
						const Vector<Real>&     a_weights,
						const int               a_order);

  /*!
    @brief Compute a least squares interpolation to a specified order. 
    @details This is the general version which lets the caller specify why derivatives he wants out of the interpolation. The user
    specifies this via a_derivs where each IntVect represent a differentiation. E.g. IntVect(0,0,0) is just f, IntVect(1,0,0) = d/dx,
    IntVect(1,1,2) = (d^4)/(dx dy dz^2) and so on. The stencils for the directional derivatives are coded onto the return map. E.g.
    if a_derivs contains IntVect(1,0,0) the stencil for d/dx is found in map.at(IntVect(1,0,0)). The user can solve a smaller system
    if some of the terms in the expansion are known, typically used when computing an approximation to the gradient at a point in space
    where the scalar value is known. Note that this will modify the right hand side of the system, and the user will have to figure out 
    how to make sense of the stencil. 
    @param[in] a_derives       Specification of which terms to obtain from the series.
    @param[in] a_knownTerms    Which terms in the Taylor series are known, and will be eliminated from the least squares system. 
    @param[in] a_allVofs       Vofs to include in the interpolation method. 
    @param[in] a_displacements Displacement vectors from the Vofs to the interpolation point
    @param[in] a_p             Weighting order for the least squares system. Must be > 0 to have an effect. 
    @param[in] a_order         Order of the interpolation. 
    @return Stencils for each of the derivatives specified in the input argument. 
    @note This calls the more general version. 
  */
  static std::map<IntVect, VoFStencil> computeInterpolationStencil(const IntVectSet&       a_derivs,
								   const IntVectSet&       a_knownTerms,
								   const Vector<VolIndex>& a_allVofs,
								   const Vector<RealVect>& a_displacements,
								   const int               a_p,
								   const int               a_order);

  /*!
    @brief Compute a least squares interpolation to a specified order. 
    @details This is the general version which lets the caller specify why derivatives he wants out of the interpolation. The user
    specifies this via a_derivs where each IntVect represent a differentiation. E.g. IntVect(0,0,0) is just f, IntVect(1,0,0) = d/dx,
    IntVect(1,1,2) = (d^4)/(dx dy dz^2) and so on. The stencils for the directional derivatives are coded onto the return map. E.g.
    if a_derivs contains IntVect(1,0,0) the stencil for d/dx is found in map.at(IntVect(1,0,0)). The user can solve a smaller system
    if some of the terms in the expansion are known, typically used when computing an approximation to the gradient at a point in space
    where the scalar value is known. Note that this will modify the right hand side of the system, and the user will have to figure out 
    how to make sense of the stencil. 
    @param[in] a_derives       Specification of which unknowns in the Taylor series will be returned.
    @param[in] a_knownTerms    Which terms in the Taylor series are known. 
    @param[in] a_allVofs       Vofs to include in the interpolation method. 
    @param[in] a_displacements Displacement vectors from the Vofs to the interpolation point
    @param[in] a_weights       Weights for the least squares system. 
    @param[in] a_order         Order of the returned interpolation. 
    @return Stencils for each of the derivatives specified in the input argument. 
    @note This will throw an error if you don't have enough equations for obtaining the specified order. 
    @note The return map is always initialized with empty stencils. If computing the pseudoinverse fails, the returned stencils are all empty. 
  */
  static std::map<IntVect, VoFStencil> computeInterpolationStencil(const IntVectSet&       a_derivs,
								   const IntVectSet&       a_knownTerms,
								   const Vector<VolIndex>& a_allVofs,
								   const Vector<RealVect>& a_displacements,
								   const Vector<Real>&     a_weights,
								   const int               a_order);



  /*!
    @brief Trim a system of equations down to a specified size. 
    @details This will discard displacement vectors with the longest distance. If the system size is already equal to or smaller then this routine does nothing. 
    @param[inout] a_vofs          Grid cells
    @param[inout] a_displacements Displacement vectors
    @param[in]    a_systemSize    Desired system size. 
  */
  static void trimSystem(Vector<VolIndex>& a_vofs,
			 Vector<RealVect>& a_displacements,
			 const int         a_systemSize);

};

#include <CD_NamespaceFooter.H>

#endif
