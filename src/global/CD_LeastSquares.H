/*!
  @file   CD_LeastSquares.H
  @brief  Agglomeration of various routines that use least squares computations
  @author Robert Marskar
  @date   May 2016
*/

#ifndef _CD_LeastSquares_H_
#define _CD_LeastSquares_H_

#include "Stencils.H"
#include "EBISBox.H"
#include "IntVectSet.H"
#include "RealVect.H"
#include "ProblemDomain.H"
#include "VolIndex.H"
#include "PolyGeom.H"
#include "EBISLayout.H"

namespace ChomboDischarge {
  namespace LeastSquares {

    enum class CellPosition {
      Center,
      Centroid,
      Boundary
    };

    /*!
      @brief Compute the position (ignoring the "origin) of a Vof. 
      @param[in] a_position  Position within the grid cell (e.g. center, centroid, boundary).
      @param[in] a_vof       VolIndex to compute the position of
      @param[in] a_ebisbox   EBISBox
      @param[in] a_dx        Grid resolution
      @return Tosition of a Vof, ignoring the "origin" of the coordinate system since this is only used together with LeastSquares::displacement(...)
    */
    RealVect position(const CellPosition a_position,
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
    RealVect displacement(const CellPosition a_from,
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
    RealVect displacement(const CellPosition a_from,
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
    Vector<RealVect> getDisplacements(const CellPosition      a_from,
				      const CellPosition      a_to,
				      const VolIndex&         a_fromVof,
				      const Vector<VolIndex>& a_toVofs,
				      const EBISBox&          a_ebisbox,
				      const Real&             a_dx);

    /*!
      @brief Removes equations from LSQ system if displacement length is shorter than a_tolerance.
    */
    void removeEquations(Vector<VolIndex>& a_allVofs, Vector<RealVect>& a_displacements, const Real a_tolerance);

    /*!
      @brief Create a list of weights. This routine returns a list of diagonal weights for a least squares system. The weights
      are given as 1/||x1-x0||^a where a is the power. 
      @param[in] a_displacements Displacement vectors, i.e. (x1-x0,x2-x0, x3-x0,...,(xN-x0)^T. 
      @param[in] a_power Power factor for weights.
    */
    Vector<Real> makeDiagWeights(const Vector<RealVect>& a_displacements, const int a_pow);

    /*!
      @brief Get a least squares gradient stencil solution for a system of first-order Taylor extrapolation for computing the gradient at a point.
      @param[in] a_allVofs       Vofs for which the displacement vectors were computed. 
      @param[in] a_displacements Displacement vectors from Vofs to the interpolation point
      @param[in] a_weights       Weights for each system of equations. 
      @note This computes weights using makeDiagWeights(a_displacements, a_p) and then calls the other version. 
      @return Stencil for obtaining the gradient in a point. Stencils for each direction are put in the VoFStencil.variable(i). 
    */
    VoFStencil computeGradStenOrderOne(const Vector<VolIndex>& a_allVofs,
				       const Vector<RealVect>& a_displacements,
				       const int&              a_p);

    /*!
      @brief Get a least squares gradient stencil solution for a system of first-order Taylor extrapolation for computing the gradient at a point.
      @param[in] a_allVofs       Vofs for which the displacement vectors were computed. 
      @param[in] a_displacements Displacement vectors from Vofs to the interpolation point
      @param[in] a_weights       Weights for each system of equations. 
      @return Stencil for obtaining the gradient in a point. Stencils for each direction are put in the VoFStencil.variable(i). 
      @note This routine will cause a run-time error if you have less than SpaceDim equations, which is the minimum number of equations needed for obtaining a gradient. 
      However, the routine does *not* check if the inputs are degenerate. 
    */
    VoFStencil computeGradStenOrderOne(const Vector<VolIndex>& a_allVofs,
				       const Vector<RealVect>& a_displacements,
				       const Vector<Real>&     a_weights);
  
    /*!
      @brief Compute a least squares gradient stencil on the EB centroid
      @param[in]  a_vof Input vof
      @param[in]  a_ebisbox ebisbox
      @param[in]  a_dx Resolution
      @param[in]  a_p Scaling
      @return     Returns a least squares stencil for evaluating the a first order approximation to the normal derivative on the boundary centroid. 
      @details    This will first try to get all Vofs defined in the quadrant defined by the normal vector, i.e. 3 cells in 2D and up to 7 cells in 3D. If that
      does not work the routine will try to fetch a symmetric stencil. If that does not work it will get all Vofs within a certain radius. 
    */
    VoFStencil getBndryGradStenOrderOne(const VolIndex& a_vof,
					const EBISBox&  a_ebisbox,
					const Real&     a_dx,
					const int       a_p);

    /*!
      @brief Assuming that a_stencil is a gradient stencil, project it along a direction. 
      @param[in] a_stencil Stencil describing a gradient. 
      @param[in] a_projection Vector to project stencil along. 
      @return  Returns gradient stencil projected along a_project. 
      @note Users need to make sure ||a_projection|| = 1 for this routine to make sense. 
    */
    VoFStencil projectGradSten(const VoFStencil& a_stencil, const RealVect& a_projection);

    /*!
      @brief Return the sum of weights in the stencil, for a specific variable. 
      @param[in] a_stencil Stencil
      @param[in] a_variable Which variable to sum. 
      @return    Returns sum of weights for a specific variable in the stencil. 
    */
    Real sumWeights(const VoFStencil& a_stencil, const int a_variable);

    /*!
      @brief Compute the sum of all weights in a stencil. 
      @param[in] a_stencil Stencil
      @return    Returns sum of weights in a_stencil (over all variables).
    */
    Real sumAllWeights(const VoFStencil& a_stencil);
  
    /*!
      @brief Get the size of a Taylor expansion for a given order
      @param[in] a_order Desired order in Taylor expansion
      @return Returns number of terms in Taylor expansion. 
    */
    int getTaylorExpansionSize(const int a_order);

    /*!
      @brief Compute an interpolation stencil to specified order by solving a least squares system
      @param[in] a_allVofs Vofs to include in the interpolation method. 
      @param[in] a_displacements Displacement vectors from the Vofs to the interpolation point
      @param[in] a_order         Order of the interpolation. 
      @param[in] a_weightP       Weighting order for the least squares system. Must be > 0 to have an effect. 
      @return Stencil for the least squares solution to the interpolation problem. 
      @note This computes weights using makeDiagWeights(a_displacements, a_pow) and then calls the other version. 
    */
    VoFStencil computeInterpolationStencil(const Vector<VolIndex>& a_allVofs,
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
    VoFStencil computeInterpolationStencil(const Vector<VolIndex>& a_allVofs,
					   const Vector<RealVect>& a_displacements,
					   const Vector<Real>&     a_weights,
					   const int               a_order);

    /*!
      @brief Compute a least squares interpolation to a specified order. 
      @details This is the general version which lets the caller specify why derivatives he wants out of the interpolation. The user
      specifies this via a_derivs where each IntVect represent a differentiation. E.g. IntVect(0,0,0) is just f, IntVect(1,0,0) = d/dx,
      IntVect(1,1,2) = (d^4)/(dx dy dz^2) and so on. The stencils for the directional derivatives are coded onto the return map. E.g.
      if a_derivs contains IntVect(1,0,0) the stencil for d/dx is found in map.at(IntVect(1,0,0)). 
      @param[in] a_derives       Specification of which derivates. 
      @param[in] a_allVofs       Vofs to include in the interpolation method. 
      @param[in] a_displacements Displacement vectors from the Vofs to the interpolation point
      @param[in] a_p             Weighting order for the least squares system. Must be > 0 to have an effect. 
      @param[in] a_order         Order of the interpolation. 
      @return Stencils for each of the derivatives specified in the input argument. 
      @note This calls the more general version. 
    */
    std::map<IntVect, VoFStencil> computeInterpolationStencil(const IntVectSet&       a_derivs,
							      const Vector<VolIndex>& a_allVofs,
							      const Vector<RealVect>& a_displacements,
							      const int               a_p,
							      const int               a_order);

    /*!
      @brief Compute a least squares interpolation to a specified order. 
      @details This is the general version which lets the caller specify why derivatives he wants out of the interpolation. The user
      specifies this via a_derivs where each IntVect represent a differentiation. E.g. IntVect(0,0,0) is just f, IntVect(1,0,0) = d/dx,
      IntVect(1,1,2) = (d^4)/(dx dy dz^2) and so on. The stencils for the directional derivatives are coded onto the return map. E.g.
      if a_derivs contains IntVect(1,0,0) the stencil for d/dx is found in map.at(IntVect(1,0,0)). 
      @param[in] a_derives       Specification of which derivates. 
      @param[in] a_allVofs       Vofs to include in the interpolation method. 
      @param[in] a_displacements Displacement vectors from the Vofs to the interpolation point
      @param[in] a_weights       Weights for the least squares system. 
      @param[in] a_order         Order of the returned interpolation. 
      @return Stencils for each of the derivatives specified in the input argument. 
      @note This will throw an error if you don't have enough equations for obtaining the specified order. 
      @note The return map is always initialized with empty stencils. If computing the pseudoinverse fails, the returned stencils are all empty. 
    */
    std::map<IntVect, VoFStencil> computeInterpolationStencil(const IntVectSet&       a_derivs,
							      const Vector<VolIndex>& a_allVofs,
							      const Vector<RealVect>& a_displacements,
							      const Vector<Real>&     a_weights,
							      const int               a_order);

    /*!
      @brief Get an interpolation stencil to a point in a cell with a specified order and radius. 
    */
    VoFStencil getInterpolationStencilUsingAllVofsInRadius(const CellPosition a_cellPos,
							   const CellPosition a_otherCellsPos,
							   const VolIndex&    a_startVof,
							   const EBISBox&     a_ebisbox,
							   const Real         a_dx,
							   const int          a_p,
							   const int          a_radius,
							   const int          a_order);
  

  };
};

#endif
