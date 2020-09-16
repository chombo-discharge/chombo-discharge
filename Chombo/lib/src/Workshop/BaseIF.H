#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEIF_H_
#define _BASEIF_H_

#include "RealVect.H"
#include "ProblemDomain.H"
#include "IndexTM.H"
#include "EB_TYPEDEFS.H"
#include "VolIndex.H"
#include "FaceIndex.H"
#include "Notation.H"
#include "GeometryService.H"
#include "IndexedMoments.H"
#include "NamespaceHeader.H"

///
/**
   This is the base class for an implicit function specification of geometry.
   All that is needed is a constructor/destructor, a method to give the value
   of the function at any point in space (in 2D or 3D), and a factory method
   (these will probably all be very similar).
 */
class BaseIF
{

public:

  /// Default constructor
  BaseIF()
  {
  }

  /// Default destructor
  virtual ~BaseIF()
  {
  }

  ///return int x^p dV for the vof
  virtual IndMomSpaceDim getExactVolumeMoments(const VolIndex    & a_vof,
                                               const Real        & a_dx) const
  {
    IndMomSpaceDim retval;
    retval.setToZero();
    MayDay::Error("not implemented");
    return retval;
  }
                                          
  
  ///return int_eb x^p dA for the eb
  virtual IndMomSpaceDim getExactEBMoments(const VolIndex   & a_vof,
                                           const Real       & a_dx) const
  {
    IndMomSpaceDim retval;
    retval.setToZero();
    MayDay::Error("not implemented");
    return retval;
  }
                                          
                                          
  ///return int_eb x^p n_i dA for the eb
  virtual IndMomSpaceDim getExactEBNormalMoments(const VolIndex   & a_vof,
                                                 const Real       & a_dx,
                                                 const int        & a_ni) const
  {
    IndMomSpaceDim retval;
    retval.setToZero();
    MayDay::Error("not implemented");
    return retval;
  }
                                          
  ///
  virtual IndMomSpaceDim getExactEBNormalPartialDerivs(const VolIndex   & a_vof,
                                                       const Real       & a_dx,
                                                       const int        & a_ni) const
  {
    IndMomSpaceDim retval;
    retval.setToZero();
    MayDay::Error("not implemented");
    return retval;
  }

  ///return int x^p  dA for the face
  virtual IndMomSDMinOne getExactFaceMoments(const FaceIndex  & a_face,
                                             const Real       & a_dx) const
  {
    IndMomSDMinOne retval;
    retval.setToZero();
    MayDay::Error("not implemented");
    return retval;
  }
                                          
  
  
  ///
  /**
   Return the value of the function at a_point.  When delineating a domain,
   the level set value=0 represents the boundary and value<0 is inside the
   fluid.
  */
  virtual Real value(const RealVect& a_point) const = 0;

  ///return the partial derivative at the point
  virtual Real derivative(const  IntVect& a_deriv,
                          const RealVect& a_point) const
  {
    MayDay::Error("not implemented");
    return 0.;
  }

  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    Real retval= 0;
    if (a_partialDerivative == IndexTM<int,GLOBALDIM>::Zero)
      {
        retval = value(a_point);
      }
    else
      {
        retval = derivative(a_partialDerivative, a_point);
      }
    return retval;
  }


  virtual Real value(const IntVect  & a_deriv,
                     const RealVect & a_point) const
  {
    IndexTM<int,GLOBALDIM>   deriv;
    IndexTM<Real,GLOBALDIM>  point;
    for(int idir = 0; idir < SpaceDim; idir ++)
      {
        deriv[idir] = a_deriv[idir];
        point[idir] = a_point[idir];
      }
    return value(deriv, point);
  }

  virtual bool fastIntersection(const Box&           a_region,
                                const ProblemDomain& a_domain,
                                const RealVect&      a_origin,
                                const Real&          a_dx) const
  {
    RealVect low, high;
    corners(a_region, a_origin, a_dx, low, high);
    return fastIntersection(low, high);
  }

  virtual bool fastIntersection(const RealVect& a_low,
                                const RealVect& a_high) const
  {return false;}

  virtual GeometryService::InOut InsideOutside(const Box&           a_region,
                                               const ProblemDomain& a_domain,
                                               const RealVect&      a_origin,
                                               const Real&          a_dx) const
  {
    RealVect low, high;
    corners(a_region, a_origin, a_dx, low, high);
    return InsideOutside(low, high);
  }

  virtual GeometryService::InOut InsideOutside(const RealVect& a_low,
                                               const RealVect& a_high) const

  {
    MayDay::Abort("This class has not implemented a fastIntersection operation");
    return GeometryService::Irregular;
  }

  ///
  /**
     Return the value of the function at a_point (of type INdexTM).
  */
  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    CH_assert(GLOBALDIM==SpaceDim);
    RealVect rvpoint;
    for(int idir = 0; idir < SpaceDim; idir++)
      {
        rvpoint[idir] = a_point[idir];
      }
    return value(rvpoint);
  };

  ///
  /**
     Return the derivative of the function at a_point (of type INdexTM).
  */
  virtual Real derivative(const IndexTM< int,GLOBALDIM>& a_deriv,
                          const IndexTM<Real,GLOBALDIM>& a_point
                          ) const
  {
    CH_assert(GLOBALDIM==SpaceDim);
    RealVect rvpoint;
    IntVect  ivderiv;
    for(int idir = 0; idir < SpaceDim; idir++)
      {
        rvpoint[idir] = a_point[idir];
        ivderiv[idir] = a_deriv[idir];
      }
    return derivative(ivderiv, rvpoint);
  };
  ///
  /**
     Return a newly allocated derived class.  The responsibility
     for deleting the memory is left to the calling function.
  */
  virtual BaseIF* newImplicitFunction() const = 0;

  virtual void print(ostream& out) const
  {
    MayDay::Abort("Print function not implemented");
  };

    ///
  /**
     An Implicit Function has three options for implementing this function
     1) do nothing, allow the empty base implementation to remain in place as a null-op
     2) take the makeGrids call as a directive: Here are the grids EBIndexSpace is wanting to use, configure
         yourself accordingly to make this efficient for you.  
     3) discard the DisjointBoxLayout EBIndexSpace would like and insert your own implementation of layout
        EBIndexSpace will faithfully use a_grids returned from this function, including it's load balance. 
  */
  virtual void makeGrids( const ProblemDomain&      a_domain,
                          DisjointBoxLayout&        a_grids,
                          const int&                a_maxGridSize,
                          const int&                a_maxIrregGridSize )
  {
    //default operation: Implicit Function and EBIndexSpace do not know about each other's data layout.
  }
  
  static void corners(const Box& a_region, const RealVect& a_origin, const Real& a_dx,
                      RealVect& a_lo, RealVect& a_hi)
  {
    a_lo = a_origin + RealVect(a_region.smallEnd())*a_dx;
    a_hi = a_origin + RealVect(a_region.bigEnd()+IntVect::Unit)*a_dx;
  }

  ///
  /**
     Inform the implicit function, IF, of a BoxLayout change in the
     application that is using it.  If the implicit function uses
     distributed data, this gives it a chance to react to the spatial
     layout of the application that is using it.

     This is a empty implementation for IFs that don't need this
     functionality.
  */
  virtual void boxLayoutChanged(const DisjointBoxLayout & a_newBoxLayout,
                                const RealVect          & a_dx)
  {
  }
};

#include "NamespaceFooter.H"

#endif
