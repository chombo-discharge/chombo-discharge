#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PYIF_H_
#define _PYIF_H_

#ifdef CH_USE_EB

#include "BaseIF.H"
#include <Python.h>

//! \class PyIF
//! This class allows one to construct implicit functions using constructive
//! solid geometry. Set operations and primitives are supported.
class PyIF: public BaseIF
{

  public:

  //! Returns true if the given Python object can be interpreted as an
  //! implicit function, false otherwise.
  //! \param a_pyObject The Python object in question.
  static bool isValid(PyObject* a_pyObject);

  //! Construct an implicit function from the given Python object.
  //! \param a_pyObject The Python object used to construct the implicit function.
  explicit PyIF(PyObject* a_pyObject);

  //! Destructor.
  ~PyIF();

  //! Overridden methods.

  Real value(const RealVect& a_point) const;

  Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
             const IndexTM<Real,GLOBALDIM>& a_point) const;

  bool fastIntersection(const Box&           a_region,
                        const ProblemDomain& a_domain,
                        const RealVect&      a_origin,
                        const Real&          a_dx) const;

  bool fastIntersection(const RealVect& a_low,
                        const RealVect& a_high) const;

  GeometryService::InOut InsideOutside(const Box&           a_region,
                                       const ProblemDomain& a_domain,
                                       const RealVect&      a_origin,
                                       const Real&          a_dx) const;

  GeometryService::InOut InsideOutside(const RealVect& a_low,
                                       const RealVect& a_high) const;

  Real value(const IndexTM<Real,GLOBALDIM>& a_point) const;

  BaseIF* newImplicitFunction() const;

  void print(ostream& out) const;

  void boxLayoutChanged(const DisjointBoxLayout & a_newBoxLayout,
                        const RealVect          & a_dx);

  private:

  // The underlying C++ implicit function.
  RefCountedPtr<BaseIF> m_IF;
};

#ifndef PyIF_Check
#define PyIF_Check(x) PyIF::isValid(x)
#endif

#endif

#endif
