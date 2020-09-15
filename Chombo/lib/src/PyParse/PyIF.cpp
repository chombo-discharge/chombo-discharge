#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef CH_USE_EB

#include "PyIF.H"
#include "PyWorkshop.H"

//-----------------------------------------------------------------------
bool
PyIF::
isValid(PyObject* a_pyObject)
{
  return PyImplicitFunction_Check(a_pyObject);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyIF::
PyIF(PyObject* a_pyObject):
  m_IF()
{
  CH_assert(PyIF::isValid(a_pyObject));
  m_IF = PyImplicitFunction_IF(a_pyObject);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyIF::
~PyIF()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
PyIF::
value(const RealVect& a_point) const
{
  return m_IF->value(a_point);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
PyIF::
value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
      const IndexTM<Real,GLOBALDIM>& a_point) const
{
  return m_IF->value(a_partialDerivative, a_point);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool
PyIF::
fastIntersection(const Box&           a_region,
                 const ProblemDomain& a_domain,
                 const RealVect&      a_origin,
                 const Real&          a_dx) const
{
  return m_IF->fastIntersection(a_region, a_domain, a_origin, a_dx);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool
PyIF::
fastIntersection(const RealVect& a_low,
                 const RealVect& a_high) const
{
  return m_IF->fastIntersection(a_low, a_high);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
GeometryService::InOut
PyIF::
InsideOutside(const Box&           a_region,
              const ProblemDomain& a_domain,
              const RealVect&      a_origin,
              const Real&          a_dx) const
{
  return m_IF->InsideOutside(a_region, a_domain, a_origin, a_dx);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
GeometryService::InOut
PyIF::
InsideOutside(const RealVect& a_low,
              const RealVect& a_high) const
{
  return m_IF->InsideOutside(a_low, a_high);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
PyIF::
value(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  return m_IF->value(a_point);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
BaseIF*
PyIF::
newImplicitFunction() const
{
  // If someone asks for a new implicit function, we just hand them one
  // created by the one that we're pointed at.
  return m_IF->newImplicitFunction();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PyIF::
print(ostream& out) const
{
  return m_IF->print(out);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PyIF::
boxLayoutChanged(const DisjointBoxLayout & a_newBoxLayout,
                 const RealVect          & a_dx)
{
  return m_IF->boxLayoutChanged(a_newBoxLayout, a_dx);
}
//-----------------------------------------------------------------------

#endif
