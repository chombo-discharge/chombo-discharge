#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTISPHEREIF_H_
#define _MULTISPHEREIF_H_

#include "MayDay.H"
#include "RealVect.H"
#include "RefCountedPtr.H"

#include "BaseIF.H"
#include "ComplementIF.H"

#include "NamespaceHeader.H"

class RealBox
{
public:
  RealBox()
  {
    m_lo = -RealVect::Unit;
    m_hi =  RealVect::Zero;
  };

  RealBox(const RealVect& a_lo,
          const RealVect& a_hi)
  {
    m_lo = a_lo;
    m_hi = a_hi;
  };

  void operator=(const RealBox& a_realBox)
  {
    m_lo = a_realBox.m_lo;
    m_hi = a_realBox.m_hi;
  };

  virtual ~RealBox()
  {
  };

  virtual void define(const RealVect& a_lo,
                      const RealVect& a_hi)
  {
    m_lo = a_lo;
    m_hi = a_hi;
  };

  bool contains(const RealVect& a_point)
  {
    bool inside = true;

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (a_point[idir] < m_lo[idir] || a_point[idir] > m_hi[idir])
      {
        inside = false;
        break;
      }
    }

    return inside;
  };

  RealVect size() const
  {
    return (m_hi - m_lo);
  };

  const RealVect& getLo() const
  {
    return m_lo;
  };

  RealVect getLo()
  {
    return m_lo;
  };

  const RealVect& getHi() const
  {
    return m_hi;
  };

  RealVect getHi()
  {
    return m_hi;
  };

  void setLo(const RealVect& a_lo)
  {
    m_lo = a_lo;
  };

  void setHi(const RealVect& a_hi)
  {
    m_hi = a_hi;
  };

protected:
  RealVect m_lo;
  RealVect m_hi;
};

class SphereTree
{
public:
  SphereTree(const RealBox&          a_bbox,
             const Vector<Real>&     a_radii,
             const Vector<RealVect>& a_centers);

  virtual ~SphereTree()
  {
    if (m_left != NULL) {
      delete m_left;
    }

    if (m_right != NULL) {
      delete m_right;
    }
  };

  const SphereTree& findNode(const RealVect& a_point) const
  {
    if (m_numSpheres != 0)
    {
      return *this;
    }
    else if (m_left->m_bbox.contains(a_point))
    {
      return m_left->findNode(a_point);
    }
    else if (m_right->m_bbox.contains(a_point))
    {
      return m_right->findNode(a_point);
    }
    else
    {
      pout() << "SphereTree::findNode - " << a_point
             << " not contained in " << m_bbox.getLo()
             << " - " << m_bbox.getHi()
             << endl;

      MayDay::Error("SphereTree::findNode - Couldn't find a node containing the given point");
    }
  };

  const int& getNumSpheres() const
  {
    return m_numSpheres;
  };

  int getNumSpheres()
  {
    return m_numSpheres;
  };

  const Vector<Real>& getRadii() const
  {
    return m_radii;
  };

  Vector<Real> getRadii()
  {
    return m_radii;
  };

  const Vector<RealVect>& getCenters() const
  {
    return m_centers;
  };

  Vector<RealVect> getCenters()
  {
    return m_centers;
  };

protected:
  RealBox          m_bbox;

  int              m_numSpheres;
  Vector<Real>     m_radii;
  Vector<RealVect> m_centers;

  SphereTree*      m_left;
  SphereTree*      m_right;

private:
  SphereTree()
  {
    MayDay::Abort("SphereTree uses strong construction");
  }

  void operator=(const SphereTree& a_inputIF)
  {
    MayDay::Abort("SphereTree doesn't allow assignment");
  }
};

///
/**
    This implicit function specifies a union of spheres.
 */
class MultiSphereIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying sphere radii (a_radii), centers (a_centers), and
      whether the domain is on the inside (a_inside).
   */
  MultiSphereIF(const Vector<Real>&     a_radii,
                const Vector<RealVect>& a_centers,
                const bool&             a_inside);

  MultiSphereIF(const int&                a_numSpheres,
                const bool&               a_inside,
                const RealBox&            a_bbox,
                RefCountedPtr<SphereTree> a_sphereTree);

  /// Destructor
  virtual ~MultiSphereIF();

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM>& a_point) const;

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  void partitionSpace(const Vector<Real>&     a_radii,
                      const Vector<RealVect>& a_centers);

  int              m_numSpheres;  // number of spheres
  bool             m_inside;      // inside flag

  RealBox          m_bbox;

  RefCountedPtr<SphereTree> m_sphereTree;

private:
  MultiSphereIF()
  {
    MayDay::Abort("MultiSphereIF uses strong construction");
  }

  MultiSphereIF(const MultiSphereIF& a_inputIF)
  {
    MayDay::Abort("MultiSphereIF doesn't allow copy construction");
  }

  void operator=(const MultiSphereIF& a_inputIF)
  {
    MayDay::Abort("MultiSphereIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
