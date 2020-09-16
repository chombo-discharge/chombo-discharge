#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SphereIF.H"
#include "IntersectionIF.H"
#include "MultiSphereIF.H"
#include "PolyGeom.H"

#include "NamespaceHeader.H"

MultiSphereIF::MultiSphereIF(const Vector<Real>&     a_radii,
                             const Vector<RealVect>& a_centers,
                             const bool&             a_inside)
{
  CH_TIME("MultiSphereIF::MultiSphereIF");

  // Remember the number of spheres
  m_numSpheres = a_radii.size();

  // Remember the parameters
  m_inside  = a_inside;

  partitionSpace(a_radii,a_centers);
}

MultiSphereIF::MultiSphereIF(const int&                a_numSpheres,
                             const bool&               a_inside,
                             const RealBox&            a_bbox,
                             RefCountedPtr<SphereTree> a_sphereTree)
{
  CH_TIME("MultiSphereIF::MultiSphereIF2");

  // Remember the number of spheres
  m_numSpheres = a_numSpheres;

  // Remember the parameters
  m_inside  = a_inside;

  m_bbox       = a_bbox;
  m_sphereTree = a_sphereTree;
}

MultiSphereIF::~MultiSphereIF()
{
}

Real MultiSphereIF::value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                          const IndexTM<Real,GLOBALDIM>& a_point) const
{
  bool found = false;

  RealVect nearPoint = a_point;

  nearPoint.max(m_bbox.getLo());
  nearPoint.min(m_bbox.getHi());
#if 0
  pout() << "value: " << a_point << " -> " << nearPoint << endl;
#endif
  const SphereTree& node = m_sphereTree->findNode(nearPoint);

  int numSpheres = node.getNumSpheres();

  const Vector<RealVect>& centers = node.getCenters();
  const Vector<Real>&     radii   = node.getRadii  ();
#if 0
  pout() << "  # = " << numSpheres << endl;
#endif
  Real minValue = 1.0;
  int whichSphere = -1;

  for (int isphere = 0; isphere < numSpheres; isphere++)
  {
    RealVect vectDist = a_point - centers[isphere];
    Real curValue = PolyGeom::dot(vectDist,vectDist) - radii[isphere]*radii[isphere];

    if (isphere == 0 || (curValue < minValue))
    {
      found = true;

      minValue = curValue;
      whichSphere = isphere;
    }
  }
#if 0
  pout() << "  which = " << whichSphere << ", min = " << minValue << endl;
#endif
  if (!found)
  {
    MayDay::Error("logic error in MultiSphere::value");
  }

  RealVect center = centers[whichSphere];

  Real retval = minValue;
  int maxDir  = a_partialDerivative.maxDir(false);

  int derivativeOrder = a_partialDerivative.sum();

  if (derivativeOrder == 0)
  {
    retval = minValue;
  }
  else if (derivativeOrder == 1)
  {
    retval = 2.0*(a_point[maxDir] - center[maxDir]);
  }
  else if (derivativeOrder == 2)
  {
    if (a_partialDerivative[maxDir] == 2)
    {
      // unmixed second partial = 2.0
      retval = 2.0;
    }
    else
    {
      // mixed partials = 0.0
      retval = 0.0;
    }
  }
  else
  {
    // higher partials = 0.0
    retval = 0.0;
  }

  // Change the sign to change inside to outside
  if (!m_inside && derivativeOrder > 0)
  {
    retval = -retval;
  }

  return retval;
}

Real MultiSphereIF::value(const RealVect& a_point) const
{
  Real retval;

  IndexTM<Real,GLOBALDIM> point;

  for (int idir = 0; idir < GLOBALDIM; idir++)
  {
    if (idir < SpaceDim)
    {
      point[idir] = a_point[idir];
    }
    else
    {
      point[idir] = 0.0;
    }
  }

  retval = value(IndexTM<int,GLOBALDIM>::Zero,point);

  return retval;
}

BaseIF* MultiSphereIF::newImplicitFunction() const
{
  CH_TIME("MultiSphereIF::newImplicitFunction");

  MultiSphereIF* spherePtr = new MultiSphereIF(m_numSpheres,
                                               m_inside,
                                               m_bbox,
                                               m_sphereTree);

  return static_cast<BaseIF*>(spherePtr);
}

void MultiSphereIF::partitionSpace(const Vector<Real>&     a_radii,
                                   const Vector<RealVect>& a_centers)
{
  CH_TIME("MultiSphereIF::partitionSpace");

  RealVect lo = RealVect::Zero;
  RealVect hi = RealVect::Zero;

  for (int i = 0; i < m_numSpheres; i++)
  {
    RealVect r = a_radii[i]*RealVect::Unit;
    RealVect c = a_centers[i];

    if (i == 0)
    {
      lo = c - r;
      hi = c + r;
    }
    else
    {
      lo.min(c-r);
      hi.max(c+r);
    }
  }

  m_bbox.define(lo,hi);

  RefCountedPtr<SphereTree> curTree(new SphereTree(m_bbox,a_radii,a_centers));

  m_sphereTree = curTree;
}

SphereTree::SphereTree(const RealBox&          a_bbox,
                       const Vector<Real>&     a_radii,
                       const Vector<RealVect>& a_centers)
{
  RealVect boxSize = a_bbox.size();
  int numSpheres = a_radii.size();
#if 0
  pout() << "SphereTree: " << a_bbox.getLo() 
         << " - " << a_bbox.getHi()
         << ", # = " << a_centers.size();
#endif
  m_bbox = a_bbox;

  m_numSpheres = 0;

  m_left = NULL;
  m_right = NULL;

  if (numSpheres <= 100 || boxSize[boxSize.maxDir(false)] < a_radii[0])
  {
#if 0
    pout() << " ==> Leaf" << endl;
#endif
    m_numSpheres = numSpheres;
    m_radii      = a_radii;
    m_centers    = a_centers;
  }
  else
  {
#if 0
    pout() << " ==> Node" << endl;
#endif
    int maxDir = boxSize.maxDir(false);

    RealVect bloat = 2 * a_radii[0] * RealVect::Unit;
    
    RealVect leftLo = a_bbox.getLo();
    RealVect leftHi = a_bbox.getHi();
    leftHi[maxDir] = leftLo[maxDir] + boxSize[maxDir]/2.0;

    RealBox leftBox (leftLo ,leftHi );

    leftLo -= bloat;
    leftHi += bloat;

    RealBox leftBoxBloated (leftLo ,leftHi );

    Vector<Real>     leftRadii;
    Vector<RealVect> leftCenters;

    RealVect rightLo = a_bbox.getLo();
    rightLo[maxDir] = rightLo[maxDir] + boxSize[maxDir]/2.0;
    RealVect rightHi = a_bbox.getHi();

    RealBox rightBox(rightLo,rightHi);

    rightLo -= bloat;
    rightHi += bloat;

    RealBox rightBoxBloated (rightLo ,rightHi );

    Vector<Real>     rightRadii;
    Vector<RealVect> rightCenters;

    for (int i = 0; i < numSpheres; i++)
    {
      const RealVect& curCenter = a_centers[i];
      const Real&     curRadius = a_radii  [i];

      if (leftBoxBloated.contains(curCenter))
      {
        leftRadii  .push_back(curRadius);
        leftCenters.push_back(curCenter);
      }

      if (rightBoxBloated.contains(curCenter))
      {
        rightRadii  .push_back(curRadius);
        rightCenters.push_back(curCenter);
      }
    }

    m_left  = new SphereTree(leftBox ,leftRadii ,leftCenters );
    m_right = new SphereTree(rightBox,rightRadii,rightCenters);
  }
}
#include "NamespaceFooter.H"
