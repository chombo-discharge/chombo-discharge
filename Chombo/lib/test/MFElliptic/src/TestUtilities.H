#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TESTUTILITIES_H_
#define _TESTUTILITIES_H_

#include "SphereIF.H"
#include "PlaneIF.H"
#include "BaseBCValue.H"
#include "MFCellFAB.H"
#include "MFIndexSpace.H"
#include "GeometryShop.H"
#include "PolyGeom.H"

#include "UsingNamespace.H"

class ConstantDirichletBC: public BaseBCValue
{
public:
  ConstantDirichletBC(Real a_constant):
    m_constant(a_constant)
  {}

  Real value (const RealVect &a_point,
              const RealVect &a_n,
              const Real&     a_time,
              const int&      a_comp) const;

private:
  Real m_constant;
};

class TrilinearDirichletBC: public BaseBCValue
{
public:
  TrilinearDirichletBC(Real     a_constant,
                       RealVect a_linearCoefficients):
    m_constant(a_constant),
    m_linearCoefficients(a_linearCoefficients)
  {}

  Real value (const RealVect &a_point,
              const RealVect &a_n,
              const Real&     a_time,
              const int&      a_comp) const;

private:
  Real m_constant;
  RealVect m_linearCoefficients;
};

class TriquadraticDirichletBC: public BaseBCValue
{
public:
  TriquadraticDirichletBC(Real     a_constant,
                          RealVect a_linearCoefficients,
                          RealVect a_xxCoefficients,
                          RealVect a_xyCoefficients):
    m_constant(a_constant),
    m_linearCoefficients(a_linearCoefficients),
    m_xxCoefficients(a_xxCoefficients),
    m_xyCoefficients(a_xyCoefficients)
  {}

  Real value (const RealVect &a_point,
              const RealVect &a_n,
              const Real&     a_time,
              const int&      a_comp) const;

private:
  Real m_constant;
  RealVect m_linearCoefficients;
  RealVect m_xxCoefficients;
  RealVect m_xyCoefficients;
};

extern
int readDomainInfo(int&       a_nlevels,
                   Box&       a_domain,
                   Real&      a_dx,
                   RealVect&  a_origin,
                   int&       a_maxBoxSize);

extern
void readPlaneInfo(RealVect& a_normal,
                   RealVect& a_point);

extern
void readSphereInfo(RealVect& a_center,
                    Real&     a_radius);

extern
void makeHierarchy(      Vector<DisjointBoxLayout>& a_dbl,
                   const ProblemDomain&             a_baseDomain,
                   const IntVectSet&                a_baseTags,
                   const Vector<int>&               a_refRatio,
                   const int                        a_nlevels,
                         int                        a_maxSize);

extern
void zeroAdjacentToBoundary(LevelData<MFCellFAB>& a_data,
                            const ProblemDomain&  a_domain);

extern
void setValue(      LevelData<EBCellFAB>& a_phase,
              const BaseBCValue&          a_bc,
              const Box&                  a_domain,
              const RealVect&             a_dx,
              const RealVect&             a_origin,
                    bool                  a_useKappa = false);

#endif
