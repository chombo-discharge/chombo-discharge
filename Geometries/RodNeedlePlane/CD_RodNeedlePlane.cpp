/* chombo-discharge
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RodNeedlePlane.cpp
  @brief  Implementation pf CD_RodNeedlePlane.H
  @author Fanny Skirbekk
*/

// Chombo includes
#include <ParmParse.H>
#include <SmoothUnion.H>
#include <SmoothIntersection.H>

// Our includes
#include <CD_RodNeedlePlane.H>
#include <CD_RodIF.H>
#include <CD_NeedleIF.H>
#include <CD_RoundedCylinderIF.H>
#include <CD_NamespaceHeader.H>

RodNeedlePlane::RodNeedlePlane()
{
  this->setGasPermittivity(1.0);

  ParmParse ppRod("RodNeedlePlane.rod");
  ParmParse ppNeedle("RodNeedlePlane.needle");
  ParmParse ppPlane("RodNeedlePlane.plane");

  bool useRod;
  bool useNeedle;
  bool usePlane;

  ppRod.get("on", useRod);
  ppNeedle.get("on", useNeedle);
  ppPlane.get("on", usePlane);

  // Choose electrodes
  if (useRod && useNeedle)
    this->defineRodWNeedle();
  else if (useRod)
    this->defineRod();
  else if (useNeedle)
    this->defineNeedle();
  if (usePlane)
    this->definePlane();
}

RodNeedlePlane::~RodNeedlePlane() {}

void
RodNeedlePlane::defineRodWNeedle()
{
  Vector<Real> v(SpaceDim);
  Real         rodRadius;
  RealVect     rodE1, rodE2;
  bool         rodLive;

  Real length;
  int  tipDir;
  Real needleRadius;
  Real tipRadius;
  Real angle;
  bool needleLive;

  Real cornerCurve;

  ParmParse ppRod("RodNeedlePlane.rod");
  ParmParse ppNeedle("RodNeedlePlane.needle");

  ppRod.get("radius", rodRadius);
  ppRod.get("live", rodLive);

  ppNeedle.get("length", length);
  ppNeedle.get("radius", needleRadius);
  ppNeedle.get("tipRadius", tipRadius);
  ppNeedle.get("angle", angle);
  ppNeedle.get("live", needleLive);

  ppNeedle.get("cornerCurve", cornerCurve);

  ppRod.getarr("endpoint1", v, 0, SpaceDim);
  rodE1 = RealVect(D_DECL(v[0], v[1], v[2]));
  ppRod.getarr("endpoint2", v, 0, SpaceDim);
  rodE2 = RealVect(D_DECL(v[0], v[1], v[2]));

  Vector<BaseIF*> electrodeParts;
  electrodeParts.push_back((BaseIF*)new RodIF(rodE1, rodE2, rodRadius, false));
  electrodeParts.push_back((BaseIF*)new NeedleIF(length, needleRadius, false, tipRadius, angle, cornerCurve));

  RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF>(new SmoothIntersection(electrodeParts, cornerCurve));

  m_electrodes.push_back(Electrode(bif, (rodLive && needleLive)));

  for (int i = 0; i < electrodeParts.size(); ++i) {
    delete electrodeParts[i];
  }
}

void
RodNeedlePlane::defineRod()
{
  Vector<Real> v(SpaceDim);
  Real         radius;
  RealVect     e1, e2;
  bool         live;

  ParmParse pp("RodNeedlePlane.rod");

  pp.get("radius", radius);
  pp.get("live", live);

  pp.getarr("endpoint1", v, 0, SpaceDim);
  e1 = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("endpoint2", v, 0, SpaceDim);
  e2 = RealVect(D_DECL(v[0], v[1], v[2]));

  RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF>(new RodIF(e1, e2, radius, false));

  m_electrodes.push_back(Electrode(bif, live));
}

void
RodNeedlePlane::defineNeedle()
{
  Vector<Real> v(SpaceDim);
  Real         length;
  int          tipDir;
  Real         radius;
  Real         tipRadius;
  Real         angle;
  Real         cornerCurve;
  bool         live;

  ParmParse pp("RodNeedlePlane.needle");

  pp.get("length", length);
  pp.get("radius", radius);
  pp.get("tipRadius", tipRadius);
  pp.get("angle", angle);
  pp.get("cornerCurve", cornerCurve);
  pp.get("live", live);

  RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF>(new NeedleIF(length, radius, false, tipRadius, angle, cornerCurve));

  m_electrodes.push_back(Electrode(bif, live));
}

void
RodNeedlePlane::definePlane()
{
  Vector<Real> v(SpaceDim);
  Real         radius;
  Real         thickness;
  Real         curve;
  RealVect     normal, point;
  bool         live;

  ParmParse pp("RodNeedlePlane.plane");

  pp.get("radius", radius);
  pp.get("thickness", thickness);
  pp.get("curve", curve);
  pp.get("live", live);
  pp.getarr("normal", v, 0, SpaceDim);
  normal = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("point", v, 0, SpaceDim);
  point = RealVect(D_DECL(v[0], v[1], v[2]));

  RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF>(
    new RoundedCylinderIF(point, point - thickness * normal, radius, curve, false));

  m_electrodes.push_back(Electrode(bif, live));
}

#include <CD_NamespaceFooter.H>
