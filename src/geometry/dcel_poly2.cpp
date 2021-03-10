/*!
  @file   dcel_poly2.cpp
  @brief  Implementation of dcel_poly2.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_vertex.H"
#include "dcel_face.H"
#include "dcel_poly2.H"
#include "dcel_iterator.H"

using namespace dcel;

Polygon2D::Point2D::Point2D(const Real a_x, const Real a_y){
  x = a_x;
  y = a_y;
}

Polygon2D::Point2D& Polygon2D::Point2D::operator=(const Point2D& a_other) noexcept {
  x = a_other.x;
  y = a_other.y;

  return *this;
}

Polygon2D::Point2D& Polygon2D::Point2D::operator+=(const Point2D& a_other) noexcept {
  x += a_other.x;
  y += a_other.y;

  return *this;
}

Polygon2D::Point2D& Polygon2D::Point2D::operator-=(const Point2D& a_other) noexcept {
  x -= a_other.x;
  y -= a_other.y;

  return *this;
}

Polygon2D::Point2D Polygon2D::Point2D::operator+(const Point2D& a_other) const noexcept {
  return Point2D(x+a_other.x, y+a_other.y);
}

Polygon2D::Point2D Polygon2D::Point2D::operator-(const Point2D& a_other) const noexcept {
  return Point2D(x-a_other.x, y-a_other.y);
}

Real Polygon2D::Point2D::dotProduct(const Point2D& a_other) noexcept{
  return x*a_other.x + y*a_other.y;
}

Real Polygon2D::Point2D::length() const noexcept {
  return sqrt(x*x + y*y);
}

Real Polygon2D::Point2D::length2() const noexcept {
  return x*x + y*y;
}

Polygon2D::Polygon2D(const face& a_face){
  this->define(a_face);
}

bool Polygon2D::isPointInside(const RealVect& a_point, const InsideOutsideAlgorithm a_algorithm) {
  bool ret;
  
  switch(a_algorithm){
  case InsideOutsideAlgorithm::SubtendedAngle:
    ret = this->isPointInsidePolygonSubtend(a_point);
    break;
  case InsideOutsideAlgorithm::CrossingNumber:
    ret = this->isPointInsidePolygonCrossingNumber(a_point);
    break;
  case InsideOutsideAlgorithm::WindingNumber:
    ret = this->isPointInsidePolygonWindingNumber(a_point);
    break;
  default:
    std::cerr << "In file 'dcel_face.cpp' function dcel::Polygon2D::isPointInside - unsupported algorithm requested.\n";
  }

  return ret;
}

bool Polygon2D::isPointInsidePolygonWindingNumber(const RealVect& a_point) const noexcept {
  const Point2D p = this->projectPoint(a_point);
  
  const int wn = this->computeWindingNumber(p);

  return wn != 0;
}

bool Polygon2D::isPointInsidePolygonCrossingNumber(const RealVect& a_point) const noexcept {
  const Point2D p = this->projectPoint(a_point);
  
  const int cn  = this->computeCrossingNumber(p);
    
  const bool ret = (cn&1);

  return ret;
}

bool Polygon2D::isPointInsidePolygonSubtend(const RealVect& a_point) const noexcept {
  const Point2D p = this->projectPoint(a_point);

  Real sumTheta = this->computeSubtendedAngle(p); // Should be = 2pi if point is inside. 

  sumTheta = std::abs(sumTheta)/(2.*M_PI); 

  const bool ret = (round(sumTheta) == 1); // 2PI if the polygon is inside. 

  return ret;
}

bool Polygon2D::isPointOnEdge(const Point2D& a_point, const Point2D& a_endPoint1, const Point2D& a_endPoint2, const Real a_thresh) const noexcept{
  const Real AB2 = (a_endPoint2 - a_endPoint1).length2();
  const Real CA2 = (a_point     - a_endPoint1).length2();
  const Real CB2 = (a_point     - a_endPoint2).length2();

  bool ret = false;
  
  if(CA2 + CB2 <= AB2*(1. + a_thresh)) ret = true;

  return ret;
}

bool Polygon2D::isPointOnBoundary(const Point2D& a_point, const Real a_thresh) const noexcept {
  const int N = m_points.size();

  bool ret = false;
  for (int i = 0; i < N; i++){
    const Point2D& p1 = m_points[i];
    const Point2D& p2 = m_points[(i+1)%N];

    const bool pointOnLine = this->isPointOnEdge(a_point, p1, p2, a_thresh);

    if(pointOnLine){
      ret = true;
      break;
    }
  }

  return ret;
}

Polygon2D::Point2D Polygon2D::projectPoint(const RealVect& a_point) const noexcept {
  return Point2D(a_point[m_xDir], a_point[m_yDir]);
}

void Polygon2D::define(const face& a_face) noexcept {
  m_points.resize(0);

  const RealVect& normal = a_face.getNormal();
  m_ignoreDir = normal.maxDir(true);
  
  // for (int dir = 0; dir < SpaceDim; dir++){
  //   m_ignoreDir = (std::abs(normal[dir]) > std::abs(normal[m_ignoreDir])) ? dir : m_ignoreDir;
  // }

  m_xDir = 3;
  m_yDir = -1;

  for (int dir = 0; dir < SpaceDim; dir++){
    if(dir != m_ignoreDir){
      m_xDir = std::min(m_xDir, dir);
      m_yDir = std::max(m_yDir, dir);
    }
  }

  if(normal[m_ignoreDir] > 0.0){
    int tmp = m_xDir;
    m_xDir = m_yDir;
    m_yDir = tmp;
  }

  // Ignore coordinate with biggest normal component
  for (const auto& v : a_face.getVertices()){
    const RealVect& p = v->getPosition();
    
    m_points.emplace_back(projectPoint(p));
  }
}

int Polygon2D::computeWindingNumber(const Point2D& P) const noexcept {
  int wn = 0;    // the  winding number counter

  const int N = m_points.size();

  auto isLeft = [](const Point2D& P0, const Point2D& P1, const Point2D& P2){
    return (P1.x - P0.x)*(P2.y - P0.y) - (P2.x -  P0.x)*(P1.y - P0.y);
  };

  // loop through all edges of the polygon
  for (int i = 0; i < N; i++) {   // edge from V[i] to  V[i+1]

    const Point2D& P1 = m_points[i];
    const Point2D& P2 = m_points[(i+1)%N];

    const int res = int(isLeft(P1, P2, P));
    
    if (P1.y <= P.y) {          // start y <= P.y
      if (P2.y  > P.y)      // an upward crossing
	if (res > 0)  // P left of  edge
	  ++wn;            // have  a valid up intersect
    }
    else {                        // start y > P.y (no test needed)
      if (P2.y <= P.y)     // a downward crossing
	if (res < 0)  // P right of  edge
	  --wn;            // have  a valid down intersect
    }
  }
  
  return wn;
}

int Polygon2D::computeCrossingNumber(const Point2D& P) const noexcept {
  int cn = 0; 

  const int N = m_points.size();

  constexpr Real thresh = 1.E-6;
  
  for (int i = 0; i < N; i++) {    // edge from V[i]  to V[i+1]
    const Point2D& P1 = m_points[i];
    const Point2D& P2 = m_points[(i+1)%N];

    const bool upwardCrossing   = (P1.y <= P.y) && (P2.y >  P.y);
    const bool downwardCrossing = (P1.y >  P.y) && (P2.y <= P.y);

    if(upwardCrossing || downwardCrossing){
      const Real t = (P.y - P1.y)/(P2.y - P1.y);

      if (P.x <  P1.x + t * (P2.x - P1.x)) {// P.x < intersect
	  cn += 1;   // a valid crossing of y=P.y right of P.x
      }
    }
  }

  return cn;
}

Real Polygon2D::computeSubtendedAngle(const Point2D& p) const noexcept {
  Real sumTheta = 0.0;
  
  const int N = m_points.size();

  constexpr Real thresh = 1.E-6;
  
  for (int i = 0; i < N; i++){
    const Point2D p1 = m_points[i]       - p;
    const Point2D p2 = m_points[(i+1)%N] - p;
    
    const Real theta1 = atan2(p1.y, p1.x);
    const Real theta2 = atan2(p2.y, p2.x);

    Real dTheta = theta2 - theta1;

    while (dTheta > M_PI)
      dTheta -= 2.0*M_PI;
    while (dTheta < -M_PI)
      dTheta += 2.0*M_PI;

    sumTheta += dTheta;
  }

  return sumTheta;
}
