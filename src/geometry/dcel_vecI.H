/*!
  @file   dcel_vecI.H
  @brief  Implementation of dcel_vec.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_VECI_H_
#define _DCEL_VECI_H_

#include "dcel_vec.H"

#include <math.h>

namespace dcel {

  inline
  Vec2::Vec2(const double& a_x, const double& a_y){
    x = a_x;
    y = a_y;
  }

  inline
  Vec2& Vec2::operator=(const Vec2& u) noexcept {
    x = u.x;
    y = u.y;

    return *this;
  }

  inline
  Vec2& Vec2::operator+=(const Vec2& u) noexcept {
    x += u.x;
    y += u.y;

    return *this;
  }

  inline
  Vec2& Vec2::operator-=(const Vec2& u) noexcept {
    x -= u.x;
    y -= u.y;

    return *this;
  }

  inline
  Vec2 Vec2::operator+(const Vec2& u) const noexcept {
    return Vec2(x+u.x, y+u.y);
  }

  inline
  Vec2 Vec2::operator-(const Vec2& u) const noexcept {
    return Vec2(x-u.x, y-u.y);
  }

  inline
  double Vec2::dot(const Vec2& u) noexcept{
    return x*u.x + y*u.y;
  }

  inline
  double Vec2::length() const noexcept {
    return sqrt(x*x + y*y);
  }

  inline
  double Vec2::length2() const noexcept {
    return x*x + y*y;
  }

  template <class T>
  inline
  Vec3<T>::Vec3(const T& a_x, const T& a_y, const T& a_z){
    X[0] = a_x;
    X[1] = a_y;
    X[2] = a_z;
  }

  template <class T>
  inline
  Vec3<T>& Vec3<T>::operator=(const Vec3<T>& u) noexcept {
    X[0] = u[0];
    X[1] = u[1];
    X[2] = u[2];

    return *this;
  }

  template <class T>
  inline
  Vec3<T>& Vec3<T>::operator+=(const Vec3<T>& u) noexcept {
    X[0] += u[0];
    X[1] += u[1];
    X[2] += u[2];

    return *this;
  }

  template <class T>
  inline
  Vec3<T>& Vec3<T>::operator-=(const Vec3<T>& u) noexcept {
    X[0] -= u[0];
    X[1] -= u[1];
    X[2] -= u[2];

    return *this;
  }

  template <class T>
  inline
  Vec3<T> Vec3<T>::operator+(const Vec3<T>& u) const noexcept {
    return Vec3<T>(X[0]+u[0], X[1]+u[1], X[2]+u[2]);
  }

  template <class T>
  inline
  Vec3<T> Vec3<T>::operator-(const Vec3<T>& u) const noexcept {
    return Vec3<T>(X[0]-u[0], X[1]-u[1], X[2]-u[2]);
  }

  template <class T>
  inline
  Vec3<T> Vec3<T>::cross(const Vec3<T>& u) const noexcept {
    return Vec3<T>(X[1]*u[2] - X[2]*u[1],
		X[2]*u[1] - X[1]*u[2],
		X[0]*u[1] - X[1]*u[0]);
  }

  template <class T>
  inline
  T& Vec3<T>::operator[](int i) noexcept {
    return X[i];
  }

  template <class T>
  inline
  const T& Vec3<T>::operator[](int i) const noexcept {
    return X[i];
  }

  template <class T>
  inline
  T Vec3<T>::dot(const Vec3<T>& u) noexcept {
    return X[0]*u[0] + X[1]*u[1] + X[2]*u[2];
  }

  template <class T>
  inline
  T Vec3<T>::length() const noexcept {
    return sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
  }

  template <class T>
  inline
  T Vec3<T>::length2() const noexcept {
    return X[0]*X[0] + X[1]*X[1] + X[2]*X[2];
  }
}

#endif
