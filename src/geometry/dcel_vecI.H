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

  template <class T>
  inline
  Vec2<T>::Vec2(const T& a_x, const T& a_y){
    x = a_x;
    y = a_y;
  }

  template <class T>
  inline
  Vec2<T>& Vec2<T>::operator=(const Vec2<T>& u) noexcept {
    x = u.x;
    y = u.y;

    return *this;
  }

  template <class T>
  inline
  Vec2<T> Vec2<T>::operator+(const Vec2<T>& u) const noexcept {
    return Vec2<T>(x+u.x, y+u.y);
  }

  template <class T>
  inline
  Vec2<T> Vec2<T>::operator-(const Vec2<T>& u) const noexcept {
    return Vec2<T>(x-u.x, y-u.y);
  }

  template <class T>
  inline
  Vec2<T> Vec2<T>::operator*(const T& s) const noexcept {
    return Vec2<T>(x*s, y*s);
  }

  template <class T>
  inline
  Vec2<T> Vec2<T>::operator/(const T& s) const noexcept {
    const T is = 1./s;
    return Vec2<T>(x*is, y*is);
  }

  template <class T>
  inline
  Vec2<T>& Vec2<T>::operator+=(const Vec2<T>& u) noexcept {
    x += u.x;
    y += u.y;

    return *this;
  }

  template <class T>
  inline
  Vec2<T>& Vec2<T>::operator-=(const Vec2<T>& u) noexcept {
    x -= u.x;
    y -= u.y;

    return *this;
  }

  template <class T>
  inline
  Vec2<T>& Vec2<T>::operator*=(const T& s) noexcept {
    x *= s;
    y *= s;

    return *this;
  }

  template <class T>
  inline
  Vec2<T>& Vec2<T>::operator/=(const T& s) noexcept {
    const T is = 1./s;
    
    x *= is;
    y *= is;

    return *this;
  }

  template <class T>
  inline
  T Vec2<T>::dot(const Vec2<T>& u) const noexcept{
    return x*u.x + y*u.y;
  }

  template <class T>
  inline
  T Vec2<T>::length() const noexcept {
    return sqrt(x*x + y*y);
  }

  template <class T>
  inline
  T Vec2<T>::length2() const noexcept {
    return x*x + y*y;
  }

  template <class T>
  inline
  Vec2<T> operator*(const T& s, const Vec2<T>& a_other) noexcept {
    return a_other*s;
  }

  template <class T>
  inline
  Vec2<T> operator/(const T& s, const Vec2<T>& a_other) noexcept {
    return a_other/s;
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
  Vec3<T> Vec3<T>::operator*(const T& s) const noexcept {
    return Vec3(s*X[0], s*X[1], s*X[2]);
  }

  template <class T>
  inline
  Vec3<T> Vec3<T>::operator/(const T& s) const noexcept {
    const T is = 1./s;
    return Vec3(is*X[0], is*X[1], is*X[2]);
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
  Vec3<T>& Vec3<T>::operator*=(const T& s) noexcept {
    X[0] *= s;
    X[1] *= s;
    X[2] *= s;

    return *this;
  }

  template <class T>
  inline
  Vec3<T>& Vec3<T>::operator/=(const T& s) noexcept {
    const T is = 1./s;
    
    X[0] *= is;
    X[1] *= is;
    X[2] *= is;

    return *this;
  }

  template <class T>
  inline
  Vec3<T> Vec3<T>::cross(const Vec3<T>& u) const noexcept {
    return Vec3<T>(X[1]*u[2] - X[2]*u[1],
		   X[2]*u[0] - X[0]*u[2],
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
  int Vec3<T>::maxDir(const bool a_doAbs) const noexcept {
    int mDir = 0;
    
    for (int dir = 0; dir < 3; dir++){
      if (a_doAbs){
	if (std::abs(X[dir]) > std::abs(X[mDir])) {
	  mDir = dir;
	}
      }
      else{
	if (X[dir] > X[mDir]){
	  mDir = dir;
	}
      }
    }
    
    return mDir;
  }

  template <class T>
  inline
  bool Vec3<T>::operator==(const Vec3<T>& u) const noexcept {
    return (X[0]==u[0] && X[1]==u[1] && X[2]==u[2]);
  }

  template <class T>
  inline
  bool Vec3<T>::operator<(const Vec3<T>& u) const noexcept {
    return (X[0]<u[0] && X[1]<u[1] && X[2]<u[2]);
  }

  template <class T>
  inline
  bool Vec3<T>::operator>(const Vec3<T>& u) const noexcept {
    return (X[0]>u[0] && X[1]>u[1] && X[2]>u[2]);
  }

  template <class T>
  inline
  bool Vec3<T>::operator<=(const Vec3<T>& u) const noexcept {
    return (X[0]<=u[0] && X[1]<=u[1] && X[2]<=u[2]);
  }

  template <class T>
  inline
  bool Vec3<T>::operator>=(const Vec3<T>& u) const noexcept {
    return (X[0]>=u[0] && X[1]>=u[1] && X[2]>=u[2]);
  }

  template <class T>
  inline
  T Vec3<T>::dot(const Vec3<T>& u) const noexcept {
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

  template <class T>
  inline
  Vec3<T> operator*(const T& s, const Vec3<T>& a_other) noexcept {
    return a_other*s;
  }

  template <class T>
  inline
  Vec3<T> operator/(const T& s, const Vec3<T>& a_other) noexcept {
    return a_other/s;
  }
}

#endif
