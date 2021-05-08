/*!
  @file   vec.H
  @brief  Declaration of 2D and 3D point/vector classes with templated precision (so that we can decouple from Chombo). 
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _VEC_H_
#define _VEC_H_

#include "CD_NamespaceHeader.H"

template<class T>
class Vec2T {
public:

  Vec2T();
  Vec2T(const Vec2T& u);
  Vec2T(const T& a_x, const T& a_y);
  ~Vec2T() = default;
    
  T x;
  T y;

  inline
  static Vec2T<T> zero() noexcept;

  inline
  static Vec2T<T> one() noexcept;

  inline
  static Vec2T<T> min() noexcept;

  inline
  static Vec2T<T> max() noexcept;

  inline
  static Vec2T<T> infinity() noexcept;

  inline
  Vec2T<T>& operator=(const Vec2T& a_other) noexcept;

  inline
  Vec2T<T> operator+(const Vec2T& a_other) const noexcept;

  inline
  Vec2T<T> operator-(const Vec2T& a_other) const noexcept;
    
  inline
  Vec2T<T> operator-() const noexcept;

  inline
  Vec2T<T> operator*(const T& s) const noexcept;

  inline
  Vec2T<T> operator/(const T& s) const noexcept;

  inline
  Vec2T<T>& operator+=(const Vec2T& a_other) noexcept;

  inline
  Vec2T<T>& operator-=(const Vec2T& a_other) noexcept;

  inline
  Vec2T<T>& operator*=(const T& s) noexcept;

  inline
  Vec2T<T>& operator/=(const T& s) noexcept;

  inline
  T dot(const Vec2T& a_other) const noexcept;

  inline
  T length() const noexcept;

  inline
  T length2() const noexcept;
};

template <class T>
inline
Vec2T<T> operator*(const T& s, const Vec2T<T>& a_other) noexcept;

template <class T>
inline
Vec2T<T> operator/(const T& s, const Vec2T<T>& a_other) noexcept;

template <class T>
class Vec3T {
public:
    
  Vec3T();

    
  Vec3T(const Vec3T<T>& a_u);
  Vec3T(const T& a_x, const T& a_y, const T& a_z);
  ~Vec3T() = default;

  inline
  static Vec3T<T> zero() noexcept;

  inline
  static Vec3T<T> one() noexcept;

  inline
  static Vec3T<T> min() noexcept;

  inline
  static Vec3T<T> max() noexcept;

  inline
  static Vec3T<T> infinity() noexcept;

  inline
  T& operator[](int i) noexcept;

  inline
  const T& operator[](int i) const noexcept;
    
  inline
  bool operator==(const Vec3T<T>& u) const noexcept;
    
  inline
  bool operator<(const Vec3T<T>& u) const noexcept;
    
  inline
  bool operator>(const Vec3T<T>& u) const noexcept;
    
  inline
  bool operator<=(const Vec3T<T>& u) const noexcept;
    
  inline
  bool operator>=(const Vec3T<T>& u) const noexcept;
    
  inline
  Vec3T<T>& operator=(const Vec3T<T>& u) noexcept;
    
  inline
  Vec3T<T> operator+(const Vec3T<T>& u) const noexcept;
    
  inline
  Vec3T<T> operator-(const Vec3T<T>& u) const noexcept;

  inline
  Vec3T<T> operator-() const noexcept;
    
  inline
  Vec3T<T> operator*(const T& s) const noexcept;
    
  inline
  Vec3T<T> operator/(const T& s) const noexcept;
    
  inline
  Vec3T<T>& operator+=(const Vec3T<T>& u) noexcept;
    
  inline
  Vec3T<T>& operator-=(const Vec3T<T>& u) noexcept;
    
  inline
  Vec3T<T>& operator*=(const T& s) noexcept;
    
  inline
  Vec3T<T>& operator/=(const T& s) noexcept;
    
  inline
  Vec3T<T> min(const Vec3T<T>& u) noexcept;
    
  inline
  Vec3T<T> max(const Vec3T<T>& u) noexcept;
    
  inline
  Vec3T<T> cross(const Vec3T<T>& u) const noexcept;
    
  inline
  T dot(const Vec3T<T>& u) const noexcept;

  inline
  int maxDir(const bool a_doAbs) const noexcept;

  inline
  T length() const noexcept;

  inline
  T length2() const noexcept;

protected:
  T X[3];
};

template <class T>
inline
Vec3T<T> operator*(const T& s, const Vec3T<T>& a_other) noexcept;

template <class T>
inline
Vec3T<T> operator/(const T& s, const Vec3T<T>& a_other) noexcept;

template <class T>
inline
Vec3T<T> min(const Vec3T<T>& u, const Vec3T<T>& v) noexcept;

template <class T>
inline
Vec3T<T> max(const Vec3T<T>& u, const Vec3T<T>& v) noexcept;
#include "CD_NamespaceFooter.H"

#include "VecI.H"

#endif
