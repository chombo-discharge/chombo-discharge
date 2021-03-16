/*!
  @file   dcel_vec.H
  @brief  Declaration of 2D and 3D point/vector classes (so that we can decouple from Chombo). 
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_VEC_H_
#define _DCEL_VEC_H_

namespace dcel {

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
    static decltype(auto) zero() noexcept;

    inline
    static decltype(auto) one() noexcept;

    inline
    decltype(auto) operator=(const Vec2T& a_other) noexcept;

    inline
    decltype(auto) operator+(const Vec2T& a_other) const noexcept;

    inline
    decltype(auto) operator-(const Vec2T& a_other) const noexcept;

    inline
    decltype(auto) operator*(const T& s) const noexcept;

    inline
    decltype(auto) operator/(const T& s) const noexcept;

    inline
    decltype(auto) operator+=(const Vec2T& a_other) noexcept;

    inline
    decltype(auto) operator-=(const Vec2T& a_other) noexcept;

    inline
    decltype(auto) operator*=(const T& s) noexcept;

    inline
    decltype(auto) operator/=(const T& s) noexcept;

    inline
    decltype(auto) dot(const Vec2T& a_other) const noexcept;

    inline
    decltype(auto) length() const noexcept;

    inline
    decltype(auto) length2() const noexcept;
  };

  template <class T>
  inline
  decltype(auto) operator*(const T& s, const Vec2T<T>& a_other) noexcept;

  template <class T>
  inline
  decltype(auto) operator/(const T& s, const Vec2T<T>& a_other) noexcept;

  template <class T>
  class Vec3T {
  public:
    
    Vec3T();
    Vec3T(const Vec3T& a_u);
    Vec3T(const T& a_x, const T& a_y, const T& a_z);
    ~Vec3T() = default;

    inline
    static decltype(auto) zero() noexcept;

    inline
    static decltype(auto) one() noexcept;

    inline
    decltype(auto) operator=(const Vec3T& u) noexcept;

    inline
    decltype(auto) operator+(const Vec3T& u) const noexcept;

    inline
    decltype(auto) operator-(const Vec3T& u) const noexcept;

    inline
    decltype(auto) operator*(const T& s) const noexcept;

    inline
    decltype(auto) operator/(const T& s) const noexcept;

    inline
    decltype(auto) operator+=(const Vec3T& u) noexcept;

    inline
    decltype(auto) operator-=(const Vec3T& u) noexcept;

    inline
    decltype(auto) operator*=(const T& s) noexcept;

    inline
    decltype(auto) operator/=(const T& s) noexcept;

    inline
    decltype(auto) operator[](int i) noexcept;

    inline
    decltype(auto) operator[](int i) const noexcept;

    inline
    decltype(auto) operator==(const Vec3T& u) const noexcept;

    inline
    decltype(auto) operator<(const Vec3T& u) const noexcept;

    inline
    decltype(auto) operator>(const Vec3T& u) const noexcept;

    inline
    decltype(auto) operator<=(const Vec3T& u) const noexcept;

    inline
    decltype(auto) operator>=(const Vec3T& u) const noexcept;

    inline
    decltype(auto) min(const Vec3T& u) noexcept;

    inline
    decltype(auto) max(const Vec3T& u) noexcept;

    inline
    decltype(auto) maxDir(const bool a_doAbs) const noexcept;

    inline
    decltype(auto) cross(const Vec3T& u) const noexcept;

    inline
    decltype(auto) dot(const Vec3T& u) const noexcept;

    inline
    decltype(auto) length() const noexcept;

    inline
    decltype(auto) length2() const noexcept;

  protected:
    T X[3];
  };

  template <class T>
  inline
  decltype(auto) operator*(const T& s, const Vec3T<T>& a_other) noexcept;

  template <class T>
  inline
  decltype(auto) operator/(const T& s, const Vec3T<T>& a_other) noexcept;

  template <class T>
  inline
  decltype(auto) min(const Vec3T<T>& u, const Vec3T<T>& v) noexcept;

  template <class T>
  inline
  decltype(auto) max(const Vec3T<T>& u, const Vec3T<T>& v) noexcept;
}

#include "dcel_vecI.H"

#endif
