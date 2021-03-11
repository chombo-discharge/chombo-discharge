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

    static const Vec2T Zero;
    static const Vec2T Unit;

    Vec2T() = default;
    Vec2T(const T& a_x, const T& a_y);
    ~Vec2T() = default;
    
    T x;
    T y;

    inline
    auto operator=(const Vec2T& a_other) noexcept -> Vec2T<T>&;

    inline
    auto operator+(const Vec2T& a_other) const noexcept -> Vec2T<T>;

    inline
    auto operator-(const Vec2T& a_other) const noexcept -> Vec2T<T>;

    inline
    auto operator*(const T& s) const noexcept -> Vec2T<T>;

    inline
    auto operator/(const T& s) const noexcept -> Vec2T<T>;

    inline
    auto operator+=(const Vec2T& a_other) noexcept -> Vec2T<T>&;

    inline
    auto operator-=(const Vec2T& a_other) noexcept -> Vec2T<T>&;

    inline
    auto operator*=(const T& s) noexcept -> Vec2T<T>&;

    inline
    auto operator/=(const T& s) noexcept -> Vec2T<T>&;

    inline
    auto dot(const Vec2T& a_other) const noexcept -> T;

    inline
    auto length() const noexcept -> T;

    inline
    auto length2() const noexcept ->T;
  };

  template <class T>
  inline
  auto operator*(const T& s, const Vec2T<T>& a_other) noexcept -> Vec2T<T>;

  template <class T>
  inline
  auto operator/(const T& s, const Vec2T<T>& a_other) noexcept -> Vec2T<T>;

  template <class T>
  class Vec3T {
  public:

    static const Vec3T Zero;
    static const Vec3T Unit;
    
    Vec3T() = default;
    Vec3T(const T& a_x, const T& a_y, const T& a_z);
    ~Vec3T() = default;

    inline
    auto operator=(const Vec3T& u) noexcept -> Vec3T<T>&;

    inline
    auto operator+(const Vec3T& u) const noexcept -> Vec3T<T>;

    inline
    auto operator-(const Vec3T& u) const noexcept -> Vec3T<T>;

    inline
    auto operator*(const T& s) const noexcept -> Vec3T<T>;

    inline
    auto operator/(const T& s) const noexcept -> Vec3T<T>;

    inline
    auto operator+=(const Vec3T& u) noexcept -> Vec3T<T>&;

    inline
    auto operator-=(const Vec3T& u) noexcept -> Vec3T<T>&;

    inline
    auto operator*=(const T& s) noexcept -> Vec3T<T>&;

    inline
    auto operator/=(const T& s) noexcept -> Vec3T<T>&;

    inline
    auto operator[](int i) noexcept -> T&;

    inline
    auto operator[](int i) const noexcept -> const T& ;

    inline
    auto operator==(const Vec3T& u) const noexcept -> bool;

    inline
    auto operator<(const Vec3T& u) const noexcept -> bool;

    inline
    auto operator>(const Vec3T& u) const noexcept -> bool;

    inline
    auto operator<=(const Vec3T& u) const noexcept -> bool;

    inline
    auto operator>=(const Vec3T& u) const noexcept -> bool;

    inline
    auto maxDir(const bool a_doAbs) const noexcept -> int;

    inline
    auto cross(const Vec3T& u) const noexcept -> Vec3T<T>;

    inline
    auto dot(const Vec3T& u) const noexcept -> T;

    inline
    auto length() const noexcept -> T;

    inline
    auto length2() const noexcept -> T;

  protected:
    T X[3];
  };


  template <class T>
  inline
  auto operator*(const T& s, const Vec3T<T>& a_other) noexcept -> Vec3T<T>;

  template <class T>
  inline
  auto operator/(const T& s, const Vec3T<T>& a_other) noexcept -> Vec3T<T>;
}

#include "dcel_vecI.H"

#endif
