/*!
  @file   dcel_vec.H
  @brief  Declaration of 2D and 3D point/vector classes (so that we can decouple from Chombo). 
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_VEC_H_
#define _DCEL_VEC_H_

#include <memory>
#include <vector>

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
    Vec3T& operator=(const Vec3T& u) noexcept;

    inline
    Vec3T operator+(const Vec3T& u) const noexcept;

    inline
    Vec3T operator-(const Vec3T& u) const noexcept;

    inline
    Vec3T operator*(const T& s) const noexcept;

    inline
    Vec3T operator/(const T& s) const noexcept;

    inline
    Vec3T& operator+=(const Vec3T& u) noexcept;

    inline
    Vec3T& operator-=(const Vec3T& u) noexcept;

    inline
    Vec3T& operator*=(const T& s) noexcept;

    inline
    Vec3T& operator/=(const T& s) noexcept;

    inline
    T& operator[](int i) noexcept;

    inline
    const T& operator[](int i) const noexcept;

    inline
    int maxDir(const bool a_doAbs) const noexcept;

    inline
    bool operator==(const Vec3T& u) const noexcept;

    inline
    bool operator<(const Vec3T& u) const noexcept;

    inline
    bool operator>(const Vec3T& u) const noexcept;

    inline
    bool operator<=(const Vec3T& u) const noexcept;

    inline
    bool operator>=(const Vec3T& u) const noexcept;

    inline
    Vec3T cross(const Vec3T& u) const noexcept;

    inline
    T dot(const Vec3T& u) const noexcept;

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
}

#include "dcel_vecI.H"

#endif
