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
  class Vec2 {
  public:

    static const Vec2 Zero;
    static const Vec2 Unit;

    Vec2() = default;
    Vec2(const T& a_x, const T& a_y);
    ~Vec2() = default;
    
    T x;
    T y;

    inline
    Vec2& operator=(const Vec2& a_other) noexcept;


    inline
    Vec2 operator+(const Vec2& a_other) const noexcept;

    inline
    Vec2 operator-(const Vec2& a_other) const noexcept;

    inline
    Vec2 operator*(const T& s) const noexcept;

    inline
    Vec2 operator/(const T& s) const noexcept;

    inline
    Vec2& operator+=(const Vec2& a_other) noexcept;

    inline
    Vec2& operator-=(const Vec2& a_other) noexcept;

    inline
    Vec2& operator*=(const T& s) noexcept;

    inline
    Vec2& operator/=(const T& s) noexcept;

    inline
    T dot(const Vec2& a_other) const noexcept;

    inline
    T length() const noexcept;

    inline
    T length2() const noexcept;
  };

  template <class T>
  inline
  Vec2<T> operator*(const T& s, const Vec2<T>& a_other) noexcept;

  template <class T>
  inline
  Vec2<T> operator/(const T& s, const Vec2<T>& a_other) noexcept;

  template <class T>
  class Vec3 {
  public:

    static const Vec3 Zero;
    static const Vec3 Unit;
    
    Vec3() = default;
    Vec3(const T& a_x, const T& a_y, const T& a_z);
    ~Vec3() = default;

    inline
    Vec3& operator=(const Vec3& u) noexcept;

    inline
    Vec3 operator+(const Vec3& u) const noexcept;

    inline
    Vec3 operator-(const Vec3& u) const noexcept;

    inline
    Vec3 operator*(const T& s) const noexcept;

    inline
    Vec3 operator/(const T& s) const noexcept;

    inline
    Vec3& operator+=(const Vec3& u) noexcept;

    inline
    Vec3& operator-=(const Vec3& u) noexcept;

    inline
    Vec3& operator*=(const T& s) noexcept;

    inline
    Vec3& operator/=(const T& s) noexcept;

    inline
    T& operator[](int i) noexcept;

    inline
    const T& operator[](int i) const noexcept;

    inline
    int maxDir(const bool a_doAbs) const noexcept;

    inline
    bool operator==(const Vec3& u) const noexcept;

    inline
    bool operator<(const Vec3& u) const noexcept;

    inline
    bool operator>(const Vec3& u) const noexcept;

    inline
    bool operator<=(const Vec3& u) const noexcept;

    inline
    bool operator>=(const Vec3& u) const noexcept;

    inline
    Vec3 cross(const Vec3& u) const noexcept;

    inline
    T dot(const Vec3& u) const noexcept;

    inline
    T length() const noexcept;

    inline
    T length2() const noexcept;

  protected:
    T X[3];
  };


  template <class T>
  inline
  Vec3<T> operator*(const T& s, const Vec3<T>& a_other) noexcept;

  template <class T>
  inline
  Vec3<T> operator/(const T& s, const Vec3<T>& a_other) noexcept;
}

#include "dcel_vecI.H"

#endif
