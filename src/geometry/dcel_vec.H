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

  class Vec2 {
  public:

    static const Vec2 Zero;
    static const Vec2 Unit;
    
    Vec2() = default;
    Vec2(const double a_x, const double a_y);
    ~Vec2() = default;
    
    double x;
    double y;

    inline
    Vec2& operator=(const Vec2& a_other) noexcept;

    inline
    Vec2& operator+=(const Vec2& a_other) noexcept;

    inline
    Vec2& operator-=(const Vec2& a_other) noexcept;

    inline
    Vec2 operator+(const Vec2& a_other) const noexcept;

    inline
    Vec2 operator-(const Vec2& a_other) const noexcept;

    inline
    double dot(const Vec2& a_other) noexcept;

    inline
    double length() const noexcept;

    inline
    double length2() const noexcept;
  };

  class Vec3 {
  public:

    static const Vec3 Zero;
    static const Vec3 Unit;
    
    Vec3() = default;
    Vec3(const Vec3& u) = default;
    Vec3(const double a_x, const double a_y, const double a_z);
    ~Vec3() = default;

    inline
    Vec3& operator=(const Vec3& u) noexcept;

    inline
    Vec3& operator+=(const Vec3& u) noexcept;

    inline
    Vec3& operator-=(const Vec3& u) noexcept;

    inline
    Vec3 operator+(const Vec3& u) const noexcept;

    inline
    Vec3 operator-(const Vec3& u) const noexcept;

    inline
    double& operator[](int i) noexcept;

    inline
    const double& operator[](int i) const noexcept;

    inline
    Vec3 cross(const Vec3& u) const noexcept;

    inline
    double dot(const Vec3& u) noexcept;

    inline
    double length() const noexcept;

    inline
    double length2() const noexcept;

  protected:
    double X[3];
  };
}

#include "dcel_vecI.H"

#endif
