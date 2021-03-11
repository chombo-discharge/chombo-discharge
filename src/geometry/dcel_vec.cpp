/*!
  @file   dcel_vec.cpp
  @brief  Implementation of dcel_vec.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_vec.H"

using namespace dcel;

template class Vec2<float>;
template class Vec2<double>;

template <>
const Vec2<float> Vec2<float>::Zero(0., 0.);

template <>
const Vec2<float> Vec2<float>::Unit(1., 1.);

template <>
const Vec2<double> Vec2<double>::Zero(0., 0.);

template <>
const Vec2<double> Vec2<double>::Unit(1., 1.);


template class Vec3<float>;
template class Vec3<double>;

template <>
const Vec3<float> Vec3<float>::Zero(0., 0., 0.);

template <>
const Vec3<float> Vec3<float>::Unit(1., 1., 1.);

template <>
const Vec3<double> Vec3<double>::Zero(0., 0., 0.);

template <>
const Vec3<double> Vec3<double>::Unit(1., 1., 1.);
