/*!
  @file   dcel_vec.cpp
  @brief  Implementation of dcel_vec.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_vec.H"

using namespace dcel;

template class Vec2T<float>;
template class Vec2T<double>;

template <>
const Vec2T<float> Vec2T<float>::Zero(0., 0.);

template <>
const Vec2T<float> Vec2T<float>::Unit(1., 1.);

template <>
const Vec2T<double> Vec2T<double>::Zero(0., 0.);

template <>
const Vec2T<double> Vec2T<double>::Unit(1., 1.);


template class Vec3T<float>;
template class Vec3T<double>;

template <>
const Vec3T<float> Vec3T<float>::Zero(0., 0., 0.);

template <>
const Vec3T<float> Vec3T<float>::Unit(1., 1., 1.);

template <>
const Vec3T<double> Vec3T<double>::Zero(0., 0., 0.);

template <>
const Vec3T<double> Vec3T<double>::Unit(1., 1., 1.);
