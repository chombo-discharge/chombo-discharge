/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RealProxy.H
  @brief  Lightweight stand-ins for Chombo's Real/RealVect so the SoA prototype can
          be compiled and unit-tested in isolation (no Chombo build required).
  @author Robert Marskar
  @details This is ONLY for the Dev/ prototype. In the real implementation these
           types are replaced by Chombo's <REAL.H> and <RealVect.H>. The proxy
           deliberately mimics the small subset of the RealVect API that the
           particle code touches (operator[], arithmetic, component count).
*/

#ifndef CD_REALPROXY_H
#define CD_REALPROXY_H

// Std includes
#include <array>
#include <cstddef>
#include <ostream>

namespace proxy {

  /**
    @brief Spatial dimension. Set to 2 or 3 to mirror Chombo's CH_SPACEDIM.
  */
#ifndef PROXY_SPACEDIM
#define PROXY_SPACEDIM 3
#endif

  static constexpr int SpaceDim = PROXY_SPACEDIM;

  /**
    @brief Floating-point precision stand-in for Chombo's Real.
  */
  using Real = double;

  /**
    @brief Minimal stand-in for Chombo's RealVect.
    @details Trivially copyable so it can be memcpy'd during linearization, exactly
    like the real RealVect.
  */
  class RealVect
  {
  public:
    RealVect() noexcept
    {
      m_data.fill(0.0);
    }

    explicit RealVect(const Real a_val) noexcept
    {
      m_data.fill(a_val);
    }

    Real&
    operator[](const int a_dir) noexcept
    {
      return m_data[a_dir];
    }

    const Real&
    operator[](const int a_dir) const noexcept
    {
      return m_data[a_dir];
    }

    RealVect&
    operator+=(const RealVect& a_other) noexcept
    {
      for (int d = 0; d < SpaceDim; d++) {
        m_data[d] += a_other.m_data[d];
      }
      return *this;
    }

    RealVect&
    operator*=(const Real a_scale) noexcept
    {
      for (int d = 0; d < SpaceDim; d++) {
        m_data[d] *= a_scale;
      }
      return *this;
    }

    static RealVect
    Zero() noexcept
    {
      return RealVect(0.0);
    }

  private:
    std::array<Real, SpaceDim> m_data;
  };

  inline std::ostream&
  operator<<(std::ostream& a_os, const RealVect& a_vec)
  {
    a_os << "(";
    for (int d = 0; d < SpaceDim; d++) {
      a_os << a_vec[d] << (d < SpaceDim - 1 ? ", " : "");
    }
    a_os << ")";
    return a_os;
  }

} // namespace proxy

#endif
