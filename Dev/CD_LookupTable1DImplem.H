/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LookupTable1D.H
  @brief  Declaration of a lookup table in one independent variable.
  @author Robert Marskar
*/

#ifndef CD_LookupTableImplem_H
#define CD_LookupTableImplem_H

// Std includes
#include <algorithm>
#include <math.h>

// Our includes
#include <CD_LookupTable1D.H>

template <typename T, size_t N, typename I>
LookupTable1D<T, N, I>::LookupTable1D() noexcept
{
  this->reset();
};

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::reset() noexcept
{
  m_isGood          = false;
  m_grid            = std::make_tuple(CoordinateSystem::Uniform, -1.0, -1.0, -1.0);
  m_rangeStrategyLo = OutOfRangeStrategy::Constant;
  m_rangeStrategyHi = OutOfRangeStrategy::Constant;

  m_rawData.clear();
  m_structuredData.clear();
}

template <typename T, size_t N, typename I>
template <typename... Ts>
inline void
LookupTable1D<T, N, I>::addData(const Ts&... x)
{
  std::array<T, sizeof...(Ts)> arr = {(T)x...};

  m_rawData.emplace_back(arr);
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::swap(const size_t a_columnOne, const size_t a_columnTwo) noexcept
{
  for (auto& r : m_rawData) {
    const auto tmp = r[a_columnOne];

    r[a_columnOne] = r[a_columnTwo];
    r[a_columnTwo] = tmp;
  }
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::truncate(const T& a_min, const T& a_max, const size_t a_column) noexcept
{
  std::vector<std::array<T, N + 1>> truncatedData;

  for (const auto& xy : m_rawData) {
    if (xy[a_column] >= a_min || xy[a_column] <= a_max) {
      truncatedData.emplace_back(xy);
    }
  }

  m_rawData(std::move(truncatedData));
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::setRangeStrategyLo(const OutOfRangeStrategy& a_strategy) noexcept
{
  m_rangeStrategyLo = a_strategy;
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::setRangeStrategyHi(const OutOfRangeStrategy& a_strategy) noexcept
{
  m_rangeStrategyHi = a_strategy;
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::prepareTable(const size_t&           a_independentVariable,
                                     const size_t&           a_numPoints,
                                     const CoordinateSystem& a_coordSys)
{

  const std::string baseError = "LookupTable1D<T,N,I>::prepareTable";

  if (a_numPoints <= 1) {
    throw std::runtime_error(baseError + " - must have 'a_numPoints > 1'");
  }

  // Sort the raw data along the specified column.
  using Values                   = std::array<T, N + 1>;
  std::vector<Values> sortedData = m_rawData;

  std::sort(sortedData.begin(),
            sortedData.end(),
            [i = a_independentVariable](const Values& X1, const Values& X2) -> bool {
              return X1[i] < X2[i];
            });

  // Set up the 1D grid
  T xmin;
  T xmax;
  T delta;

  xmin = sortedData.front()[a_independentVariable];
  xmax = sortedData.back()[a_independentVariable];

  std::vector<T> coords(a_numPoints);

  switch (a_coordSys) {
  case CoordinateSystem::Uniform: {
    delta = (xmax - xmin) / (a_numPoints - 1);

    for (size_t i = 0; i < a_numPoints; i++) {
      coords[i] = xmin + i * delta;
    }

    break;
  }
  case CoordinateSystem::Logarithmic: {
    if (xmin <= std::numeric_limits<T>::min()) {
      throw std::runtime_error(baseError + " - but can't have negativ values for logarithmic grid");
    }
    else {
      delta = log(xmax / xmin) / (a_numPoints - 1);

      for (int i = 0; i < a_numPoints; i++) {
        coords[i] = xmin * exp(i * delta);
      }

      break;
    }
  }
  default: {
    throw std::runtime_error(baseError + " - logic bust (unsupported coordinate system)");

    break;
  }
  }

  m_structuredData.clear();

  size_t curRow = 0;

  for (const auto& x : coords) {

    while (true) {
      if (curRow >= m_rawData.size()) {
        std::cout << curRow << std::endl;
        throw std::runtime_error(baseError + " - logic bust (array out of range)");
      }

      const auto& row1 = m_rawData[curRow];
      const auto& row2 = m_rawData[curRow + 1];

      const auto x1 = row1[a_independentVariable];
      const auto x2 = row2[a_independentVariable];

      if (x >= x1 && x <= x2) {
        const auto t = (x - x1) / (x2 - x1);

        std::array<T, N + 1> interpRow;

        for (size_t i = 0; i <= N; i++) {
          //          interpRow[i] = row1[i] + t * (row2[i] - row1[i]);
        }

        m_structuredData.emplace_back(interpRow);

        break;
      }
      else {
        curRow++;
      }
    }
  }

  m_grid = std::make_tuple(a_coordSys, xmin, xmax, delta);
}

#endif
