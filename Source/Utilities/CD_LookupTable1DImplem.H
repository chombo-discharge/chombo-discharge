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
#include <iomanip>
#include <fstream>
#include <limits>

// Our includes
#include <CD_LookupTable1D.H>

template <typename T, size_t N, typename I>
LookupTable1D<T, N, I>::LookupTable1D() noexcept
{
  this->reset();
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::reset() noexcept
{
  m_isGood          = false;
  m_grid            = std::make_tuple(LookupTable::Spacing::Uniform, 0, -1.0, -1.0, -1.0);
  m_rangeStrategyLo = LookupTable::OutOfRangeStrategy::Constant;
  m_rangeStrategyHi = LookupTable::OutOfRangeStrategy::Constant;

  m_rawData.clear();
  m_structuredData.clear();
}

template <typename T, size_t N, typename I>
template <typename... Ts>
inline void
LookupTable1D<T, N, I>::addData(const Ts&... x) noexcept
{
  std::array<T, sizeof...(Ts)> arr = {(T)x...};

  m_rawData.emplace_back(arr);
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::addData(const std::array<T, N + 1>& x) noexcept
{
  m_rawData.emplace_back(x);
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

  for (auto& r : m_structuredData) {
    const auto tmp = r[a_columnOne];

    r[a_columnOne] = r[a_columnTwo];
    r[a_columnTwo] = tmp;
  }

  if (m_isGood) {
    if (a_columnOne == std::get<4>(m_grid)) {
      std::get<4>(m_grid) = a_columnOne;
    }
    else if (a_columnTwo == std::get<4>(m_grid)) {
      std::get<4>(m_grid) = a_columnTwo;
    }
  }
}

template <typename T, size_t N, typename I>
template <size_t K>
inline void
LookupTable1D<T, N, I>::scale(const T& a_scale) noexcept
{
  for (auto& r : m_rawData) {
    r[K] *= a_scale;
  }

  for (auto& r : m_structuredData) {
    r[K] *= a_scale;
  }

  if (m_isGood && K == std::get<2>(m_grid)) {
    std::get<4>(m_grid) *= a_scale;
  }
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::truncate(const T& a_min, const T& a_max, const size_t a_column) noexcept
{
  std::vector<std::array<T, N + 1>> truncatedData;

  for (const auto& xy : m_rawData) {
    if (xy[a_column] >= a_min && xy[a_column] <= a_max) {
      truncatedData.emplace_back(xy);
    }
  }

  m_rawData = truncatedData;
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::setRangeStrategyLo(const LookupTable::OutOfRangeStrategy& a_strategy) noexcept
{
  m_rangeStrategyLo = a_strategy;
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::setRangeStrategyHi(const LookupTable::OutOfRangeStrategy& a_strategy) noexcept
{
  m_rangeStrategyHi = a_strategy;
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::prepareTable(const size_t&               a_independentVariable,
                                     const size_t&               a_numPoints,
                                     const LookupTable::Spacing& a_spacing)
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

  switch (a_spacing) {
  case LookupTable::Spacing::Uniform: {
    delta = (xmax - xmin) / (a_numPoints - 1);

    for (size_t i = 0; i < a_numPoints; i++) {
      coords[i] = xmin + i * delta;
    }

    break;
  }
  case LookupTable::Spacing::Exponential: {
    if (xmin <= std::numeric_limits<T>::min()) {
      throw std::runtime_error(baseError + " - but must have all 'x > 0.0' for logarithmic grid");
    }

    else {
      delta = log10(xmax / xmin) / (a_numPoints - 1);

      for (int i = 0; i < a_numPoints; i++) {
        coords[i] = xmin * std::pow(10.0, i * delta);
      }

      break;
    }
  }
  default: {
    throw std::runtime_error(baseError + " - logic bust (unsupported table spacing system)");

    break;
  }
  }

  m_structuredData.clear();

  size_t curRow = 0;

  for (const auto& x : coords) {
    while (curRow < m_rawData.size() - 1) {
      const auto& row1 = m_rawData[curRow];
      const auto& row2 = m_rawData[curRow + 1];

      const auto x1 = row1[a_independentVariable];
      const auto x2 = row2[a_independentVariable];

      if (x >= x1 && x <= x2) {
        const auto t = (x - x1) / (x2 - x1);

        std::array<T, N + 1> interpRow;

        for (size_t i = 0; i <= N; i++) {
          interpRow[i] = row1[i] + t * (row2[i] - row1[i]);
        }

        m_structuredData.emplace_back(interpRow);

        break;
      }

      curRow++;
    }
  }

  m_grid = std::make_tuple(a_spacing, a_independentVariable, xmin, xmax, delta);

  m_isGood = true;
}

template <typename T, size_t N, typename I>
inline std::vector<std::array<T, N + 1>>&
LookupTable1D<T, N, I>::getRawData() noexcept
{
  return (m_rawData);
}

template <typename T, size_t N, typename I>
inline std::vector<std::array<T, N + 1>>&
LookupTable1D<T, N, I>::getStructuredData() noexcept
{
  return (m_structuredData);
}

template <typename T, size_t N, typename I>
inline const std::vector<std::array<T, N + 1>>&
LookupTable1D<T, N, I>::getRawData() const noexcept
{
  return (m_rawData);
}

template <typename T, size_t N, typename I>
inline const std::vector<std::array<T, N + 1>>&
LookupTable1D<T, N, I>::getStructuredData() const noexcept
{
  return (m_structuredData);
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::writeRawData(const std::string& a_file) const noexcept
{
  this->writeToFile(a_file, m_rawData);
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::writeStructuredData(const std::string& a_file) const noexcept
{
  this->writeToFile(a_file, m_structuredData);
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::outputRawData(std::ostream& a_ostream) const noexcept
{
  this->outputData(a_ostream, m_rawData);
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::outputStructuredData(std::ostream& a_ostream) const noexcept
{
  this->outputData(a_ostream, m_structuredData);
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::outputData(std::ostream&                            a_ostream,
                                   const std::vector<std::array<T, N + 1>>& a_data) const noexcept
{
  for (const auto& r : a_data) {
    for (const auto& c : r) {
      a_ostream << std::left << std::setw(14) << c;
    }
    a_ostream << "\n";
  }
}

template <typename T, size_t N, typename I>
inline void
LookupTable1D<T, N, I>::writeToFile(const std::string&                       a_file,
                                    const std::vector<std::array<T, N + 1>>& a_data) const noexcept
{

#ifdef CH_MPI
  if (procID() == 0) {
#endif
    std::ofstream file;

    file.open(a_file);
    this->outputData(file, a_data);
    file.close();
#ifdef CH_MPI
  }
#endif
}

template <typename T, size_t N, typename I>
inline size_t
LookupTable1D<T, N, I>::getIndexLo(const T& a_x) const
{
  const LookupTable::Spacing& spacing = std::get<0>(m_grid);
  const T&                    x0      = std::get<2>(m_grid);
  const T&                    delta   = std::get<4>(m_grid);

  size_t idxLo;

  switch (spacing) {
  case LookupTable::Spacing::Uniform: {
    idxLo = std::floor((a_x - x0) / delta);

    break;
  }
  case LookupTable::Spacing::Exponential: {
    idxLo = std::floor(log10(a_x / x0) / delta);

    break;
  }
  default: {
    throw std::runtime_error("LookupTable1D<T,N,I>::getIndexLo - logic bust (unsupported table spacing system)");
  }
  }

  return idxLo;
}

template <typename T, size_t N, typename I>
inline std::array<T, N + 1>
LookupTable1D<T, N, I>::interpolate(const T& a_x) const
{
  if (!m_isGood) {
    throw std::runtime_error("LookupTable1D<T, N, I>::interpolate(array) but need to call 'prepareTable first'");
  }

  const size_t& indVar = std::get<1>(m_grid);
  const T&      xmin   = std::get<2>(m_grid);
  const T&      xmax   = std::get<3>(m_grid);

  std::array<T, N + 1> ret;

  if (a_x < xmin) {
    switch (m_rangeStrategyLo) {
    case LookupTable::OutOfRangeStrategy::Constant: {
      ret = m_structuredData.front();

      break;
    }
    case LookupTable::OutOfRangeStrategy::Interpolate: {
      const T diffDx   = m_structuredData[1][indVar] - m_structuredData[0][indVar];
      const T extrapDx = m_structuredData[0][indVar] - a_x;

      ret = m_structuredData.front();

      for (size_t i = 0; i < N + 1; i++) {
        const T slope = (m_structuredData[1][i] - m_structuredData[0][i]) / diffDx;

        ret[i] -= slope * extrapDx;
      }

      break;
    }
    default: {
      throw std::runtime_error("LookupTable1D<T, N, I>::interpolate(array) unsupported range strategy at low end");

      break;
    }
    }
  }
  else if (a_x > xmax) {
    switch (m_rangeStrategyHi) {
    case LookupTable::OutOfRangeStrategy::Constant: {
      ret = m_structuredData.back();

      break;
    }
    case LookupTable::OutOfRangeStrategy::Interpolate: {
      const size_t& size = m_structuredData.size();

      const T diffDx   = m_structuredData[size - 1][indVar] - m_structuredData[size - 2][indVar];
      const T extrapDx = a_x - m_structuredData[size - 1][indVar];

      ret = m_structuredData.back();

      for (size_t i = 0; i < N + 1; i++) {
        const T slope = (m_structuredData[size - 1][i] - m_structuredData[size - 2][i]) / diffDx;

        ret[i] += slope * extrapDx;
      }

      break;
    }
    default: {
      throw std::runtime_error("LookupTable1D<T, N, I>::interpolate(array) unsupported range strategy at low end");

      break;
    }
    }
  }
  else {
    const size_t idx = this->getIndexLo(a_x);

    const T diffDx   = m_structuredData[idx + 1][indVar] - m_structuredData[idx][indVar];
    const T interpDx = a_x - m_structuredData[idx][indVar];

    ret = m_structuredData[idx];

    for (size_t i = 0; i < N + 1; i++) {
      const T slope = (m_structuredData[idx + 1][i] - m_structuredData[idx][i]) / diffDx;

      ret[i] += slope * interpDx;
    }
  }

  return ret;
}

template <typename T, size_t N, typename I>
template <size_t K>
inline T
LookupTable1D<T, N, I>::interpolate(const T& a_x) const
{
  return std::get<K>(this->interpolate(a_x));
}

template <typename T, size_t N, typename I>
template <size_t M, typename... Ts>
LookupTable1D<T, M, I>
LookupTable1D<T, N, I>::slice(const Ts&... a_columns) const noexcept
{
  static_assert(M <= N, "LookupTable1D<T, N,I>::slice must have have 'M <= N'");

  std::array<size_t, sizeof...(Ts)> arr = {(size_t)a_columns...};

  LookupTable1D<T, M, I> ret;

  for (const auto& r : m_rawData) {
    std::array<T, M + 1> slicedData;

    for (size_t i = 0; i < M + 1; i++) {
      slicedData[i] = r[arr[i]];
    }

    ret.addData(slicedData);
  }

  return ret;
}

#endif
