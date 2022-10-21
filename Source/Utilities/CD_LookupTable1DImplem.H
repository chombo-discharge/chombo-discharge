/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LookupTable1DImplem.H
  @brief  Implementation of CD_LookupTable1D.H
  @author Robert Marskar
*/

#ifndef CD_LookupTable1DImplem_H
#define CD_LookupTable1DImplem_H

// Std includes
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>

// Our includes
#include <CD_LookupTable1D.H>
#include <CD_NamespaceHeader.H>

template <int N>
inline LookupTable1D<N>::LookupTable1D()
{

  // Default settings
  m_spacing = TableSpacing::Uniform;
  m_state   = std::make_tuple(-1, false, false);

  m_strategyLo = OutOfRangeStrategy::Constant;
  m_strategyHi = OutOfRangeStrategy::Constant;
}

template <int N>
inline LookupTable1D<N>::~LookupTable1D()
{}

template <int N>
inline LookupTable1D<N>::LookupTable1D(const LookupTable1D& a_table)
{
  m_data    = a_table.m_data;
  m_state   = a_table.m_state;
  m_delta   = a_table.m_delta;
  m_spacing = a_table.m_spacing;

  m_strategyLo = a_table.m_strategyLo;
  m_strategyHi = a_table.m_strategyHi;
}

template <int N>
template <typename... Ts>
inline void
LookupTable1D<N>::addEntry(const Ts&... x)
{
  std::array<Real, sizeof...(Ts)> arr = {x...};
  m_data.emplace_back(arr);

  std::get<Sorted>(m_state)  = false;
  std::get<Uniform>(m_state) = false;
}

template <int N>
inline void
LookupTable1D<N>::sort(const int a_sortedColumn)
{

  auto comparator = [idx = a_sortedColumn](const std::array<Real, N>& a, const std::array<Real, N>& b) -> bool {
    return a[idx] < b[idx];
  };

  std::sort(m_data.begin(), m_data.end(), comparator);

  std::get<Column>(m_state) = a_sortedColumn;
  std::get<Sorted>(m_state) = true;
}

template <int N>
template <int K>
inline void
LookupTable1D<N>::scale(const Real a_scale)
{
  for (auto& r : m_data) {
    r[K] *= a_scale;
  }

  if (K == std::get<Column>(m_state)) {
    m_delta *= a_scale;
  }
}

template <int N>
inline void
LookupTable1D<N>::swap(const int a_column1, const int a_column2)
{
  for (auto& r : m_data) {
    const Real tmp = r[a_column1];

    r[a_column1] = r[a_column2];
    r[a_column2] = tmp;
  }

  if (std::get<Sorted>(m_state)) {

    if (a_column1 == std::get<Column>(m_state)) {
      std::get<Column>(m_state) = a_column1;
    }
    else if (a_column2 == std::get<Column>(m_state)) {
      std::get<Column>(m_state) = a_column2;
    }
  }
}

template <int N>
inline void
LookupTable1D<N>::clear()
{
  m_data.resize(0);
}

template <int N>
inline void
LookupTable1D<N>::setRange(const Real a_min, const Real a_max, const int a_colRange)
{

  // Do a backup of the data and reset data holder.
  const std::vector<std::array<Real, N>> oldData = m_data;
  this->clear();

  for (const auto& d : oldData) {
    if (d[a_colRange] >= a_min && d[a_colRange] <= a_max) {
      m_data.emplace_back(d);
    }
  }

  std::get<Sorted>(m_state)  = false;
  std::get<Uniform>(m_state) = false;
}

template <int N>
inline void
LookupTable1D<N>::setOutOfRangeStrategyLow(const OutOfRangeStrategy a_strategy)
{
  m_strategyLo = a_strategy;
}

template <int N>
inline void
LookupTable1D<N>::setOutOfRangeStrategyHigh(const OutOfRangeStrategy a_strategy)
{
  m_strategyHi = a_strategy;
}

template <int N>
inline void
LookupTable1D<N>::setTableSpacing(const TableSpacing a_spacing)
{

  // If we change the table spacing, we are no longer in "uniform" mode.
  if (a_spacing != m_spacing) {
    std::get<Uniform>(m_state) = false;
  }

  m_spacing = a_spacing;
}

template <int N>
inline void
LookupTable1D<N>::dumpTable(std::ostream& a_outputStream) const
{
  for (const auto& r : m_data) {
    for (const auto& c : r) {
      a_outputStream << c << "\t";
    }
    a_outputStream << "\n";
  }
}

template <int N>
inline void
LookupTable1D<N>::dumpTable(const std::string a_fileName) const
{
#ifdef CH_MPI
  if (procID() == 0) {
#endif
    ofstream outputFile;

    outputFile.open(a_fileName);

    // Run through the data and write our stuff.
    for (const auto& row : m_data) {
      for (const auto& entry : row) {
        outputFile << entry << "\t\t";
      }
      outputFile << "\n";
    }

    outputFile.close();
#ifdef CH_MPI
  }
#endif
}

template <int N>
inline int
LookupTable1D<N>::getIndexLo(const Real a_x) const
{
  int index = 0;

  // This is the column where sorted along.
  const int sortedColumn = std::get<Column>(m_state);

  // This is the lowest point in the table.
  const Real x0 = m_data.front()[sortedColumn];

  switch (m_spacing) {
  case TableSpacing::Uniform: {
    index = std::floor((a_x - x0) / m_delta);

    break;
  }
  case TableSpacing::Exponential: {
    index = std::floor(log(a_x / x0) / m_delta);

    break;
  }
  default: {
    MayDay::Error("LookupTable1D<N>::getIndexLo - logic bust");
  }
  }

  return index;
}

template <int N>
template <int K>
inline Real
LookupTable1D<N>::getEntry(const Real a_x) const
{

  // TLDR: This is the main routine for fetching a column of data. This will switch between constant/linear interpolation outside of
  //       the stored data range.

  Real ret = 0.0;

  if (m_data.size() == 1) {
    ret = m_data.front()[K];
  }
  else {

    // We have sorted our data from smallest to largest along some column -- this is that column.
    const int sortedColumn = std::get<Column>(m_state);

    // This is the LOWER index among the two indices we will interpolate between.
    const int idx = this->getIndexLo(a_x);

    // See how we should enforce fetching the data when retrieving data outside the stored data range.
    if (idx < 0) { // Input point is outside the lower range

      // Figure out how we should interpolate on the low side.
      switch (m_strategyLo) {
      case OutOfRangeStrategy::Constant: {
        ret = m_data.front()[K];

        break;
      }
      case OutOfRangeStrategy::Linear: {
        const int size = m_data.size();

        const Real idxLo = 0;
        const Real idxHi = 1;

        const Real xLo = m_data[idxLo][sortedColumn];
        const Real xHi = m_data[idxHi][sortedColumn];

        // Figure out distance between the two points
        const Real deltaX = xHi - xLo;

        // Figure out the distance between the low index and our own input point.
        const Real dist = xLo - a_x;

        // Compute the slope.
        const Real slope = (m_data[idxHi][K] - m_data[idxLo][K]) / deltaX;

        // Do the interpolation.
        ret = m_data.front()[K] - dist * slope;

        break;
      }
      default:
        MayDay::Error("LookupTable1D<N>::getEntry -- logic bust when fetching out-of-range endpoint on the low side");

        break;
      }
    }
    else if (idx >= m_data.size() - 1) { // Input point is outside the range on the high end.

      // Figure out how to get the data outside the range.
      switch (m_strategyHi) {
      case OutOfRangeStrategy::Constant: {
        ret = m_data.back()[K];

        break;
      }
      case OutOfRangeStrategy::Linear: {
        const int size = m_data.size();

        const int idxLo = m_data.size() - 2;
        const int idxHi = m_data.size() - 1;

        const Real xLo = m_data[idxLo][sortedColumn];
        const Real xHi = m_data[idxHi][sortedColumn];

        // Compute the spacing between the two points
        const Real deltaX = xHi - xLo;

        // Compute the slope.
        const Real slope = (m_data[size - 1][K] - m_data[size - 2][K]) / deltaX;

        // Compute the distance to the endpoint.
        const Real dist = a_x - xHi;

        // Linearly interpolate.
        ret = m_data.back()[K] + dist * slope;

        break;
      }
      default:
        MayDay::Error("LookupTable1D<N>::getEntry -- logic bust when fetching out-of-range endpoint on the high side");

        break;
      }
    }
    else {

      // Compute spacing between the two points.
      const Real d = (m_data[idx + 1][sortedColumn] - m_data[idx][sortedColumn]);

      // Compute the slope between the two points.
      const Real t = (a_x - m_data[idx][sortedColumn]) / d;

      // Linearly interpolate.
      ret = m_data[idx][K] + t * (m_data[idx + 1][K] - m_data[idx][K]);
    }
  }

  return ret;
}

template <int N>
inline std::array<Real, N>
LookupTable1D<N>::getData(const Real a_x) const
{

  std::array<Real, N> ret;

  if (m_data.size() == 1) {
    ret = m_data.front();
  }
  else {
    const int sortedColumn = std::get<Column>(m_state);

    // This is the LOWER index among the two indices we will interpolate between.
    const int idx = this->getIndexLo(a_x);

    // See how we should enforce fetching the data when retrieving data outside the stored data range.
    if (idx < 0) { // On the low side.

      // Figure out how we should interpolate on the low side.
      switch (m_strategyLo) {
      case OutOfRangeStrategy::Constant: {
        ret = m_data.front();

        break;
      }
      case OutOfRangeStrategy::Linear: {
        ret = m_data.front();

        // Distance from first point in the table to the interpolation point. This
        // will be a positive number.
        const Real dist = m_data.front()[sortedColumn] - a_x;

        // Interpolate each dependent and independent variable.
        for (int i = 0; i < N; i++) {

          // This is the slope for variable i.
          const Real slope = (m_data[1][i] - m_data[0][i]) / m_delta;

          // Extrapolate beyond endpoint
          ret[i] -= dist * slope;
        }

        break;
      }
      default:
        MayDay::Error("LookupTable1D<N>::getData -- logic bust when fetching out-of-range endpoint on the low side");

        break;
      }
    }
    else if (idx >= m_data.size() - 1) {

      // Figure out how to get the data outside the range.
      switch (m_strategyHi) {
      case OutOfRangeStrategy::Constant: {
        ret = m_data.back();

        break;
      }
      case OutOfRangeStrategy::Linear: {
        const int size = m_data.size();

        ret = m_data.back();

        // This is the distance from the interpolation to the endpoint in the table. This will be a positive number.
        const Real dist = a_x - m_data.back()[sortedColumn];

        // Interpolate each dependent and independent variable.
        for (int i = 0; i < N; i++) {

          // This is the slope for variable i
          const Real slope = (m_data[size - 1][i] - m_data[size - 2][i]) / m_delta;

          ret[i] += dist * slope;
        }

        break;
      }
      default:
        MayDay::Error("LookupTable1D<N>::getEntry -- logic bust when fetching out-of-range endpoint on the high side");

        break;
      }
    }
    else {
      MayDay::Error("LookupTable1D<N>::getEntry (array version -- spacing is wrong");
      const Real t = (a_x - m_data[idx][sortedColumn]) / m_delta;

      // Interpolate all variables.
      for (int i = 0; i < N; i++) {
        ret[i] = m_data[idx][i] + t * (m_data[idx + 1][i] - m_data[idx][i]);
      }
    }
  }

  return ret;
}

template <int N>
inline std::array<Real, N>
LookupTable1D<N>::directLookup(const Real a_x) const
{

  std::array<Real, N> ret;
  if (!std::get<Sorted>(m_state)) {
    std::cerr << "LookupTable1D<N>::directLookup -- logic bust. Data must be sorted along the input column!\n";
  }
  else {
    const int column = std::get<Column>(m_state);

    if (a_x <= m_data.front()[column]) {
      ret = m_data.front();
    }
    else if (a_x >= m_data.back()[column]) {
      ret = m_data.back();
    }
    else { // Do linear interpolation. We look through the rows and find two entries that bracket a_x. We then linearly interpolate between those entries.
      for (int irow = 0; irow < m_data.size() - 1; irow++) {
        const std::array<Real, N>& row1 = m_data[irow];
        const std::array<Real, N>& row2 = m_data[irow + 1];

        const Real x1 = row1[column];
        const Real x2 = row2[column];

        if (a_x >= x1 && a_x <= x2) {
          const Real t = (a_x - x1) / (x2 - x1);

          // Now interpolate data in each column.
          for (int j = 0; j < N; j++) {
            ret[j] = row1[j] + t * (row2[j] - row1[j]);
          }
        }
      }
    }
  }

  return ret;
}

template <int N>
inline void
LookupTable1D<N>::makeUniform(const int a_numRows)
{

  if (std::get<Sorted>(m_state)) {
    const int sortedColumn = std::get<Column>(m_state);

    if (a_numRows > 1) {

      // Find minimum and maximum values for the independent variable
      const Real xMin = m_data.front()[sortedColumn];
      const Real xMax = m_data.back()[sortedColumn];

      // Create a new vector of x-values
      std::vector<Real> xValues(a_numRows, 0.0);

      switch (m_spacing) {
      case TableSpacing::Uniform: {
        // Figure out the spacing
        m_delta = (xMax - xMin) / (a_numRows - 1);

        // Create uniformly spaced points
        for (int i = 0; i < a_numRows; i++) {
          xValues[i] = xMin + i * m_delta;
        }

        break;
      }
      case TableSpacing::Exponential: {
        // Figure out the spacing
        m_delta = log(xMax / xMin) / a_numRows;

        // Create exponentially spaced points
        for (int i = 0; i < a_numRows; i++) {
          xValues[i] = xMin * exp(i * m_delta);
        }

        break;
      }
      default: {
        MayDay::Error("LookupTable1D<N>::makeUniform - logic bust");
      }
      }

      // Now linearly interpolate the data.
      std::vector<std::array<Real, N>> newData;

      for (const auto& x : xValues) {
        newData.emplace_back(this->directLookup(x));
      }

      // Now replace the data with the freshly interpolated data.
      m_data = newData;

      // Set the state flag for uniform spacing to true.
      std::get<Uniform>(m_state) = true;
    }
    else {
      std::cerr << "LookupTable1D<N>::makeUniform -- must have at least two rows!\n";
    }
  }
  else {
    std::cerr << "LookupTable1D<N>::makeUniform -- must sort the column along one of the coordinates first!\n";
  }
}

template <int N>
inline int
LookupTable1D<N>::getNumEntries() const
{
  return m_data.size();
}

#include <CD_NamespaceFooter.H>

#endif
