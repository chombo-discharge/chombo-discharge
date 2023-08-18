/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LookupTable1D.H
  @brief  Declaration of a lookup table in one independent variable.
  @author Robert Marskar
*/

#ifndef CD_LookupTable1D_H
#define CD_LookupTable1D_H

// Std includes
#include <iostream>
#include <vector>
#include <array>
#include <type_traits>

// Our includes
#include <CD_LookupTable.H>

/*!
  @brief Class for interpolation of f = f(x) data in one independent variable x.
  @details N is the number of additional dependent variables. T must be a floating point type (e.g., float, double)
*/
template <typename T = Real, size_t N = 1, typename I = std::enable_if_t<std::is_floating_point<T>::value>>
class LookupTable1D
{
public:
  static_assert(N >= 1, "LookupTable1D<N> must have N >= 1");

  /*!
    @brief Default constructor. Creates a table without any entries. 
  */
  LookupTable1D() noexcept;

  /*!
    @brief Destructor (does nothing). 
  */
  virtual ~LookupTable1D() noexcept = default;

  /*!
    @brief Reset everything.
  */
  inline void
  reset() noexcept;

  /*!
    @brief Slice this table, keeping only the user-specified columns.
    @details This returns a new table with a reduced number of columns.
    @param[in] 
  */
  template <size_t M, typename... Ts>
  LookupTable1D<T, M, I>
  slice(const Ts&... a_columns) const noexcept;

  /*!
    @brief Add entry. 
    @param[in] x Entry to add. For example addData(1,1,1,1). Number of elements in x must be N+1
  */
  template <typename... Ts>
  inline void
  addData(const Ts&... x) noexcept;

  /*!
    @brief Add entry. 
    @param[in] x Entry to add. 
  */
  inline void
  addData(const std::array<T, N + 1>& x) noexcept;

  /*!
    @brief Utility function for swapping columns. 
    @details This is done on the raw data -- if the user wants to swap columns in the resampled/structured data
    then he needs to call this function first and then resample the table through prepareTable
    @param[in] a_columnOne Column to swap.
    @param[in] a_columnTwo Column to swap.
  */
  inline void
  swap(const size_t a_columnOne, const size_t a_columnTwo) noexcept;

  /*!
    @brief Utility function which scales one of the columns (either dependent or independent variable)
  */
  template <size_t K>
  inline void
  scale(const T& a_scale) noexcept;

  /*!
    @brief Utility function for truncating raw data along one of the variables (either dependent or independent).
    @details This will discard (from the raw data) all data that fall outside the input interval. This is done on the raw data -- the user will
    need to call prepareTable if the result should propagate into the resampled/structured data.
    @param[in] a_min Minimum value represented.
    @param[in] a_max Maximum value represented.
  */
  inline void
  truncate(const T& a_min, const T& a_max, const size_t a_column) noexcept;

  /*!
    @brief Set the out-of-range strategy on the low end
    @param[in] a_strategy Out-of-range strategy on low end
  */
  inline void
  setRangeStrategyLo(const LookupTable::OutOfRangeStrategy& a_strategy) noexcept;

  /*!
    @brief Set the out-of-range strategy on the high end
    @param[in] a_strategy Out-of-range strategy on the high end
  */
  inline void
  setRangeStrategyHi(const LookupTable::OutOfRangeStrategy& a_strategy) noexcept;

  /*!
    @brief Turn the raw data into uniform data for fast lookup.
    @param[in] a_independentVariable The independent variable (i.e., column in the input data).
    @param[in] a_numPoints Number of points in the input table. 
    @param[in] a_spacing Table spacing
  */
  inline void
  prepareTable(const size_t& a_independentVariable, const size_t& a_numPoints, const LookupTable::Spacing& a_spacing);

  /*!
    @brief Interpolation function for specific dependent variable K
    @param[in] a_x Independent variable x
  */
  template <size_t K>
  inline T
  interpolate(const T& a_x) const;

  /*!
    @brief Interpolate whole table.
    @param[in] a_x Independent variable x
  */
  inline std::array<T, N + 1>
  interpolate(const T& x) const;

  /*!
    @brief Access function for raw data. 
    @return Returns m_rawData
  */
  inline std::vector<std::array<T, N + 1>>&
  getRawData() noexcept;

  /*!
    @brief Access function for structured data. 
    @return Returns m_structuredData
  */
  inline std::vector<std::array<T, N + 1>>&
  getStructuredData() noexcept;

  /*!
    @brief Access function for raw data.
    @return Returns m_rawData
  */
  inline const std::vector<std::array<T, N + 1>>&
  getRawData() const noexcept;

  /*!
    @brief Access function for structured data.
    @return Returns m_structuredData
  */
  inline const std::vector<std::array<T, N + 1>>&
  getStructuredData() const noexcept;

  /*!
    @brief Dump raw table data to file
    @param[in] a_file File name
  */
  inline void
  writeRawData(const std::string& a_file) const noexcept;

  /*!
    @brief Dump structured table data to file
    @param[in] a_file File name
  */
  inline void
  writeStructuredData(const std::string& a_file) const noexcept;

  /*!
    @brief Dump raw table data to output stream.
    @param[in] a_ostream Output stream
  */
  inline void
  outputRawData(std::ostream& a_ostream = std::cout) const noexcept;

  /*!
    @brief Dump structured table data to file.
    @param[in] a_ostream Output stream
  */
  inline void
  outputStructuredData(std::ostream& a_ostream = std::cout) const noexcept;

protected:
  /*!
    @brief Check if data can be interpolated
  */
  bool m_isGood;

  /*!
    @brief Out-of-range strategy on low end
  */
  LookupTable::OutOfRangeStrategy m_rangeStrategyLo;

  /*!
    @brief Out-of-range strategy on high end
  */
  LookupTable::OutOfRangeStrategy m_rangeStrategyHi;

  /*!
    @brief Underlying 1D grid. This is populated when calling prepareTable.
  */
  std::tuple<LookupTable::Spacing, size_t, T, T, T> m_grid;

  /*!
    @brief Raw data
  */
  std::vector<std::array<T, N + 1>> m_rawData;

  /*!
    @brief Structured data. This is populated when calling prepareTable.
  */
  std::vector<std::array<T, N + 1>> m_structuredData;

  /*!
    @brief Utility function for outputting data.
    @param[in] a_ostream Output stream
    @param[in] a_data Data to be sent to output stream
  */
  inline void
  outputData(std::ostream& a_ostream, const std::vector<std::array<T, N + 1>>& a_data) const noexcept;

  /*!
    @brief Utility function for outputting data to a file
    @param[in] a_file File name
    @param[in] a_data Data to be sent to output stream
  */
  inline void
  writeToFile(const std::string& a_file, const std::vector<std::array<T, N + 1>>& a_data) const noexcept;

  /*!
    @brief Get the lower index that brackets the input variable between two data points in the structured grid.
    @param[in] a_x Independent variable.
  */
  inline size_t
  getIndexLo(const T& a_x) const;
};

#include <CD_LookupTable1DImplem.H>

#endif
