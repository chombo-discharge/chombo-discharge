/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LookupTable1D.H
  @brief  Declaration of a table for looking up coefficients etc. 
  @author Robert Marskar
*/

#ifndef CD_LookupTable1D_H
#define CD_LookupTable1D_H

// Std includes
#include <iostream>
#include <vector>
#include <array>
#include <tuple>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>
#include <CD_LookupTable.H>

/*!
  @brief Class for interpolation of f = f(x) data. 
*/
template <int N>
class LookupTable1D
{
public:
  static_assert(N > 1, "LookupTable1D<N> must have N > 1");
  /*!
    @brief Default constructor. Creates a table without any entries. 
  */
  LookupTable1D();

  /*!
    @brief Copy constructor. 
  */
  LookupTable1D(const LookupTable1D& a_table);

  /*!
    @brief Destructor (does nothing). 
  */
  virtual ~LookupTable1D();

  /*!
    @brief Clear all data
  */
  inline void
  clear();

  /*!
    @brief Make the table into a uniform table with the number of entries. 
    @details This will make the data uniform with constant spacing along the specified column data. 
    @note The spacing that is used is based on the smallest and largest value in the specified column. 
  */
  inline void
  makeUniform(const int a_numRows);

  /*!
    @brief Scale data in column by a specified factor. 
    @details The template parameter is the column to be scaled. 
    @param[in] a_scale       Factor to scale by
    @param[in] a_whichColumn Column to scale
  */
  template <int K>
  inline void
  scale(const Real a_scale);

  /*!
    @brief Swap two columns
    @param[in] a_column1 First column
    @param[in] a_column2 Second column
  */
  inline void
  swap(const int a_column1, const int a_column2);

  /*!
    @brief Restrict table range to maximum and minimum. Data outside these ranges are discarded. 
    @param[in] a_min         Minimum x-value in the table. 
    @param[in] a_max         Maximum x-value in the table. 
    @param[in] a_whichColumn Ranges based on this column. 
  */
  inline void
  setRange(const Real a_min, const Real a_max, const int a_whichColumn = 0);

  /*!
    @brief Restrict table range to minimum value. Data outside these ranges are discarded. 
    @param[in] a_min         Minimum x-value in the table. 
    @param[in] a_whichColumn Ranges based on this column. 
  */
  inline void
  setMinRange(const Real a_min, const int a_whichColumn = 0);

  /*!
    @brief Restrict table range to maximum value. Data outside this range are discarded. 
    @param[in] a_max         Maximum x-value in the table. 
    @param[in] a_whichColumn Ranges based on this column. 
  */
  inline void
  setMaxRange(const Real a_max, const int a_whichColumn = 0);

  /*!
    @brief Set the strategy when retrieving data that is out of range on the low side
    @details Set the input parameter a_strategy to Constant if you want to virtually pad the data with the last endpoint or Linear if you want linear interpolation.
    @param[in] a_strategy Range strategy. 
  */
  inline void
  setOutOfRangeStrategyLow(const OutOfRangeStrategy a_strategy);

  /*!
    @brief Set the strategy when retrieving data that is out of range on the high side
    @details Set the input parameter a_strategy to Constant if you want to virtually pad the data with the last endpoint or Linear if you want linear interpolation.
    @param[in] a_strategy Range strategy. 
  */
  inline void
  setOutOfRangeStrategyHigh(const OutOfRangeStrategy a_strategy);

  /*!
    @brief Set the table spacing mode
    @param[in] a_spacing Spacing mode.
  */
  inline void
  setTableSpacing(const TableSpacing a_spacing);

  /*!
    @brief Sort the table. Useful if the input data was not already sorted. 
    @param[in] a_whichColumn Sort along this column. 
    @details On output, the data in the table is sorted -- smallest x goes first. 
  */
  inline void
  sort(const int a_whichColumn = 0);

  /*!
    @brief Add entry. 
    @param[in] x Entry to add. You can add an entire array through the variadic template. 
  */
  template <typename... Ts>
  inline void
  addEntry(const Ts&... x);

  /*!
    @brief Dump table to pout()
  */
  inline void
  dumpTable(std::ostream& a_outputStream = std::cout) const;

  /*!
    @brief Dump table to a file. 
    @param[in] a_fileName File name
  */
  inline void
  dumpTable(const std::string a_fileName) const;

  /*!
    @brief Get entry
    @note When calling this, make sure the data is sorted AND uniform. For performance reasons this is not checked when calling this routine. 
    @param[in] a_x Interpolation point (independent variable)
  */
  template <int K>
  inline Real
  getEntry(const Real a_x) const;

  /*!
    @brief Get data.
    @note When calling this, make sure the data is sorted AND uniform. For performance reasons this is not checked when calling this routine. 
    @param[in] a_x Interpolation point (independent variable)
  */
  inline std::array<Real, N>
  getData(const Real a_x) const;

  /*!
    @brief Get size of the table.
    @returns m_data.size()
  */
  inline int
  getNumEntries() const;

protected:
  /*!
    @brief End point stategy on the low end of the data 
  */
  OutOfRangeStrategy m_strategyLo;

  /*!
    @brief End point stategy on the high end of the data 
  */
  OutOfRangeStrategy m_strategyHi;

  /*!
    @brief For setting the table spacing
  */
  TableSpacing m_spacing;

  /*!
    @brief Table state. First entry is the column along which we sort and make the table uniform. The second entry is 
    whether or not the table is actually sorted, and the third entry is whether or not the table is uniform along the sorted column.
  */
  std::tuple<int, bool, bool> m_state;

  /*!
    @brief A base resolution for the electric field tabulations (we assumed that the rates were mapped out linearly)
  */
  Real m_delta;

  /*!
    @brief Data
  */
  std::vector<std::array<Real, N>> m_data;

  /*!
    @brief Direct lookup, assumes non-uniform ordering. 
    @details This is very slow, because it will look through the entire table. 
    @note The data MUST be sorted along the input column for this to work. 
  */
  inline std::array<Real, N>
  directLookup(const Real a_x) const;

  /*!
    @brief Get the lower index during a lookup
    @param[in] a_x Interpolation point. 
  */
  inline int
  getIndexLo(const Real a_x) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_LookupTable1DImplem.H>

#endif
