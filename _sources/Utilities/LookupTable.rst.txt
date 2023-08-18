.. _Chap:LookupTable:

Lookup tables
*************

LookupTable1D
-------------

``LookupTable1D`` is a class for interpolating data stored in a row-column format.
It is used in order to easily retrieve input data that can be stored in table formats.

.. important::

   LookupTable1D is used for data lookup *in one independent variable*.
   It does not support higher-dimensional data interpolation. 

The class is templated as

.. code-block::

   template <typename T = Real, size_t N = 1>
   class LookupTable1D

where the template parameter ``N`` indicates the number of dependent variables (``N=1`` yields a compile-time error).
Internally, the data is stored as an ``std::vector<std::array<T, N + 1> >`` where the vector entries are rows and the ``std::array<T, N + 1>`` are the column entries in that row.

.. tip::

   The internal floating point representation defaults to ``T`` and the number of dependent variables to 1. 

For performance reasons, the ``LookupTable1D`` class is used on regularly spaced data (either uniformly or logarithmically spaced).
The user can still pass irregularly spaced data into ``LookupTable1D``, and regularize the data later (i.e., interpolate it onto a regularly spaced 1D grid).
Usage of ``LookupTable1D`` will therefore consist of the following:

#. Add data rows into the table.
#. Swap columns or scale data if necessary.
#. Truncate data ranges if necessary. 
#. Regularize the table with a specified number of grid points.
#. *Retrieve data*.

These steps are discussed below.

Inserting data
______________

To add data to the table, use the member function

.. code-block:: c++

   template <typename... Ts>
   inline
   void addData(const Ts&... x);

where the parameter pack must have ``N+1`` entries.
For example, to add two rows of data to a table with two dependent variables:

.. code-block:: c++

   LookupTable1D<Real, 2> myTable;

   myTable.addData(4.0, 5.0, 6.0);
   myTable.addData(1.0, 2.0, 3.0);   

This will insert two new rows at the end up the table.

.. important::

   Input data points do not need to be uniformly spaced, or even sorted.
   Users will insert rows one by one; ``LookupTable1D`` has functions for sorting and regularizing the table. 

Data modification
_________________

Scaling
^^^^^^^

To scale data in a particular column, use

.. code-block:: c++

   template <size_t K>
   void scale(const T& a_scale) noexcept;

where ``K`` is the column to be scaled.

Column swapping
^^^^^^^^^^^^^^^

To swap columns, use

.. code-block:: c++

   void swap(const size_t a_columnOne, const size_t a_columnTwo) noexcept;

.. code-block:: text

   2.0  2.0  3.0
   1.0  5.0  6.0

and one calls ``swap(1,2)`` the final table becomes

.. code-block:: text

   2.0  3.0  2.0
   1.0  6.0  5.0

.. warning::

   Note that swapping two columns destroys the sorting and one will need to set the independent variable again afterwards.   

Range truncation
^^^^^^^^^^^^^^^^

To restrict the data range, call

.. code-block:: c++

   void truncate(const T a_min, const T a_max, const size_t a_variable);

where ``a_min`` and ``a_max`` are the permissible ranges for data in the input column (``a_variable``).
Data outside these ranges is discarded from the table. 

Regularize table
________________

When regularizing the table, the potentially irregularly spaced raw data is interpolated onto a regular grid.
The user must specify:

#. The independent variable.
#. Number of grid points in regular grid.
#. Grid point spacing.

A ``LookupTable1D`` is regularized through

.. code-block:: c++

   inline void
   prepareTable(const size_t& a_independentVariable, const size_t& a_numPoints, const LookupTable::Spacing& a_spacing);

Here, ``a_independentVariable`` is the independent variable and ``a_numPoints`` is the number of grid points in the regularized table.
Two different spacings are supported: ``LookupTable::Spacing::Uniform`` and ``LookupTable::Spacing::Exponential``.

**Uniform spacing**

With uniform spacing, grid points in the table are spaced as

.. math::

   x_i = x_{\textrm{min}} + \frac{i}{N-1}\left(x_{\textrm{max}} - x_{\textrm{min}}\right),\quad i\in[0,N-1]

where :math:`x_{\textrm{min}}` and :math:`x_{\textrm{max}}` is the minimum and maximum data range for the independent variable (i.e., column).

**Exponential spacing**

If grid points are exponentially spaced then

.. math::

   x_i = x_{\textrm{min}}\left(\frac{x_{\textrm{max}}}{x_{\textrm{min}}}\right)^{\frac{i}{N-1}}, \quad i\in[0,N-1].

.. warning::

   Note that one must have :math:`x_{\textrm{min}} > 0` when using exponentially spaced points. 

Data interpolation
__________________

To retrieve data from one of the columns, one can fetch either a specific value in a row, or the entire row. 

.. code-block:: c++

   // For fetching column K
   template<size_t K>
   T interpolate(const T& a_x) const noexcept;

   // For fetching the entire row
   std::array<T, N> interpolate(const const T& a_x) const noexcept;

In the above, the template parameter ``K`` is the column to retrieve and ``a_x`` is the value of the independent variable.

.. important::

   ``LookupTable1D`` will *always* use piecewise linear interpolation between two grid points.

For example, consider table regularized and sorted along the middle column:

.. code-block:: text

   2.0  1.0  3.0
   1.0  3.0  6.0
   1.0  5.0  4.0

To retrieve an interpolated value for ``x=2.0`` in the third column we call

.. code-block:: c++

   LookupTable1D<Real, 2> myTable,

   const T val = myTable.interpolate<2>(2.0);

which will return a value of 4.5 (linearly interpolated).

Out-of-range strategy
_____________________

Extrapolation outside the valid data range is determined by a user-specified strategy.
When calling the ``interpolate`` function with an argument that exceeds the bounds of the raw or regular data, the range strategy is either:

* Return the value at the endpoint.
* Extrapolate from the endpoint.

To set the range strategy one can use

.. code-block:: c++

   void setRangeStrategyLo(const LookupTable::OutOfRangeStrategy& a_strategy) noexcept;
   void setRangeStrategyHi(const LookupTable::OutOfRangeStrategy& a_strategy) noexcept;   

where ``a_strategy`` must be either of

* ``LookupTable::OutOfRangeStrategy::Constant``.
* ``LookupTable::OutOfRangeStrategy::Interpolate``.

The default behavior is ``LookupTable::OutOfRangeStrategy::Constant``.

Viewing tables
______________

For debugging purposes, ``LookupTable1D`` can write the internal data to an output stream or a file through various member functions:

.. code-block:: c++

   void writeRawData(const std::string& a_file) const noexcept;
   void writeStructuredData(const std::string& a_file) const noexcept;

   void outputRawData(std::ostream& a_ostream = std::cout) const noexcept;
   void outputStructuredData(std::ostream& a_ostream = std::cout) const noexcept;		

These functions will print the table (either raw or regularized) to an output stream or file.


