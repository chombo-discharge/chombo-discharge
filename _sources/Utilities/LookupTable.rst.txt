.. _Chap:LookupTable:

LookupTable1D
*************

``LookupTable1D`` is a class for storing and interpolating structured data for lookup in one independent variable.
It is used in order to easily retrieve input data that can be stored in table formats.

.. important::

   LookupTable1D is used for data lookup *in one independent variable*.
   It does not support higher-dimensional data interpolation. 

The class is templated as

.. literalinclude:: ../../../../Source/Utilities/CD_LookupTable1D.H
   :language: c++
   :lines: 28-29

where the template parameter ``N`` indicates the number of dependent variables (``N=1`` yields a compile-time error).
Internally, the data is stored as an ``std::vector<std::array<T, N + 1> >`` where the vector entries are rows and the ``std::array<T, N + 1>`` are the column entries in that row.

.. tip::

   The internal floating point representation defaults to ``Real`` and the number of dependent variables to 1. 

For performance reasons, the ``LookupTable1D`` class is designed to only be used on regularly spaced data (either uniformly or logarithmically spaced).
The user can still pass irregularly spaced data into ``LookupTable1D``, and then regularize the data later (i.e., interpolate it onto a regularly spaced 1D grid).
Usage of ``LookupTable1D`` will therefore consist of the following steps:

#. Add data rows into the table.
#. Swap columns or scale data if necessary.
#. Truncate data ranges if necessary, and specify what happens if the user tries to fetch data outside the valid range.
#. Regularize the table with a specified number of grid points and specified grid spacing.
#. *Retrieve data*, i.e., interpolate data from the table.

All of these steps are discussed below.

Inserting data
==============

To add data to the table, use the member function

.. literalinclude:: ../../../../Source/Utilities/CD_LookupTable1D.H
   :language: c++
   :lines: 59-65
   :dedent: 2

where the parameter pack must have ``N+1`` entries.
For example, to add two rows of data to a table with two dependent variables:

.. code-block:: c++

   LookupTable1D<Real, 2> myTable;

   myTable.addData(4.0, 5.0, 6.0);
   myTable.addData(1.0, 2.0, 3.0);   

This will insert two new rows at the end up the table.

.. important::

   Input data points do not need to be uniformly spaced, or even sorted.
   While users will insert rows one by one; ``LookupTable1D`` has functions for sorting and regularizing the table later.

Data modification
=================

Scaling
-------

To scale data in a particular column, use

.. literalinclude:: ../../../../Source/Utilities/CD_LookupTable1D.H
   :language: c++
   :lines: 84-90
   :dedent: 2

where ``K`` is the column to be scaled.

Column swapping
---------------

To swap columns in the table, use

.. literalinclude:: ../../../../Source/Utilities/CD_LookupTable1D.H
   :language: c++
   :lines: 74-82
   :dedent: 2

For example, if the input data looked like

.. code-block:: text

   2.0  2.0  3.0
   1.0  5.0  6.0

and one calls ``swap(1,2)`` the final table becomes

.. code-block:: text

   2.0  3.0  2.0
   1.0  6.0  5.0

.. warning::

   The swapping function only swaps *raw* data.
   Usually, one must later re-regularize the table, especially if the independent variable was swapped.   

Range truncation
----------------

One can restrict the data range of the table by calling

.. literalinclude:: ../../../../Source/Utilities/CD_LookupTable1D.H
   :language: c++
   :lines: 92-100
   :dedent: 2

where ``a_min`` and ``a_max`` are the permissible ranges for data in the input column (``a_column``).
Data outside these ranges is discarded from the table.
This applies regardless of whether or not ``a_column`` indicates the independent or dependent variables.

.. _Chap:LookupTableRegularize:

Regularize table
================

When regularizing the table, the raw-data (which can be irregularly spaced) is interpolated onto a regular grid.
The user must specify:

#. The independent variable in which one will later interpolate. 
#. Number of grid points in the regular grid.
#. How to space the grid points.

A ``LookupTable1D`` is regularized through

.. literalinclude:: ../../../../Source/Utilities/CD_LookupTable1D.H
   :language: c++
   :lines: 116-123
   :dedent: 2

Here, ``a_independentVariable`` is the independent variable and ``a_numPoints`` is the number of grid points in the regularized table.
Two different spacings are supported: ``LookupTable::Spacing::Uniform`` and ``LookupTable::Spacing::Exponential``.

**Uniform spacing**

With uniform spacing, grid points in the table are spaced as

.. math::

   x_i = x_{\textrm{min}} + \frac{i}{N-1}\left(x_{\textrm{max}} - x_{\textrm{min}}\right),\quad i\in[0,N-1]

where :math:`x_{\textrm{min}}` and :math:`x_{\textrm{max}}` are the minimum and maximum values in the independent variable, and :math:`N` is the number of grid points.

**Exponential spacing**

If grid points are exponentially spaced then the spacing follows a power law:

.. math::

   x_i = x_{\textrm{min}}\left(\frac{x_{\textrm{max}}}{x_{\textrm{min}}}\right)^{\frac{i}{N-1}}, \quad i\in[0,N-1].

.. warning::

   One must have :math:`x_{\textrm{min}} > 0` when using exponentially spaced points.

An example code for regularizing a table is given below:

.. code-block:: c++

   LookupTable1D<Real, 2> myTable;

   myTable.prepareTable(0, 100, LookupTable::Spacing::Exponential);

Data interpolation
==================

To interpolate data from the table, one can fetch either a specific value in a row, or the entire row.
In any case, the values that are returned are linearly interpolated between grid points (in the independent variable).
The function signatures are

.. literalinclude:: ../../../../Source/Utilities/CD_LookupTable1D.H
   :language: c++
   :lines: 125-138
   :dedent: 2

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
=====================

Extrapolation outside the valid data range is determined by a user-specified strategy.
I.e., when calling the ``interpolate`` function with an argument that exceeds the bounds of the raw or regular data, the range strategy is either:

* Return the value at the endpoint of the table.
* Extrapolate from the endpoint.

To set the range strategy one can use

.. literalinclude:: ../../../../Source/Utilities/CD_LookupTable1D.H
   :language: c++
   :lines: 102-114
   :dedent: 2

where ``a_strategy`` must be either of

* ``LookupTable::OutOfRangeStrategy::Constant``.
* ``LookupTable::OutOfRangeStrategy::Interpolate``.

The default behavior is ``LookupTable::OutOfRangeStrategy::Constant``.

Viewing tables
==============

For debugging purposes, ``LookupTable1D`` can write the internal data to an output stream or a file through various member functions:

.. literalinclude:: ../../../../Source/Utilities/CD_LookupTable1D.H
   :language: c++
   :lines: 168-194
   :dedent: 2

These functions will print the table (either raw or regularized) to an output stream or file, and the user can later plot the data in an external plotting tool.


