.. _Chap:LookupTable1D:

Lookup tables
=============

``LookupTable1D`` is a class for looking up and interpolation data stored in a row-column format.
It is used in order to easily retrieve input data that can be stored in table formats.

.. important::

   LookupTable1D is used for data lookup *in one independent variable*.
   It does not support higher-dimensional data interpolation. 

The class is templated as

.. code-block::

   template <int N>
   class LookupTable1D

where the template parameter ``N`` indicates the number of columns in the data holders.
Internally, the data is stored as an ``std::vector<std::array<Real, N> >`` where the vector entries are rows and the ``std::array<Real, N>`` are data in each row.
Thus, a table ``LookupTable1D<2>`` always has two columns.

The ``LookupTable1D`` is used on regularly spaced data (for performance reasons).
Although the user can fill irregularly spaced data into ``LookupTable1D``, the class has routines for making that data regularly spaced and sorted along of its columns.
Usage of ``LookupTable1D`` will therefore consist of the following:

#. Add data rows into the table.
#. Swap columns if necessary.
#. Restrict data ranges if necessary. 
#. Sort the table along of it's columns, smallest to largest.
#. Regularize the table with a specified number of grid points.
#. *Retrieve data*.

The steps for these processes are explained in detail below.

Inserting data
--------------

To add data to the table, one will use the member function

.. code-block:: c++

   template <typename... Ts>
   inline
   void addEntry(const Ts&... x);

where the parameter pack must have ``N`` entries.
For example, to add two rows of data to a table:

.. code-block:: c++

   LookupTable1D<3> myTable;

   myTable.add(4.0, 5.0, 6.0);
   myTable.add(1.0, 2.0, 3.0);   

This will insert two new rows at the end up the table.

.. important::

   Input data points do not need to be uniformly spaced, or even sorted.
   Users will insert rows one by one; ``LookupTable1D`` has functions for sorting and regularizing the table. 

Restricting ranges
------------------

To restrict the data range, call

.. code-block:: c++

   void setRange(const Real a_min, const Real a_max, const int a_independentVariable)

where ``a_min`` and ``a_max`` are the permissible ranges for data in the input column (``a_independentVariable``).
Data outside these ranges are discarded from the table. 



Independent variable
--------------------

To select an independent variable, the user will select one of the columns and sort the data along that column.
The C++ code for this is

.. code-block:: c++

   template <int N>
   void LookupTable1D<N>::sort(const int a_independentVariable);

where the input integer indicates the column (i.e., independent variable) used for sorting.
Note that the sorting *is always from smallest to largest value*.
Thus, if one has two rows

.. code-block:: text

   1.0  5.0  6.0
   2.0  2.0  3.0

and one calls ``LookupTable1D<N>::sort(1)`` the final table becomes

.. code-block:: text

   2.0  2.0  3.0
   1.0  5.0  6.0

Note that the second column now becomes the independent variable.

Swapping columns
----------------

Columns can be swapped by calling ``LookupTable1D<N>::swap(int, int)``, which will swap two of the columns.
For example if the original data is

.. code-block:: text

   2.0  2.0  3.0
   1.0  5.0  6.0

and one calls ``swap(1,2)`` the final table becomes

.. code-block:: text

   2.0  3.0  2.0
   1.0  6.0  5.0

Note that swapping two columns destroys the sorting and one will need to set the independent variable again afterwards.

Regularize table
----------------

To regularize the table the user must first determine if the grid points should be uniformly spaced or exponentially spaced in the independent variable.

Setting grid point spacing
__________________________

The user can set the spacing by calling

.. code-block:: c++

   void setTableSpacing(const TableSpacing a_spacing);

where ``TableSpacing::Uniform`` and ``TableSpacing::Exponential`` are supported.
For example, to use uniformly or exponentially spaced grid points:

.. code-block:: c++

   LookupTable1D<2> myTable;

   myTable.setTableSpacing(TableSpacing::Uniform);     // For uniformly spaced points
   myTable.setTableSpacing(TableSpacing::Exponential); // For exponentially spaced points   

**Uniform spacing**

With uniform spacing, grid points in the table are spaced as

.. math::

   x_i = x_{\textrm{min}} + \frac{i}{N-1}\left(x_{\textrm{max}} - x_{\textrm{min}}\right),\quad i\in[0,N-1]

where :math:`x_{\textrm{min}}` and :math:`x_{\textrm{max}}` is the minimum and maximum data range for the independent variable (i.e., column).

**Exponential spacing**

If grid points are exponentially spaced then

.. math::

   x_i = x_{\textrm{min}}\left(\frac{x_{\textrm{max}}}{x_{\textrm{min}}}\right)^{\frac{i-1}{N}}, \quad i\in[0,N-1].

Regularizing table
__________________

.. code-block:: c++

   void regularize(const int a_numRows)

which will make the table into a regularly spaced table with ``a_numRows`` rows.
``LookupTable1D`` will always use piecewise linear interpolation when regularizing the table.
Specifying a number of rows that is smaller/larger than the original number of rows will downsample/upsample the table.

.. important::

   When regularizing a table through ``regularize``, the original data table is destroyed. 

Retrieving data
---------------

To retrieve data from one of the columns, one can fetch either a specific value in a row, or the entire row. 

.. code-block:: c++

   // For fetching column K
   template<int K>
   Real getEntry(const Real a_x);

   // For fetching the entire row
   std::array<Real, N> getData(const Real a_x);

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

   LookupTable1D<3> myTable,

   const Real val = myTable.getEntry<2>(2.0);

which will return a value of 4.5 (linearly interpolated). 

Viewing tables
--------------

For debugging purposes, ``LookupTable1D`` can write the internal data to an output stream or a file through the two member functions:

.. code-block:: c++
		
  void dumpTable(std::ostream& a_outputStream = std::cout) const;

  void dumpTable(const std::string a_fileName) const;

For example:

.. code-block:: c++

   LookupTable1D<10> myTable;

   // Dump table to terminal window
   myTable.dumpTable();

   // Dump table to file.
   myTable.dumpTable("myTable.dat");
