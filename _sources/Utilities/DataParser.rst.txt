.. _Chap:DataParser:

Data parsing
************

Routines for reading column data into ``LookupTable`` are defined under the ``DataParser`` namespace.
For example, for reading two columns into a lookup table (see :ref:`Chap:LookupTable`):

.. code-block:: c++

  // Read all rows
  LookupTable1D<Real, 1>
  simpleFileReadASCII(const std::string       a_fileName,
                      const int               a_xColumn     = 0,
                      const int               a_yColumn     = 1,
                      const std::vector<char> a_ignoreChars = {'#', '/'});

  // Specify rows where to start and stop reading data
  LookupTable1D<Real, 1>
  fractionalFileReadASCII(const std::string       a_fileName,
                          const std::string       a_startRead,
                          const std::string       a_stopRead,
                          const int               a_xColumn     = 0,
                          const int               a_yColumn     = 1,
                          const std::vector<char> a_ignoreChars = {'#', '/'});


.. tip::

   The ``DataParser`` C++ API is found at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/namespaceDataParser.html>`_.
