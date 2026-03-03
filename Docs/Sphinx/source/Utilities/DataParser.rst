.. _Chap:DataParser:

Data parsing
************

Routines for reading simple column data into :ref:`Chap:LookupTable` are available.
This is typically used, e.g., when parsing transport coefficients or other types of data requires by a computer simulation.

.. tip::

   The ``DataParser`` C++ API is found at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/namespaceDataParser.html>`_.

Currently, three types of file reads are supported:

#. Read column data, where the user can specify the columns that are parsed into a lookup table.

   These files can be arranged in the form

   .. code-block:: text

      0.0 1.0
      1.0 2.0
      2.0 3.0

   It is also possible to restrict which rows that are read, by specifying string identifiers on the line preceding the data and on the line immediately after the data.

   .. important::

      After parsing into a :ref:`Chap:LookupTable`, the user must regularize the table in order to use the interpolation functions.
      See :ref:`Chap:LookupTableRegularize`.
		   
#. Read particle data, where the particle data is in a :math:`(x,y,z,w)` file format, where :math:`w` is the particle weight.
   
   The function signatures for reading this type of data are:

   .. literalinclude:: ../../../../Source/Utilities/CD_DataParser.H
      :language: c++
      :lines: 29-85
      :dedent: 2

#. Read triangle mesh data from PLY or VTK files, returning a list of ``Triangle`` objects
   with per-vertex scalar metadata.
   The file type is inferred from the extension (``.ply``/``.PLY`` or ``.vtk``/``.VTK``).
   A named per-vertex scalar property (e.g., a PLY vertex property or VTK point-data scalar)
   is read as the vertex metadata; if the identifier is not found, a warning is issued and
   metadata defaults to zero.
   Non-triangular facets are skipped with a warning.

   The function signature is:

   .. literalinclude:: ../../../../Source/Utilities/CD_DataParser.H
      :language: c++
      :lines: 87-105
      :dedent: 2

   .. note::

      ``Triangle`` objects are returned in file order.
      To perform efficient spatial queries (e.g., finding the closest triangle to a point),
      pass the result to a ``TriangleCollection``.

