.. _Chap:Contributions:

Contributions
=============

We welcome feedback, bug reports, and code contributions to ``chombo-discharge``.
For questions or general discussion use the issue tracker or discussion tab at
`github.com/chombo-discharge <https://github.com/chombo-discharge/chombo-discharge>`_.

.. contents:: On this page
   :local:
   :depth: 2

Coding conventions
------------------

``chombo-discharge`` follows a consistent set of coding conventions.
A more detailed overview of the code structure is given in :ref:`Chap:CodeStandard`.

Naming
......

+--------------------+------------------+-----------------------------------+
| Entity             | Convention       | Example                           |
+====================+==================+===================================+
| Classes/structs    | ``PascalCase``   | ``AmrMesh``, ``CdrSolver``        |
+--------------------+------------------+-----------------------------------+
| Member variables   | ``m_camelCase``  | ``m_phi``, ``m_amrMesh``          |
+--------------------+------------------+-----------------------------------+
| Function arguments | ``a_camelCase``  | ``a_phi``, ``a_level``            |
+--------------------+------------------+-----------------------------------+
| Local variables    | ``camelCase``    | ``numCells``, ``loCorner``        |
+--------------------+------------------+-----------------------------------+
| Static variables   | ``s_camelCase``  | ``s_defaultOrder``                |
+--------------------+------------------+-----------------------------------+

File naming
...........

Every file must be prefixed with ``CD_``:

* Header files: ``CD_MyClass.H`` (capital ``.H``, not ``.h`` or ``.hpp``).
* Implementation files: ``CD_MyClass.cpp``.
* Template/inline implementations: ``CD_MyClassImplem.H``.

Header guards
.............

Header guards must be fully uppercased, matching the filename with dots replaced by underscores:

.. code-block:: c++

   #ifndef CD_MYCLASS_H
   #define CD_MYCLASS_H
   // ...
   #endif

Include ordering
................

Group includes into labelled sections in the following order:

.. code-block:: c++

   // Std includes
   #include <memory>
   #include <vector>

   // Chombo includes
   #include <LevelData.H>
   #include <EBISBox.H>

   // Our includes
   #include <CD_AmrMesh.H>
   #include <CD_NamespaceHeader.H>

Omit a group entirely (including its label) when it has no entries.

Namespace
.........

All ``chombo-discharge`` code lives inside the ``ChomboDischarge`` namespace, opened and closed with:

.. code-block:: c++

   #include <CD_NamespaceHeader.H>
   // ... declarations ...
   #include <CD_NamespaceFooter.H>

SPDX/REUSE compliance
---------------------

Every ``.H`` and ``.cpp`` file must carry a REUSE-compliant SPDX block at the very top, before any other content:

.. code-block:: c++

   /*
    * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
    *
    * SPDX-License-Identifier: GPL-3.0-or-later
    */

Use ``reuse annotate`` to apply the header automatically:

.. code-block:: bash

   reuse annotate --copyright "2021-2026 SINTEF Energy Research" \
                  --license GPL-3.0-or-later \
                  --style c Source/MyModule/CD_MyClass.H \
                  Source/MyModule/CD_MyClass.cpp

Do **not** add SPDX headers to ``.options`` or ``.inputs`` files; those are covered by glob
patterns in ``REUSE.toml``.

Documentation standards
-----------------------

``chombo-discharge`` uses Doxygen to generate an API reference.
The CI ``Build-documentation`` job runs ``doxygen`` with ``WARN_AS_ERROR = FAIL_ON_WARNINGS``,
so any documentation warning causes a build failure.

File-level block
................

Every ``.H`` and ``.cpp`` file must have a ``/**`` Doxygen block immediately after the SPDX
block, containing ``@file``, ``@brief``, and ``@author``:

.. code-block:: c++

   /**
     @file   CD_MyClass.H
     @brief  Declaration of MyClass.
     @author Robert Marskar
   */

For ``.cpp`` files the ``@brief`` must read ``"Implementation of CD_MyClass.H"``.

.. important::

   All Doxygen comment blocks **must** use the ``/**`` opening delimiter.

Function documentation
......................

Every function — public, protected, or private — must have at minimum a ``@brief``.
All parameters must be tagged with a direction:

* ``@param[in]``    — input-only parameter.
* ``@param[out]``   — output-only parameter.
* ``@param[in,out]`` — parameter that is read and modified (note: ``@param[inout]`` is invalid).

Non-void functions must have a ``@return`` tag.

.. code-block:: c++

   /**
     @brief Solve the system on a single level.
     @param[in,out] a_phi    Solution data; initial guess on entry.
     @param[in]     a_rhs    Right-hand side data.
     @param[in]     a_level  AMR level index.
     @return Residual norm after the solve.
   */
   Real
   solve(LevelData<EBCellFAB>& a_phi,
         const LevelData<EBCellFAB>& a_rhs,
         const int a_level);

Hiding implementation details
.............................

Template specialisations, explicit instantiations, and other implementation artefacts in
``.cpp`` files that should not appear in the API reference must be wrapped with:

.. code-block:: c++

   /// @cond DOXYGEN_SKIP
   template <>
   void MyClass::specialised<2>() { ... }
   /// @endcond

HDF5-guarded functions
......................

Documentation comments for functions inside ``#ifdef CH_USE_HDF5`` blocks must be placed
**inside** the ``#ifdef`` block.  Comments placed before the ``#ifdef`` are associated with
the next visible symbol when ``CH_USE_HDF5`` is not defined, producing spurious warnings:

.. code-block:: c++

   // Correct — comment is inside the guard
   #ifdef CH_USE_HDF5
   /**
     @brief Write checkpoint data.
     @param[out] a_handle HDF5 file handle.
     @param[in]  a_level  Grid level.
   */
   virtual void
   writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;
   #endif

Code quality tools
------------------

Pre-commit hooks
................

``chombo-discharge`` uses `pre-commit <https://pre-commit.com>`_ to enforce code quality
checks automatically before each commit.

Install and activate the hooks once per clone:

.. code-block:: bash

   pip install pre-commit
   pre-commit install

Run all hooks manually across the entire repository:

.. code-block:: bash

   pre-commit run --all-files

The ``sphinx-build`` hook is on a **manual stage** and must be invoked explicitly:

.. code-block:: bash

   pre-commit run sphinx-build --hook-stage manual

It builds the HTML documentation with ``-W --keep-going``, which treats all Sphinx warnings
as errors and collects every warning before reporting failure.  Run this before opening a
pull request if you have made changes to RST files.

The following hooks are active (see ``.pre-commit-config.yaml``):

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Hook
     - What it checks
     - Notes
   * - ``clang-format``
     - C++ source formatting
     - Enforces ``.clang-format`` style; fails if any diff is produced
   * - ``clang-tidy``
     - Static analysis of staged ``.cpp`` files
     - Requires ``compile_commands.json``.
   * - ``cppcheck``
     - Lightweight static analysis of staged ``.cpp`` files
     - No compilation database required; checks for warnings, performance, and portability issues
   * - ``reuse``
     - REUSE/SPDX licence compliance
     - Runs ``reuse lint`` on the whole repo
   * - ``codespell``
     - Spelling in source, docs, and physics files
     - Word exceptions are listed in ``.codespellignore``
   * - ``format-input-files``
     - Banner comment format in ``.options``/``.inputs`` files
     - Header must be a 100-character ``====`` rule followed by class name
   * - ``doxygen-check``
     - Runs ``doxygen Docs/doxygen.conf``
     - Fails on any warning (``WARN_AS_ERROR = FAIL_ON_WARNINGS``)
   * - ``sphinx-build``
     - Builds the HTML documentation; fails on any warning or error
     - **Manual stage only** — not triggered on commit; invoke explicitly (see below)

clang-format
............

All C++ source files are formatted with ``clang-format`` v18 using the configuration in
``.clang-format`` at the repository root.  Key settings:

* Language standard: C++11.
* Indent width: 2 spaces (no tabs).
* Column limit: 120 characters.
* Pointer alignment: left (``int* a_ptr``).
* Brace wrapping: opening braces on new lines for classes, enums, structs, and functions;
  ``else`` on new line after closing brace.

To format all source files in the repository:

.. code-block:: bash

   find Source Physics Geometries Exec \( -name "*.H" -o -name "*.cpp" \) \
       -exec clang-format -i {} +

Or rely on the pre-commit hook, which runs automatically on staged files.

clang-tidy
..........

Static analysis is performed with ``clang-tidy`` using the configuration in ``.clang-tidy``.
The active check families target C++14: ``bugprone-*``, ``clang-analyzer-*``,
``cppcoreguidelines-*``, ``misc-*``, a curated subset of ``modernize-*``
(``loop-convert``, ``make-unique``, ``use-auto``, ``use-default-member-init``,
``use-nullptr``, ``use-override``), ``performance-*``, and ``readability-*``,
with noisy or inapplicable checks disabled.
Only ``.cpp`` files under ``Source/``, ``Physics/``, ``Exec/``, and ``Geometries/``
are analysed directly; headers are checked transitively.

Generating ``compile_commands.json``
.....................................

``clang-tidy`` requires a compilation database.
``chombo-discharge`` uses a GNU Make build system, so the database must be generated
with `Bear <https://github.com/rizsotto/Bear>`_, which intercepts the compiler invocations
during a normal build:

.. code-block:: bash

   sudo apt install bear        # install Bear
   cd $DISCHARGE_HOME
   bear -- make -j$(nproc) DIM=2 discharge-lib

This produces ``compile_commands.json`` in ``$DISCHARGE_HOME``.
Regenerate it whenever you add or remove source files, or change compiler flags.

.. note::

   The compilation database must cover the files you want to analyse.
   If you are working on a single physics module, build that module explicitly:

   .. code-block:: bash

      bear -- make -j$(nproc) DIM=2 electrostatics

Running clang-tidy
..................

**Via the pre-commit hook** (staged ``.cpp`` files only):

.. code-block:: bash

   pre-commit run clang-tidy

The hook runs automatically on ``git commit``.  It passes only the staged ``.cpp``
files in the four project directories to ``clang-tidy``.
Headers are analysed transitively.
No action is taken if no ``.cpp`` files are staged.

**In the CI pipeline** (changed files in the pull request):

The ``Clang-tidy`` GitHub Actions job builds a compilation database with Bear, then
runs ``clang-tidy`` on only the ``.cpp`` files that differ between the PR branch and
``main``.  This keeps CI fast regardless of repository size.  If a pull request
touches no ``.cpp`` files the job exits immediately with success.

cppcheck
........

``cppcheck`` provides fast, lightweight static analysis that does not require a
compilation database.  It catches memory errors, use-after-free, null pointer
dereferences, performance anti-patterns, and portability issues directly from
source files in seconds.

The pre-commit hook runs automatically on staged ``.cpp`` files:

.. code-block:: bash

   pre-commit run cppcheck

Or invoke ``cppcheck`` manually on specific files or directories:

.. code-block:: bash

   cppcheck --enable=warning,performance,portability \
            --suppress=missingInclude --suppress=missingIncludeSystem \
            --inline-suppr --error-exitcode=1 \
            Source/ Physics/ Geometries/ Exec/

The ``--inline-suppr`` flag enables per-line suppression comments in source:

.. code-block:: cpp

   // cppcheck-suppress someWarningId
   risky_call();

The CI ``Cppcheck`` job runs with a full compilation database (generated by Bear), so
it has complete include paths and can report more precisely than the pre-commit hook.

codespell
.........

``codespell`` checks for common spelling mistakes in source code, documentation, and physics
configuration files.  False positives (technical terms, variable names) are listed one per
line in ``.codespellignore``.  Entries must be lowercase.

If ``codespell`` flags a word that is intentional, add it to ``.codespellignore`` in your
pull request.

CheckDocs.py
............

``CheckDocs.py`` (at the repository root) cross-references your local changes against the
Sphinx documentation.  It compares the current branch against ``main``, collects every
modified ``.H``, ``.cpp``, ``.options``, ``.inputs``, and ``.md`` file, and then scans all
RST files for ``.. literalinclude::`` directives that point to any of those files.

Run it manually at any time:

.. code-block:: bash

   python3 CheckDocs.py

Example output when a changed file is included in the docs:

.. code-block:: text

   Found literalincludes referencing changed files:

   Docs/Sphinx/source/Base/MyPage.rst:42 -> includes '../../Source/MyModule/CD_MyClass.H'

The hook **always passes** — it never blocks a commit.  Its purpose is to remind you to
review the listed documentation pages and check that the included code excerpts still render
correctly after your changes.

Options and inputs file format
..............................

Every ``.options`` and ``.inputs`` file must begin with a three-line banner:

.. code-block:: text

   # ====================================================================================================
   # ClassName class options
   # ====================================================================================================

``ClassName`` must match the C++ class name exactly.  The ``format-input-files`` hook
enforces this automatically.

Submitting a pull request
-------------------------

1. **Create a feature branch** from ``main``:

   .. code-block:: bash

      git checkout main
      git pull
      git checkout -b feature/my-feature

2. **Develop the feature.** Follow the coding conventions above.
   Add or update Sphinx documentation and Doxygen comments as appropriate.

3. **Run pre-commit** to catch formatting and documentation issues before committing:

   .. code-block:: bash

      pre-commit run --all-files

   Fix any reported issues, then stage and commit:

   .. code-block:: bash

      git add Source/MyModule/CD_MyClass.H Source/MyModule/CD_MyClass.cpp
      git commit -m "Add MyClass for ..."

4. **Run the test suite** locally in debug mode to check for assertion failures or memory
   errors (see :ref:`Chap:Testing`):

   .. code-block:: bash

      cd Exec/Tests/MyTest
      make -j$(nproc) DIM=2 DEBUG=TRUE
      ./main.*.ex regression2d.inputs

   For memory leak checking with ``valgrind``:

   .. code-block:: bash

      valgrind --leak-check=full --track-origins=yes ./main.*.ex regression2d.inputs

5. **Push the branch** to GitHub:

   .. code-block:: bash

      git push --set-upstream origin feature/my-feature

6. **Open a pull request** against ``main`` via the GitHub web interface.
   Mark it as *Draft* if the work is still in progress.
   Provide a clear description of what was changed and why.

7. **Ensure all CI checks pass** (see :ref:`Chap:CI` below).  Address any failures before
   requesting a review.

.. important::

   Squash your commits into a single commit before the pull request is merged.

.. _Chap:CI:

Continuous integration
----------------------

``chombo-discharge`` uses GitHub Actions (``CI.yml``) to validate every pull request.
All jobs must pass before a pull request can be merged into ``main``.

.. list-table::
   :header-rows: 1
   :widths: 20 15 50

   * - Job
     - Runs on
     - Purpose
   * - ``Formatting``
     - Ubuntu
     - ``clang-format`` diff check; fails if any file is not correctly formatted
   * - ``REUSE``
     - Ubuntu
     - ``reuse lint``; fails if any file is missing an SPDX header or REUSE.toml entry
   * - ``Codespell``
     - Ubuntu
     - ``codespell`` spelling check across source, docs, and physics files
   * - ``Clang-tidy``
     - Ubuntu
     - Builds a compilation database with Bear, then runs ``clang-tidy`` (via ``clang-tidy-cache``) in parallel on only the ``.cpp`` files changed in the pull request; skipped if no ``.cpp`` files changed
   * - ``Cppcheck``
     - Ubuntu
     - Builds a compilation database with Bear, then runs ``cppcheck`` on only the ``.cpp`` files changed in the pull request; skipped if no ``.cpp`` files changed
   * - ``Linux-GNU``
     - Ubuntu
     - Full 2D and 3D debug builds with GCC; depends on Formatting, REUSE, and Codespell
   * - ``Linux-oneAPI``
     - Ubuntu
     - Full 2D and 3D debug builds with Intel oneAPI; same dependencies as Linux-GNU
   * - ``Build-documentation``
     - Ubuntu
     - Doxygen HTML + Sphinx HTML; validates all ``.. literalinclude::`` paths
   * - ``CI-passed``
     - —
     - Final gate job; all other jobs must succeed before this passes

CI runs in 2-D without MPI and with ``DEBUG=TRUE`` so that all assertions are active.
Executables are not run in CI; only compilation is checked.

Upon merging to ``main``, the documentation is rebuilt and deployed to
`GitHub Pages <https://pages.github.com/>`_, keeping the online HTML, PDF, and doxygen
API always in sync with the latest release.

Bug reports
-----------

``chombo-discharge`` is probably not bug-free.
If you encounter unexpected behaviour, please open an issue at
`github.com/chombo-discharge/chombo-discharge/issues <https://github.com/chombo-discharge/chombo-discharge/issues>`_.
Include a minimal reproducible example (inputs file and build flags) where possible.
