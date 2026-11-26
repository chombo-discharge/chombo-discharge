/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ScratchPool.H
 * @brief  Declaration of the ScratchPool class and the Scratch handle it hands out.
 * @author Robert Marskar
 */

#ifndef CD_SCRATCHPOOL_H
#define CD_SCRATCHPOOL_H

// Std includes
#include <memory>
#include <ostream>
#include <tuple>
#include <utility>
#include <vector>

// Chombo includes
#include <LevelData.H>
#include <EBCellFAB.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_NamespaceHeader.H>

class ScratchPool;

/**
 * @brief Single owner of a temporary that belongs to a ScratchPool.
 * @details A Scratch is obtained from AmrMesh::fetchScratch (or ScratchPool::acquire) and returns the temporary to
 * its pool when it is destroyed. It can be moved but not copied, so a temporary has exactly one owner at a time. It
 * has no default constructor and cannot be assigned. It must never be stored in a class member (the compiler does not
 * prevent it, e.g. through a constructor's initializer list or std::optional): the entry would stay in use for the
 * object's lifetime, and the pool aborts at the next regrid. Data that must live longer than the function that
 * fetches it is fetched by a caller that spans all of its uses, and passed down by reference.
 *
 * The handle converts implicitly to T&, so it can be passed directly to non-template functions that take the data by
 * reference. A function template that deduces the data type from its argument (e.g. AmrMesh::copyData or the
 * templated DataOps functions) does not see the conversion, and takes *handle.
 * References and aliases (e.g. Vector<LevelData<T>*>) taken from it are views: they must not outlive the handle.
 * Copying the underlying EBAMRData into an object that outlives the handle is an error. The pool detects shallow
 * copies of the EBAMRData levels when the handle is destroyed. It cannot detect aliases built with AmrMesh::alias,
 * EBAMRData::slice or aliasLevelData, so those must never outlive the handle either.
 *
 * The data is uninitialized when fetched. It has exactly the number of components and ghost cells that were asked
 * for.
 * @tparam T Data type. One of the types ScratchPool stores.
 */
template <typename T>
class Scratch
{
public:
  /**
   * @brief Disallowed. A Scratch always refers to a temporary.
   */
  Scratch() = delete;

  /**
   * @brief Disallowed. A temporary has one owner.
   * @param[in] a_other Other handle
   */
  Scratch(const Scratch& a_other) = delete;

  /**
   * @brief Move constructor. Takes over the temporary, and leaves a_other empty.
   * @param[in,out] a_other Other handle. Its destructor does nothing afterwards.
   */
  Scratch(Scratch&& a_other) noexcept;

  /**
   * @brief Disallowed. A temporary has one owner.
   * @param[in] a_other Other handle
   * @return Not applicable
   */
  Scratch&
  operator=(const Scratch& a_other) = delete;

  /**
   * @brief Disallowed, so that a Scratch cannot be stored in a member that outlives the function that fetched it.
   * @param[in] a_other Other handle
   * @return Not applicable
   */
  Scratch&
  operator=(Scratch&& a_other) = delete;

  /**
   * @brief Destructor. Returns the temporary to its pool.
   */
  ~Scratch() noexcept;

  /**
   * @brief Access the data.
   * @return The temporary
   */
  operator T&() noexcept;

  /**
   * @brief Access the data.
   * @return The temporary
   */
  operator const T&() const noexcept;

  /**
   * @brief Access the data.
   * @return The temporary
   */
  T&
  operator*() noexcept;

  /**
   * @brief Access the data.
   * @return The temporary
   */
  const T&
  operator*() const noexcept;

  /**
   * @brief Access the data.
   * @return Pointer to the temporary
   */
  T*
  operator->() noexcept;

  /**
   * @brief Access the data.
   * @return Pointer to the temporary
   */
  const T*
  operator->() const noexcept;

private:
  friend class ScratchPool;

  /**
   * @brief Construct a handle to an entry of a pool. The entry is already marked in use.
   * @param[in] a_pool  Pool that owns the temporary
   * @param[in] a_entry Index of the entry in the pool's list for type T
   */
  Scratch(ScratchPool& a_pool, const int a_entry) noexcept;

  /**
   * @brief Abort if the handle was moved from. Checked in every build, not only with assertions on.
   * @param[in] a_pool Pool of the handle
   */
  static void
  assertOwner(const ScratchPool* const a_pool) noexcept;

  /**
   * @brief Pool that owns the temporary. Null after the handle was moved from.
   */
  ScratchPool* m_pool;

  /**
   * @brief Index of the entry in the pool's list for type T.
   */
  int m_entry;
};

/**
 * @brief Pool of temporaries on one Realm, reused between the functions that need them.
 * @details The pool keeps every temporary it has built, together with the key it was built for (phase, level, number
 * of components, ghost cells) and whether it is in use. acquire() hands out a free entry with a matching key, or
 * stores a new one, and returns a Scratch that gives it back when destroyed. The pool therefore grows to the largest
 * number of temporaries in use at the same time, and keeps them, Copiers and communication buffers included, until
 * the grids change.
 *
 * The pool only does the bookkeeping; the caller (AmrMesh) builds the data. Entries are built on the Realm's current
 * grids, so the owner must call clear() whenever the grids change. clear() aborts if a temporary is still in use, so
 * a Scratch can never refer to data on old grids.
 *
 * Every call aborts when made inside an OpenMP parallel region. Temporaries are fetched and returned at the level of
 * the hierarchy, never per box, so the pool needs no lock.
 */
class ScratchPool
{
public:
  /**
   * @brief Constructor. The pool is empty.
   */
  ScratchPool() noexcept;

  /**
   * @brief Disallowed. Handles point to their pool.
   * @param[in] a_other Other pool
   */
  ScratchPool(const ScratchPool& a_other) = delete;

  /**
   * @brief Disallowed. Handles point to their pool.
   * @param[in] a_other Other pool
   */
  ScratchPool(ScratchPool&& a_other) = delete;

  /**
   * @brief Disallowed. Handles point to their pool.
   * @param[in] a_other Other pool
   * @return Not applicable
   */
  ScratchPool&
  operator=(const ScratchPool& a_other) = delete;

  /**
   * @brief Disallowed. Handles point to their pool.
   * @param[in] a_other Other pool
   * @return Not applicable
   */
  ScratchPool&
  operator=(ScratchPool&& a_other) = delete;

  /**
   * @brief Destructor. Aborts if a temporary is still in use.
   */
  ~ScratchPool() noexcept;

  /**
   * @brief Hand out a temporary with the given key.
   * @details Returns a free entry of type T with the same phase, level, number of components and number of ghost cells.
   * If there is none, calls a_build and stores its result as a new entry.
   * @param[in] a_phase Phase the data lives on. Types without a phase pass phase::gas.
   * @param[in] a_level Grid level for single-level data, -1 for data over the whole hierarchy.
   * @param[in] a_nComp Number of components
   * @param[in] a_ghost Number of ghost cells
   * @param[in] a_build Builds a new temporary with this key, and exactly a_ghost ghost cells, on the current grids.
   * @return Handle that owns the temporary until it is destroyed.
   * @tparam T     Data type
   * @tparam Build Callable () -> std::unique_ptr<T>. Only called when no free entry matches.
   */
  template <typename T, typename Build>
  Scratch<T>
  acquire(const phase::which_phase a_phase,
          const int                a_level,
          const int                a_nComp,
          const int                a_ghost,
          const Build&             a_build);

  /**
   * @brief Delete all temporaries. Called whenever the grids change.
   * @details Aborts if a temporary is still in use.
   */
  void
  clear() noexcept;

  /**
   * @brief Check that no temporary is in use.
   * @return True if every temporary has been returned.
   */
  bool
  allReturned() const noexcept;

  /**
   * @brief Write the entries of the pool, one line per entry, with their keys and whether they are in use.
   * @param[in,out] a_os     Output stream
   * @param[in]     a_prefix Text written at the start of every line
   */
  void
  report(std::ostream& a_os, const std::string& a_prefix) const noexcept;

private:
  template <typename T>
  friend class Scratch;

  /**
   * @brief A temporary, the key it was built for, and whether it is in use.
   */
  template <typename T>
  struct Entry
  {
    /**
     * @brief The temporary.
     */
    std::unique_ptr<T> data;

    /**
     * @brief Phase. phase::gas for types without a phase.
     */
    phase::which_phase phase;

    /**
     * @brief Grid level for single-level data, -1 for data over the whole hierarchy.
     */
    int level;

    /**
     * @brief Number of components.
     */
    int nComp;

    /**
     * @brief Number of ghost cells.
     */
    int ghost;

    /**
     * @brief True while a Scratch owns the temporary.
     */
    bool inUse;

    /**
     * @brief Each level of the temporary as built: the LevelData it points to, its number of components and its ghost
     * cells. giveBack() checks that the data still has these levels, so that reallocating or replacing them through the
     * handle is caught.
     */
    std::vector<std::tuple<const void*, int, IntVect>> layout;
  };

  /**
   * @brief The entries, one list per type. Searched front to back.
   */
  std::tuple<std::vector<Entry<EBAMRCellData>>,
             std::vector<Entry<EBAMRFluxData>>,
             std::vector<Entry<EBAMRIVData>>,
             std::vector<Entry<EBAMRFAB>>,
             std::vector<Entry<MFAMRCellData>>,
             std::vector<Entry<MFAMRFluxData>>,
             std::vector<Entry<MFAMRIVData>>,
             std::vector<Entry<LevelData<EBCellFAB>>>>
    m_entries;

  /**
   * @brief The list of entries of type T.
   * @return m_entries' list for type T
   */
  template <typename T>
  std::vector<Entry<T>>&
  entries() noexcept;

  /**
   * @brief The list of entries of type T.
   * @return m_entries' list for type T
   */
  template <typename T>
  const std::vector<Entry<T>>&
  entries() const noexcept;

  /**
   * @brief Access the temporary of an entry.
   * @param[in] a_entry Index of the entry in the list for type T
   * @return The temporary
   */
  template <typename T>
  T&
  get(const int a_entry) noexcept;

  /**
   * @brief Mark an entry free again. Called by the Scratch that owns it.
   * @details For data over the whole hierarchy, aborts if a level of the temporary is still referenced from outside
   * the pool, i.e. if a copy of the EBAMRData outlived the handle.
   * @param[in] a_entry Index of the entry in the list for type T
   */
  template <typename T>
  void
  giveBack(const int a_entry) noexcept;

  /**
   * @brief Abort if called inside an OpenMP parallel region.
   * @param[in] a_caller Name of the calling function, for the error message
   */
  static void
  assertSerial(const char* a_caller) noexcept;

  /**
   * @brief Name of a stored type, for messages.
   * @return Type name
   */
  template <typename T>
  static const char*
  typeName() noexcept;

  /**
   * @brief Write the entries of one type, one line per entry.
   * @param[in]     a_list   Entries of type T
   * @param[in,out] a_os     Output stream
   * @param[in]     a_prefix Text written at the start of every line
   */
  template <typename T>
  void
  reportType(const std::vector<Entry<T>>& a_list, std::ostream& a_os, const std::string& a_prefix) const noexcept;

  /**
   * @brief Call a function on each level of a temporary, with the level's address, components and ghost cells.
   * @details A null level is passed as a null address, -1 components and zero ghost cells.
   * @param[in] a_data  Temporary
   * @param[in] a_visit Callable (int level, const void* address, int components, const IntVect& ghost) -> void
   */
  template <typename T, typename Visit>
  static void
  forEachLevel(const T& a_data, const Visit& a_visit) noexcept;

  /**
   * @brief Number of levels of a temporary. 1 for single-level data.
   * @param[in] a_data Temporary
   * @return Number of levels
   */
  template <typename T>
  static int
  numLevels(const T& a_data) noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_ScratchPoolImplem.H>

#endif
