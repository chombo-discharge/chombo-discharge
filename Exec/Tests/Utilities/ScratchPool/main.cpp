/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Unit test for ScratchPool: temporaries are reused only when free and when their key matches.
 * @author Robert Marskar
 */

// Std includes
#include <memory>
#include <string>
#include <utility>

// Chombo includes
#include <MayDay.H>

// Our includes
#include <CD_Driver.H>
#include <CD_ScratchPool.H>

using namespace ChomboDischarge;

/**
 * @brief Abort with a message if a condition does not hold.
 * @param[in] a_condition Condition that must hold
 * @param[in] a_message   Message describing the check
 */
static void
check(const bool a_condition, const std::string& a_message)
{
  if (!a_condition) {
    MayDay::Abort(("ScratchPool test -- " + a_message).c_str());
  }
}

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  // The pool only does the bookkeeping, so the test builds empty data and counts how often the pool asks for it.
  int numBuilt = 0;

  auto buildCell = [&numBuilt]() -> std::unique_ptr<EBAMRCellData> {
    numBuilt++;

    return std::unique_ptr<EBAMRCellData>(new EBAMRCellData(2));
  };

  auto buildLevel = [&numBuilt]() -> std::unique_ptr<LevelData<EBCellFAB>> {
    numBuilt++;

    return std::unique_ptr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());
  };

  // The pool is destroyed at the end of this block, before finalize() tears down MPI and the timers.
  {
    ScratchPool pool;

    check(pool.allReturned(), "a new pool holds no temporaries in use");

    {
      Scratch<EBAMRCellData> a = pool.acquire<EBAMRCellData>(phase::gas, -1, 1, 2, buildCell);

      check(numBuilt == 1, "the first temporary is built");
      check(!pool.allReturned(), "a held temporary is in use");

      // Same key while the first is held: a second temporary.
      Scratch<EBAMRCellData> b = pool.acquire<EBAMRCellData>(phase::gas, -1, 1, 2, buildCell);

      check(numBuilt == 2, "a held temporary is never handed out twice");
      check(&(*a) != &(*b), "two held temporaries are distinct");
    }

    check(pool.allReturned(), "temporaries are returned at scope exit");

    {
      // Both entries are free again, so these reuse them.
      Scratch<EBAMRCellData> a = pool.acquire<EBAMRCellData>(phase::gas, -1, 1, 2, buildCell);
      Scratch<EBAMRCellData> b = pool.acquire<EBAMRCellData>(phase::gas, -1, 1, 2, buildCell);

      check(numBuilt == 2, "free temporaries with the same key are reused");

      // Any difference in the key gives a new temporary.
      Scratch<EBAMRCellData> c = pool.acquire<EBAMRCellData>(phase::solid, -1, 1, 2, buildCell);
      Scratch<EBAMRCellData> d = pool.acquire<EBAMRCellData>(phase::gas, -1, 3, 2, buildCell);
      Scratch<EBAMRCellData> e = pool.acquire<EBAMRCellData>(phase::gas, -1, 1, 1, buildCell);
      Scratch<EBAMRCellData> f = pool.acquire<EBAMRCellData>(phase::gas, 0, 1, 2, buildCell);

      check(numBuilt == 6,
            "a different phase, level, number of components or number of ghost cells gives a new temporary");

      // Types are kept apart.
      Scratch<LevelData<EBCellFAB>> g = pool.acquire<LevelData<EBCellFAB>>(phase::gas, 0, 1, 2, buildLevel);

      check(numBuilt == 7, "types are kept apart");
    }

    check(pool.allReturned(), "all temporaries are returned at scope exit");

    {
      // The pool now holds two free entries with 2 ghost cells and one with 1 ghost cell (same phase, level and
      // components). Ghost cells are matched exactly: the first request for 1 ghost cell takes the 1-ghost entry, and
      // the second builds a new one rather than taking a 2-ghost entry.
      const int numBefore = numBuilt;

      Scratch<EBAMRCellData> a = pool.acquire<EBAMRCellData>(phase::gas, -1, 1, 1, buildCell);

      check(numBuilt == numBefore, "an entry with the same number of ghost cells is reused");

      Scratch<EBAMRCellData> b = pool.acquire<EBAMRCellData>(phase::gas, -1, 1, 1, buildCell);

      check(numBuilt == numBefore + 1, "an entry with more ghost cells than asked for is not used");

      Scratch<EBAMRCellData> c = pool.acquire<EBAMRCellData>(phase::gas, -1, 1, 3, buildCell);

      check(numBuilt == numBefore + 2, "an entry with fewer ghost cells than asked for is not used");
    }

    check(pool.allReturned(), "all temporaries are returned at scope exit");

    {
      // Moving a handle transfers ownership; the temporary is returned once, by the new owner.
      Scratch<EBAMRCellData> a = pool.acquire<EBAMRCellData>(phase::gas, -1, 1, 2, buildCell);

      EBAMRCellData* const data = &(*a);

      {
        Scratch<EBAMRCellData> b(std::move(a));

        check(&(*b) == data, "a moved handle refers to the same temporary");
        check(!pool.allReturned(), "the temporary is in use while the new owner lives");
      }

      check(pool.allReturned(), "the new owner returns the temporary");
    }

    check(pool.allReturned(), "the moved-from handle returns nothing");

    {
      // The implicit conversion lets a handle be passed where a reference is expected.
      Scratch<EBAMRCellData> a = pool.acquire<EBAMRCellData>(phase::gas, -1, 1, 2, buildCell);

      auto numLevels = [](const EBAMRCellData& a_data) -> int {
        return a_data.size();
      };

      check(numLevels(a) == 2, "a handle converts to a reference to its temporary");
    }

    {
      // A temporary whose levels are only referenced by the pool is returned without complaint. A copy that outlived
      // the handle would raise the level's reference count, and the return would abort.
      auto buildWithLevel = [&numBuilt]() -> std::unique_ptr<EBAMRCellData> {
        numBuilt++;

        std::unique_ptr<EBAMRCellData> data(new EBAMRCellData(1));

        (*data)[0] = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());

        return data;
      };

      {
        // A key no earlier entry has, so that the pool builds a temporary with a level.
        const int numBefore = numBuilt;

        Scratch<EBAMRCellData> a = pool.acquire<EBAMRCellData>(phase::gas, 0, 4, 2, buildWithLevel);

        check(numBuilt == numBefore + 1, "a temporary with a level is built");

        // A reference does not count as a copy.
        const EBAMRCellData& view = a;

        check(view[0].refCount() == 1, "the pool holds the only reference to a level of a held temporary");
      }

      check(pool.allReturned(), "a temporary with levels is returned at scope exit");
    }

    // Clearing deletes the temporaries, so the next request builds again.
    pool.clear();

    {
      const int numBefore = numBuilt;

      Scratch<EBAMRCellData> a = pool.acquire<EBAMRCellData>(phase::gas, -1, 1, 2, buildCell);

      check(numBuilt == numBefore + 1, "a cleared pool builds new temporaries");
    }

    check(pool.allReturned(), "all temporaries are returned at the end");
  }

  pout() << "ScratchPool test passed" << endl;

  ChomboDischarge::finalize();
}
