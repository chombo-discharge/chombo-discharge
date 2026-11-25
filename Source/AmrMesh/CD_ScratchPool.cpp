/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ScratchPool.cpp
 * @brief  Implementation of CD_ScratchPool.H
 * @author Robert Marskar
 */

// Std includes
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_ScratchPool.H>
#include <CD_NamespaceHeader.H>

ScratchPool::ScratchPool() noexcept
{}

ScratchPool::~ScratchPool() noexcept
{
  if (!(this->allReturned())) {
    std::ostringstream os;

    os << "ScratchPool::~ScratchPool -- destroying a pool while temporaries are in use:\n";

    this->report(os, "  ");

    MayDay::Abort(os.str().c_str());
  }
}

void
ScratchPool::clear() noexcept
{
  CH_TIME("ScratchPool::clear");

  ScratchPool::assertSerial("ScratchPool::clear");

  if (!(this->allReturned())) {
    std::ostringstream os;

    os << "ScratchPool::clear -- the grids change while temporaries are in use:\n";

    this->report(os, "  ");

    MayDay::Abort(os.str().c_str());
  }

  std::apply(
    [](auto&... a_list) -> void {
      (a_list.clear(), ...);
    },
    m_entries);
}

bool
ScratchPool::allReturned() const noexcept
{
  const int numInUse = this->numInUse<EBAMRCellData>() + this->numInUse<EBAMRFluxData>() +
                       this->numInUse<EBAMRIVData>() + this->numInUse<EBAMRFAB>() + this->numInUse<MFAMRCellData>() +
                       this->numInUse<MFAMRFluxData>() + this->numInUse<MFAMRIVData>() +
                       this->numInUse<LevelData<EBCellFAB>>();

  return numInUse == 0;
}

void
ScratchPool::report(std::ostream& a_os, const std::string& a_prefix) const noexcept
{
  this->reportType<EBAMRCellData>(a_os, a_prefix);
  this->reportType<EBAMRFluxData>(a_os, a_prefix);
  this->reportType<EBAMRIVData>(a_os, a_prefix);
  this->reportType<EBAMRFAB>(a_os, a_prefix);
  this->reportType<MFAMRCellData>(a_os, a_prefix);
  this->reportType<MFAMRFluxData>(a_os, a_prefix);
  this->reportType<MFAMRIVData>(a_os, a_prefix);
  this->reportType<LevelData<EBCellFAB>>(a_os, a_prefix);
}

void
ScratchPool::assertSerial(const char* a_caller) noexcept
{
#ifdef _OPENMP
  if (omp_in_parallel()) {
    const std::string err = std::string(a_caller) + " -- called inside an OpenMP parallel region";

    MayDay::Abort(err.c_str());
  }
#else
  (void)a_caller;
#endif
}

#include <CD_NamespaceFooter.H>
