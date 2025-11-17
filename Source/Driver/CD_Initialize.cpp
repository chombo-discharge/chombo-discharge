/* chombo-discharge
 * Copyright © 2025 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Initialize.cpp
  @brief  Implementation of CD_Initialize.cpp
  @author Robert Marskar
*/

// Std includes
#if defined(CH_MPI) || defined(CH_USE_PETSC)
#include <mpi.h>
#endif
#include <iostream>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <unistd.h>

// Chombo includes
#include <parstream.H>
#include <CH_Timer.H>
#include <SPMD.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_NamespaceHeader.H>

std::string dischargeInputFile = "<none>";
ParmParse*  dischargeParser    = nullptr;

#if defined(CH_USE_PETSC)
PetscErrorCode
#else
void
#endif
initialize(int argc, char* argv[])
{
#if defined(CH_USE_PETSC)
  const PetscErrorCode ierr = PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);
#else
#if defined(CH_MPI)
  MPI_Init(&argc, &argv);
#endif
#endif

  if (argc >= 2) {
    dischargeInputFile = argv[1];
    dischargeParser    = new ParmParse(argc - 2, argv + 2, nullptr, argv[1]);
  }
  else {
    dischargeParser = new ParmParse();
  }

  // Print initialization information from master rank
  bool printInitializationMessage = false;

  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != nullptr) {
    pout() << "\n";
    pout() << "================================================================================\n";
    pout() << "chombo-discharge initialization\n";
    pout() << "================================================================================\n";
    pout() << "  Working directory: " << cwd << "\n";
    pout() << "  Input file:        " << dischargeInputFile << "\n";
    pout() << "--------------------------------------------------------------------------------\n";
    pout() << "  MPI ranks:         " << numProc() << "\n";
#ifdef _OPENMP
    pout() << "  OpenMP threads:    " << omp_get_max_threads() << "\n";
#endif
#ifdef CH_USE_HDF5
    pout() << "  HDF5:              TRUE\n";
#else
    pout() << "  HDF5:              FALSE\n";
#endif
#ifdef CH_USE_PETSC
    pout() << "  PETSc:             TRUE\n";
#else
    pout() << "  PETSc:             FALSE\n";
#endif
    pout() << "================================================================================\n";
    pout() << std::endl;
  }

#if defined(CH_USE_PETSC)
  return ierr;
#endif
}

#if defined(CH_USE_PETSC)
PetscErrorCode
#else
int
#endif
finalize()
{
  delete dischargeParser;

  CH_TIMER_REPORT();

#if defined(CH_USE_PETSC)
  return PetscFinalize();
#else
#if defined(CH_MPI)
  MPI_Finalize();
#endif

  return 0;
#endif
}

#include <CD_NamespaceFooter.H>
