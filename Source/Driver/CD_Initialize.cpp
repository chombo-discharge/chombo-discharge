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

std::string dischargeInputFile;
ParmParse   dischargeParser;

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
    dischargeParser.define(argc - 2, argv + 2, nullptr, argv[1]);
  }
  else {
    dischargeInputFile = "<No file provided>";
    dischargeParser.define(argc - 2, argv + 2, nullptr, nullptr);
  }

#ifdef _OPENMP
#pragma omp flush
#endif

  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != nullptr) {
    pout() << "\n";
    pout() << "==========================================================================================" << endl;
    pout() << R"(
      _                     _                     _ _          _                          
  ___| |__   ___  _ __ ___ | |__   ___         __| (_)___  ___| |__   __ _ _ __ __ _  ___ 
 / __| '_ \ / _ \| '_ ` _ \| '_ \ / _ \ _____ / _` | / __|/ __| '_ \ / _` | '__/ _` |/ _ \
| (__| | | | (_) | | | | | | |_) | (_) |_____| (_| | \__ \ (__| | | | (_| | | | (_| |  __/
 \___|_| |_|\___/|_| |_| |_|_.__/ \___/       \__,_|_|___/\___|_| |_|\__,_|_|  \__, |\___|
)" << endl;
    pout() << "  Working directory: " << cwd << "\n";
    pout() << "  Input file:        " << dischargeInputFile << "\n";
    pout() << "  Run command:       ";
    for (int i = 0; i < argc; i++) {
      pout() << argv[i];
      if (i < argc - 1) {
        pout() << " ";
      }
    }
    pout() << "\n";
    pout() << "--------------------------------------------------------------------------------\n";
#ifdef CH_MPI
    pout() << "  MPI:    TRUE (" << numProc() << " ranks)\n";
#else
    pout() << "  MPI:    FALSE\n";
#endif
#ifdef _OPENMP
    pout() << "  OpenMP: TRUE (" << omp_get_max_threads() << " threads)\n";
#else
    pout() << "  OpenMP: FALSE\n";
#endif
#ifdef CH_USE_HDF5
    pout() << "  HDF5:   TRUE\n";
#else
    pout() << "  HDF5:   FALSE\n";
#endif
#ifdef CH_USE_PETSC
    pout() << "  PETSc:  TRUE\n";
#else
    pout() << "  PETSc:  FALSE\n";
#endif
    pout() << "==========================================================================================" << endl;
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
